/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "amrCore.H"
#include "dimensionSets.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "sigFpe.H"
#include "emptyPolyPatch.H"
#include "cloudSupport.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(amrCore, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::amrCore::amrCore(fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    error_(nullptr),
    refiner_(nullptr),
    correctFluxes_(),
    expireSampledSurfacesOnLB_(false)
{
    initializeAMR(dict);
}


Foam::amrCore::amrCore(fvMesh& mesh)
:
    mesh_(mesh),
    error_(nullptr),
    refiner_(nullptr),
    correctFluxes_(),
    expireSampledSurfacesOnLB_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::amrCore::initializeAMR(const dictionary& dict)
{
    Info<< "amrCore: Initializing AMR" << endl;

    // Create error estimator
    error_ = errorEstimator::New(mesh_, dict);

    // Create refiner
    refiner_ = fvMeshRefiner::New(mesh_, dict, false, true);

    // Read flux correction settings
    readCorrectFluxes(dict);

    // Read sampledSurface expiration setting (opt-in, default false)
    expireSampledSurfacesOnLB_ =
        dict.getOrDefault("expireSampledSurfacesOnLB", false);
}


void Foam::amrCore::initializeLB(const dictionary& dict)
{
    if (!refiner_.valid())
    {
        Info<< "amrCore: Creating refiner for LB support" << endl;

        // Create a minimal dictionary for the refiner
        dictionary minimalDict;
        minimalDict.add("refiner", "hexRefiner");
        minimalDict.add("refine", false);
        minimalDict.add("unrefine", false);

        // Copy LB settings
        if (dict.found("balance"))
        {
            minimalDict.add("balance", dict.get<bool>("balance"));
        }
        if (dict.found("allowableImbalance"))
        {
            minimalDict.add
            (
                "allowableImbalance",
                dict.get<scalar>("allowableImbalance")
            );
        }
        if (dict.found("balanceInterval"))
        {
            minimalDict.add
            (
                "balanceInterval",
                dict.get<label>("balanceInterval")
            );
        }

        refiner_ = fvMeshRefiner::New(mesh_, minimalDict, false, true);
    }

    // Configure balancer within refiner
    refiner_->balancer().read(dict);

    Info<< "amrCore: Load balancing initialized" << endl;
}


void Foam::amrCore::readDict(const dictionary& dict)
{
    if (refiner_.valid())
    {
        refiner_->readDict(dict);
    }
    if (error_.valid())
    {
        error_->read(dict);
    }
    readCorrectFluxes(dict);

    // Re-read sampledSurface expiration setting (allows runtime changes)
    expireSampledSurfacesOnLB_ =
        dict.getOrDefault("expireSampledSurfacesOnLB", false);
}


void Foam::amrCore::readCorrectFluxes(const dictionary& dict)
{
    if (!dict.found("correctFluxes"))
    {
        return;
    }

    List<Pair<word>> fluxVelocities = List<Pair<word>>
    (
        dict.lookup("correctFluxes")
    );

    // Rework into hashtable
    correctFluxes_.resize(fluxVelocities.size());
    forAll(fluxVelocities, i)
    {
        correctFluxes_.insert(fluxVelocities[i][0], fluxVelocities[i][1]);
    }

    if (debug)
    {
        Info<< "amrCore: Flux corrections configured for "
            << correctFluxes_.size() << " flux field(s)" << endl;
    }
}


bool Foam::amrCore::refine()
{
    if (!amrInitialized())
    {
        FatalErrorInFunction
            << "AMR not initialized. Call initializeAMR() first."
            << exit(FatalError);
    }

    // Check timing constraints
    if (!refiner_->canRefine(true) && !refiner_->canUnrefine(true))
    {
        return false;
    }

    // Re-read dictionary for on-the-fly modifications
    // (handled by caller if needed)

    // Store cloud positions before topology change
    cloudSupport::storeGlobalPositions(mesh_);

    // Update error field
    error_->update();
    error_->error().correctBoundaryConditions();

    // Protect patches from refinement
    label nProtected = error_->protectPatches();
    Info<< "Protecting " << returnReduce(nProtected, sumOp<label>())
        << " cells next to requested boundary patches." << endl;

    // Perform refinement
    bool changed = refiner_->refine(error_->error(), error_->maxRefinement());
    reduce(changed, orOp<bool>());

    return changed;
}


bool Foam::amrCore::balance()
{
    if (!lbInitialized())
    {
        FatalErrorInFunction
            << "Load balancing not initialized. Call initializeLB() first."
            << exit(FatalError);
    }

    // Perform balance (includes canBalance check internally)
    // Note: Do NOT call canBalance(true) here separately as that would
    // double-increment the balance iteration counter and cause the second
    // call inside refiner_->balance() to fail the interval check.
    // The storeGlobalPositions is called inside fvMeshBalance::distribute().
    bool changed = refiner_->balance();
    reduce(changed, orOp<bool>());

    return changed;
}


void Foam::amrCore::correctFluxes(const mapPolyMesh& map)
{
    const labelList& faceMap = map.faceMap();
    const labelList& reverseFaceMap = map.reverseFaceMap();

    // Storage for any master faces. These will be the original faces
    // on the coarse cell that get split into four (or rather the
    // master face gets modified and three faces get added from the master)
    labelHashSet masterFaces;

    forAll(faceMap, facei)
    {
        label oldFacei = faceMap[facei];

        if (oldFacei >= 0)
        {
            label masterFacei = reverseFaceMap[oldFacei];

            if (masterFacei < 0)
            {
                FatalErrorInFunction
                    << "Problem: should not have removed faces"
                    << " when refining."
                    << nl << "face:" << facei << abort(FatalError);
            }
            else if (masterFacei != facei)
            {
                masterFaces.insert(masterFacei);
            }
        }
    }

    if (debug)
    {
        Pout<< "Found " << masterFaces.size() << " split faces " << endl;
    }

    // Check if it's a flux field through dims
    auto isFlux = [&](const surfaceScalarField& df)
    {
        return
            df.dimensions() == dimArea*dimVelocity
         || df.dimensions() == dimArea*dimVelocity*dimDensity;
    };

    HashTable<surfaceScalarField*> fluxes
    (
        mesh_.lookupClass<surfaceScalarField>()
    );

    forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
    {
        if (!isFlux(*iter()))
        {
            continue;
        }
        if (!correctFluxes_.found(iter.key()))
        {
            continue;
        }

        const word& UName = correctFluxes_[iter.key()];

        if (UName == "none")
        {
            continue;
        }

        if (UName == "NaN")
        {
            Pout<< "Setting surfaceScalarField " << iter.key()
                << " to NaN" << endl;

            surfaceScalarField& phi = *iter();
            sigFpe::fillNan(phi.primitiveFieldRef());
            continue;
        }

        if (debug)
        {
            Pout<< "Mapping flux " << iter.key()
                << " using interpolated flux " << UName
                << endl;
        }

        surfaceScalarField& phi = *iter();
        const surfaceScalarField phiU
        (
            fvc::interpolate
            (
                mesh_.lookupObject<volVectorField>(UName)
            )
          & mesh_.Sf()
        );

        // Recalculate new internal faces
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label oldFacei = faceMap[facei];

            if (oldFacei == -1)
            {
                // Inflated/appended
                phi[facei] = phiU[facei];
            }
            else if (reverseFaceMap[oldFacei] != facei)
            {
                // face-from-masterface
                phi[facei] = phiU[facei];
            }
        }

        // Recalculate new boundary faces
        surfaceScalarField::Boundary& phiBf = phi.boundaryFieldRef();
        forAll(phiBf, patchi)
        {
            fvsPatchScalarField& patchPhi = phiBf[patchi];
            const fvsPatchScalarField& patchPhiU =
                phiU.boundaryField()[patchi];

            label facei = patchPhi.patch().start();

            forAll(patchPhi, i)
            {
                label oldFacei = faceMap[facei];

                if (oldFacei == -1)
                {
                    // Inflated/appended
                    patchPhi[i] = patchPhiU[i];
                }
                else if (reverseFaceMap[oldFacei] != facei)
                {
                    // face-from-masterface
                    patchPhi[i] = patchPhiU[i];
                }

                facei++;
            }
        }

        // Update master faces
        forAllConstIter(labelHashSet, masterFaces, iter)
        {
            label facei = iter.key();

            if (mesh_.isInternalFace(facei))
            {
                phi[facei] = phiU[facei];
            }
            else
            {
                label patchi = mesh_.boundaryMesh().whichPatch(facei);

                if (!isA<emptyPolyPatch>(mesh_.boundaryMesh()[patchi]))
                {
                    label i = facei - mesh_.boundaryMesh()[patchi].start();

                    const fvsPatchScalarField& patchPhiU =
                        phiU.boundaryField()[patchi];

                    fvsPatchScalarField& patchPhi = phiBf[patchi];

                    patchPhi[i] = patchPhiU[i];
                }
            }
        }
    }
}


void Foam::amrCore::updateMesh(const mapPolyMesh& map)
{
    // Correct fluxes on new faces
    correctFluxes(map);

    // Update refiner
    if (refiner_.valid())
    {
        refiner_->updateMesh(map);
    }

    // Note: Cloud remapping is NOT done here because this is called
    // from adaptiveFvMesh::updateMesh BEFORE fvMesh::updateMesh.
    // Calling cloud.autoMap would trigger mesh_.V() which reconstructs
    // volumes with the NEW mesh size, causing the check in
    // fvMesh::updateMesh to fail.
    //
    // Cloud remapping is done in adaptiveFvMesh::updateMesh AFTER
    // fvMesh::updateMesh completes.
}


void Foam::amrCore::distribute(const mapDistributePolyMesh& map)
{
    if (refiner_.valid())
    {
        refiner_->distribute(map);
    }
}


// ************************************************************************* //
