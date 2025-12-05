/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "fvMeshBalance.H"
#include "decompositionMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "RefineBalanceMeshObject.H"
#include "cloudSupport.H"
#include "preserveFaceZonesConstraint.H"
#include "singleProcessorFaceSetsConstraint.H"
#include "preservePatchesConstraint.H"
#include "preserveBafflesConstraint.H"
#include "sampledSurfaceWorkaround.H"
#include "dynamicMotionSolverFvMesh.H"
#include "pointIOField.H"

using namespace Foam::decompositionConstraints;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshBalance, 0);

bool fvMeshBalance::balancing = false;
}

bool Foam::fvMeshBalance::isBalancing()
{
    return balancing;
}

bool Foam::fvMeshBalance::isBalancing(const polyMesh& mesh)
{
    if (balancing)
    {
        if (!mesh.db().foundObject<polyMesh>(mesh.name()))
        {
            return false;
        }
        if (&mesh == &mesh.db().lookupObject<polyMesh>(mesh.name()))
        {
            return !Pstream::parRun();
        }
    }
    return false;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshBalance::fvMeshBalance(fvMesh& mesh)
:
    mesh_(mesh),
    decompositionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                mesh.time().system(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),
    modified_(false),
    constraintsDict_(decompositionDict_.subDictPtr("constraints")),
    //preserveFaceZonesDict_(nullptr),
    //singleProcessorFaceSetsDict_(nullptr),
    //preservePatchesDict_(nullptr),
    //preserveBafflesDict_(nullptr),
    distributor_(mesh_),
    loadPolicy_(nullptr),  // Will be set by read() with proper dict
    expireSampledSurfacesOnLB_(false)
{
    if (!constraintsDict_)
    {
        decompositionDict_.add("constraints", dictionary());
        constraintsDict_ = decompositionDict_.subDictPtr("constraints");
    }
    else
    {
        wordList toc(decompositionDict_.toc());
        forAll(toc, i)
        {
            if (decompositionDict_.isDict(toc[i]))
            {
                dictionary& dict(decompositionDict_.subDict(toc[i]));
                word type(dict.lookupOrDefault<word>("type", "none"));

                //if (type == preserveFaceZonesConstraint::typeName)
                //{
                //    preserveFaceZonesDict_ = &dict;
                //}
                //else if
                //(
                //    type == singleProcessorFaceSetsConstraint::typeName
                //)
                //{
                //    singleProcessorFaceSetsDict_ = &dict;
                //}
                //else if (type == preservePatchesConstraint::typeName)
                //{
                //    preservePatchesDict_ = &dict;
                //}
                //else if (type == preserveBafflesConstraint::typeName)
                //{
                //    preserveBafflesDict_ = &dict;
                //}
            }
        }
    }
}


Foam::fvMeshBalance::fvMeshBalance
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    decompositionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                mesh.time().system(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),
    modified_(false),
    constraintsDict_(decompositionDict_.subDictPtr("constraints")),
    //preserveFaceZonesDict_(nullptr),
    //singleProcessorFaceSetsDict_(nullptr),
    //preservePatchesDict_(nullptr),
    //preserveBafflesDict_(nullptr),
    distributor_(mesh_),
    loadPolicy_(dict.lookupOrDefault<Switch>("balance", false) ? loadPolicy::New(mesh, dict) : nullptr),
    expireSampledSurfacesOnLB_(dict.getOrDefault("expireSampledSurfacesOnLB", false))
{
    if (!constraintsDict_)
    {
        decompositionDict_.add("constraints", dictionary());
        constraintsDict_ = decompositionDict_.subDictPtr("constraints");
    }
    else
    {
        wordList toc(decompositionDict_.toc());
        forAll(toc, i)
        {
            if (decompositionDict_.isDict(toc[i]))
            {
                dictionary& dict(decompositionDict_.subDict(toc[i]));
                word type(dict.lookupOrDefault<word>("type", "none"));

                //if (type == preserveFaceZonesConstraint::typeName)
                //{
                //    preserveFaceZonesDict_ = &dict;
                //}
                //else if
                //(
                //    type == singleProcessorFaceSetsConstraint::typeName
                //)
                //{
                //    singleProcessorFaceSetsDict_ = &dict;
                //}
                //else if (type == preservePatchesConstraint::typeName)
                //{
                //    preservePatchesDict_ = &dict;
                //}
                //else if (type == preserveBafflesConstraint::typeName)
                //{
                //    preserveBafflesDict_ = &dict;
                //}
            }
        }
    }
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshBalance::~fvMeshBalance()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshBalance::read(const dictionary& balanceDict)
{
    if (!Pstream::parRun())
    {
        loadPolicy_ = nullptr;
        return;
    }

    // Only create/update loadPolicy if:
    // 1. We don't have one yet, OR
    // 2. The dict explicitly specifies a loadPolicy type
    // This prevents overwriting a configured loadPolicy with the default
    // when re-reading from a dict that doesn't have loadPolicy entry
    if (!loadPolicy_ || balanceDict.found("loadPolicy"))
    {
        loadPolicy_ = balanceDict.lookupOrDefault("balance", true)
            ? loadPolicy::New(mesh_, balanceDict)
            : nullptr;
    }

    // Read sampledSurface expiration setting (opt-in, allows runtime changes)
    expireSampledSurfacesOnLB_ =
        balanceDict.getOrDefault("expireSampledSurfacesOnLB", false);

    if (!loadPolicy_)
    {
        return;
    }

    // Change decomposition method if entry is present
    if (balanceDict.found("method"))// || balanceDict.found("decomposer"))
    {
        word method(balanceDict.lookup("method"));
    //        (
    //            {"method"}
    //             {"decomposer", "method"},
    //        );
        decompositionDict_.set("method", method);

        // Remove optional coeffs dictionary since it would override entries and
        // not necessarily be overridden
        if (decompositionDict_.isDict(method + "Coeffs"))
        {
            decompositionDict_.remove(method + "Coeffs");
        }
    }
    decompositionDict_ <<= balanceDict;
}


void Foam::fvMeshBalance::addConstraint(const word& dictName, const dictionary& dict)
{
    // Add constraints dictionary
    if (!constraintsDict_->found(dictName))
    {
        modified_ = true;

        constraintsDict_->set(dictName, dict);

        // We need to apply the updated constraints, and this can only
        // be done when the decomposer is created
        if (decomposer_.valid())
        {
            decomposer_.clear();
        }
    }
}


//void Foam::fvMeshBalance::preserveFaceZone(const wordRe& zoneName)
//{
//    if (!preserveFaceZonesDict_)
//    {
//        constraintsDict_->add("faceZones", dictionary());
//        preserveFaceZonesDict_ = constraintsDict_->subDictPtr("faceZones");
//        preserveFaceZonesDict_->set
//        (
//            "type",
//            preserveFaceZonesConstraint::typeName
//        );
//    }
//    wordReList zones
//    (
//        preserveFaceZonesDict_->lookupOrDefault("zones", wordReList())
//    );
//
//    bool exists = false;
//    forAll(zones, i)
//    {
//        if (zones[i].match(zoneName))
//        {
//            exists = true;
//            break;
//        }
//    }
//    if (!exists)
//    {
//        modified_ = true;
//        zones.append(zoneName);
//        preserveFaceZonesDict_->set("zones", zones);
//
//        if (decomposer_.valid())
//        {
//            decomposer_.clear();
//        }
//    }
//}


//void Foam::fvMeshBalance::singleProcessorFaceSet
//(
//    const word& setName,
//    const label proc
//)
//{
//    //if (!preserveFaceZonesDict_)
//    //{
//    //    constraintsDict_->add("singleProcessorFaceSets", dictionary());
//    //    singleProcessorFaceSetsDict_ =
//    //        constraintsDict_->subDictPtr("singleProcessorFaceSets");
//    //    singleProcessorFaceSetsDict_->set
//    //    (
//    //        "type",
//    //        singleProcessorFaceSetsConstraint::typeName
//    //    );
//    //}
//    List<Tuple2<word, label>> setNameAndProcs
//    (
//        singleProcessorFaceSetsDict_->lookupOrDefault
//        (
//            "singleProcessorFaceSets",
//            List<Tuple2<word, label>>()
//        )
//    );
//
//    bool exists = false;
//    forAll(setNameAndProcs, i)
//    {
//        if (setName == setNameAndProcs[i].first())
//        {
//            exists = true;
//            break;
//        }
//    }
//    if (!exists)
//    {
//        modified_ = true;
//        setNameAndProcs.append(Tuple2<word, label>(setName, proc));
//        singleProcessorFaceSetsDict_->set
//        (
//            "setNameAndProcs",
//            setNameAndProcs
//        );
//
//        if (decomposer_.valid())
//        {
//            decomposer_.clear();
//        }
//    }
//}


//void Foam::fvMeshBalance::preservePatch(const wordRe& patchName)
//{
//    if (!preservePatchesDict_)
//    {
//        constraintsDict_->add("patches", dictionary());
//        preservePatchesDict_ = constraintsDict_->subDictPtr("patches");
//        preservePatchesDict_->set
//        (
//            "type",
//            preservePatchesConstraint::typeName
//        );
//    }
//    wordReList patches
//    (
//        preservePatchesDict_->lookupOrDefault("patches", wordReList())
//    );
//
//    bool exists = false;
//    forAll(patches, i)
//    {
//        if (patches[i].match(patchName))
//        {
//            exists = true;
//            break;
//        }
//    }
//    if (!exists)
//    {
//        modified_ = true;
//        patches.append(patchName);
//        preservePatchesDict_->set("patches", patches);
//
//        if (decomposer_.valid())
//        {
//            decomposer_.clear();
//        }
//    }
//}


//void Foam::fvMeshBalance::preserveBaffles()
//{
//    if (!preserveBafflesDict_)
//    {
//        modified_ = true;
//        constraintsDict_->add("baffles", dictionary());
//        preserveBafflesDict_ = constraintsDict_->subDictPtr("baffles");
//        preserveBafflesDict_->set
//        (
//            "type",
//            preserveBafflesConstraint::typeName
//        );
//    }
//}


void Foam::fvMeshBalance::makeDecomposer() const
{
    if (decomposer_.valid())
    {
        FatalErrorInFunction
            << "Decomposer already set" << endl
            << abort(FatalError);
    }
    decomposer_ = decompositionMethod::New(decompositionDict_);

    returnReduce(1, maxOp<label>());
    if (!decomposer_->parallelAware())
    {
        WarningInFunction
            << "You have selected decomposition method "
            << decomposer_->typeName
            << " which is not parallel aware." << endl;
    }
}


Foam::decompositionMethod& Foam::fvMeshBalance::decomposer() const
{
    if (!decomposer_.valid())
    {
        makeDecomposer();
    }
    return decomposer_();
}


bool Foam::fvMeshBalance::canBalance() const
{
    if (!loadPolicy_)
    {
        return false;
    }

    if(!loadPolicy_->canBalance()) return false;

    // Decompose the mesh with uniform weights
    // The refinementHistory constraint is applied internally
    distribution_ = decomposer().decompose
    (
        mesh_,
        loadPolicy_->cellWeights()
        //scalarField(mesh_.nCells(), 1.0)
    );

    // Check if distribution will improve anything
    return loadPolicy_->willBeBeneficial(distribution_);

    //labelList procLoadNew(Pstream::nProcs(), 0);
    //forAll(distribution_, celli)
    //{
    //    procLoadNew[distribution_[celli]]++;
    //}
    //reduce(procLoadNew, sumOp<labelList>());
    //if (min(procLoadNew) == 0)
    //{
    //    DebugInfo
    //        << "New distribtion results in a load of 0. Skipping" << endl;
    //    return false;
    //}
    //scalar averageLoadNew
    //(
    //    scalar(sum(procLoadNew))/scalar(Pstream::nProcs())
    //);
    //scalar maxDevNew(max(mag(procLoadNew - averageLoadNew))/averageLoadNew);

    //if (maxDevNew > maxImbalanceRatio*0.99)
    //{
    //    Info
    //        << "    Not balancing because the new distribution does" << nl
    //        << "    not improve the load. Skipping" << nl
    //        << "    old imbalance: " << maxImbalanceRatio << nl
    //        << "    new imbalance: " << maxDevNew << nl
    //        << endl;
    //    return false;
    //}

    //return true;
}


Foam::autoPtr<Foam::mapDistributePolyMesh>
Foam::fvMeshBalance::distribute()
{
    // Synchronize oldTime fields across processors before distribution.
    // Different processors may end up with different fields after mesh changes
    // (e.g., U_0 from particle-wall interactions in kinematicCloud).
    // Hence forcing oldTime on any field that has "oldTime" on any processor
    syncOldTimeFields<volScalarField>();
    syncOldTimeFields<volVectorField>();
    syncOldTimeFields<volSphericalTensorField>();
    syncOldTimeFields<volSymmTensorField>();
    syncOldTimeFields<volTensorField>();
    syncOldTimeFields<surfaceScalarField>();
    syncOldTimeFields<surfaceVectorField>();
    syncOldTimeFields<surfaceSphericalTensorField>();
    syncOldTimeFields<surfaceSymmTensorField>();
    syncOldTimeFields<surfaceTensorField>();

    blastMeshObject::preDistribute<fvMesh>(mesh_);

    // Check if mesh uses a motion solver - special handling is required
    auto* motionMeshPtr = dynamic_cast<dynamicMotionSolverFvMesh*>(&mesh_);
    bool oldMoving = false;

    if (motionMeshPtr)
    {
        // Temporarily disable mesh motion to prevent oldPointsPtr_ recreation.
        // The oldPoints() getter creates oldPointsPtr_ on-demand when moving_ is true.
        // If any code below (cloud operations, distribute) calls oldPoints() while
        // moving_ is true, it would recreate oldPointsPtr_ after we clear it.
        // This follows the same pattern used in fvMeshDistribute::distribute().
        oldMoving = mesh_.moving(false);

        // Clear motion data BEFORE distribution.
        // Motion solvers store oldPoints and oldCellCentres for mesh motion.
        // During fvMeshAdder::add() -> polyMesh::updateMesh(), these are mapped
        // using pointMap(), but this can fail with invalid data during redistribution.
        // resetMotion() clears these pointers, preventing the mapping crash.
        // The moving(false) above prevents oldPoints() from recreating oldPointsPtr_.
        mesh_.resetMotion();
    }

    // Store global positions for all clouds before distribution
    // -- this is unified externally for all cloud types
    cloudSupport::storeGlobalPositions(mesh_);

    // Distribute clouds to new processors BEFORE mesh distribution
    // This transfers particles based on which processor their cell is going to
    cloudSupport::distributeClouds(mesh_, distribution_);

    Info<< "Distributing the mesh ..." << endl;
    balancing = true;
    autoPtr<mapDistributePolyMesh> map =
        distributor_.distribute(distribution_);
    balancing = false;

    if (motionMeshPtr)
    {
        // Restore mesh motion state after distribution
        mesh_.moving(oldMoving);

        // Reinitialize motion solver after redistribution.
        // The motion solver's points0_ field has the wrong size after redistribution
        // because it's a pointIOField (not a pointVectorField) and is not automatically
        // mapped during mesh redistribution. Calling init(false) recreates the motion
        // solver with the correct mesh topology without reinitializing the base mesh.
        // Before reinitializing, we must write the current mesh points as "points0"
        // because the motion solver constructor reads from files, and the old files
        // have the pre-redistribution point count.
        DebugInfo << "Reinitializing motion solver after redistribution" << endl;

        // Write current mesh points as "points0" so the motion solver
        // constructor can read the correct redistributed points.
        // Use the current time instance so findInstance() finds it.
        pointIOField points0
        (
            IOobject
            (
                "points0",
                mesh_.time().timeName(),
                polyMesh::meshSubDir,
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_.points()
        );
        points0.write();

        // Now reinitialize the motion solver - it will read the correct points
        motionMeshPtr->init(false);

        // Clean up the temporary points0 file to prevent interference
        // with subsequent mesh operations (e.g., refinement)
        Foam::rm(points0.objectPath());
    }

    Info << "Successfully distributed mesh" << endl;

    //label procLoadNew(mesh_.nCells());
    //label overallLoadNew(returnReduce(procLoadNew, sumOp<label>()));
    //scalar averageLoadNew(overallLoadNew/scalar(Pstream::nProcs()));
    //scalar maxDevNew
    //(
    //    returnReduce(mag(procLoadNew - averageLoadNew), maxOp<scalar>())
    //);
    //Info << "New max imbalance: " << maxDevNew/averageLoadNew*100.0 << "%"
    //    << endl;
    //if (debug)
    //{
    //    Pout<< " localImbalance = "
    //        << mag(procLoadNew - averageLoadNew)*100.0/averageLoadNew << "%, "
    //        << "Cells = " << procLoadNew
    //         << endl;
    //}

    blastMeshObject::distribute<fvMesh>(mesh_, map());

    // Expire sampled surfaces in surfaceFieldValue function objects.
    // OpenFOAM's surfaceFieldValue::updateMesh() doesn't properly expire
    // its internal sampledPtr_ cache, leading to stale face indices after
    // mesh redistribution. This workaround accesses the private member
    // directly to force expiration.
    expireSampledSurfaces(mesh_.time(), expireSampledSurfacesOnLB_);

    // Relocate particles to their new cells after mesh distribution
    cloudSupport::relocateClouds(mesh_);

    //Correct values on all coupled patches
    correctBoundaries<volScalarField>();
    correctBoundaries<volVectorField>();
    correctBoundaries<volSphericalTensorField>();
    correctBoundaries<volSymmTensorField>();
    correctBoundaries<volTensorField>();

    correctBoundaries<pointScalarField>();
    correctBoundaries<pointVectorField>();
    correctBoundaries<pointSphericalTensorField>();
    correctBoundaries<pointSymmTensorField>();
    correctBoundaries<pointTensorField>();

    return map;
}

bool Foam::fvMeshBalance::write(const bool write) const
{
    if
    (
        loadPolicy_ && modified_ && write &&
        decompositionDict_.lookupOrDefault("writeDecomposeDict", false)
    )
    {
        modified_ = false;
        IOdictionary decomposeParDict
        (
            IOobject
            (
                "decomposeParDict",
                mesh_.time().caseSystem(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            decompositionDict_
        );
        return decomposeParDict.regIOobject::write();
    }
    return true;
}


// ************************************************************************* //
