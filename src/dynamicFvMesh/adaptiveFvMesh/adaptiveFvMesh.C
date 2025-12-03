/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020 Synthetik Applied Technologies: |    Modified original
                            dynamicRefineBalanceBlastFvMesh class
                            to be more appilcable to compressible flows.
                            Improved compatibility with snappyHexMesh.
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

#include "adaptiveFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "cloudSupport.H"
#include "sampledSurfaceWorkaround.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adaptiveFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        adaptiveFvMesh,
        IOobject
    );
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::adaptiveFvMesh::dynamicMeshDictIOobject
(
    const IOobject& io
)
{
    IOobject dictHeader
    (
        "dynamicMeshDict",
        io.time().constant(),
        (io.name() == polyMesh::defaultRegion ? "" : io.name()),
        io.db(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    // defaultRegion (region0) gets loaded from constant, other ones get loaded
    // from constant/<regionname>. Normally we'd use polyMesh::dbDir() but we
    // haven't got a polyMesh yet ...
    return IOobject
    (
        "dynamicMeshDict",
        io.time().constant(),
        (io.name() == polyMesh::defaultRegion ? "" : io.name()),
        io.db(),
        (
            dictHeader.typeHeaderOk<IOdictionary>(true)
          ? IOobject::MUST_READ_IF_MODIFIED
          : IOobject::NO_READ
        ),
        IOobject::NO_WRITE
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::adaptiveFvMesh::readDict()
{
    const dictionary refineDict
    (
        dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );

    // Delegate to amrCore
    amrCore_.readDict(refineDict);
}


void Foam::adaptiveFvMesh::updateMesh(const mapPolyMesh& map)
{
    // Delegate flux correction and refiner update to amrCore
    amrCore_.updateMesh(map);

    // Call parent updateMesh
    fvMesh::updateMesh(map);

    // Remap clouds AFTER fvMesh::updateMesh completes.
    // This must happen after fvMesh::updateMesh because cloud.autoMap
    // triggers mesh_.V() which reconstructs volumes. If done before,
    // the volume size check in fvMesh::updateMesh would fail.
    cloudSupport::autoMapClouds(*this, map);

    // Expire sampled surfaces in surfaceFieldValue function objects.
    // OpenFOAM's surfaceFieldValue::updateMesh() doesn't properly expire
    // its internal sampledPtr_ cache, leading to stale face indices after
    // mesh topology changes (refinement/unrefinement).
    expireSampledSurfaces(time(), amrCore_.expireSampledSurfacesOnLB());
}


void Foam::adaptiveFvMesh::distribute
(
    const mapDistributePolyMesh& map
)
{
    // Delegate to amrCore
    amrCore_.distribute(map);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adaptiveFvMesh::adaptiveFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshDict_(dynamicMeshDictIOobject(io)),
    amrCore_(*this),
    currentTimeIndex_(-1)
{
    // Read dictionary and initialize amrCore
    const dictionary refineDict
    (
        dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );

    amrCore_.initializeAMR(refineDict);

    // Initialize load balancing if enabled
    if (refineDict.getOrDefault("balance", false))
    {
        amrCore_.initializeLB(refineDict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adaptiveFvMesh::~adaptiveFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adaptiveFvMesh::mapFields(const mapPolyMesh& mpm)
{
    dynamicFvMesh::mapFields(mpm);

    // Correct surface fields on introduced internal faces. These get
    // created out-of-nothing so get an interpolated value.
    mapNewInternalFaces<scalar>(mpm.faceMap());
    mapNewInternalFaces<vector>(mpm.faceMap());
    mapNewInternalFaces<sphericalTensor>(mpm.faceMap());
    mapNewInternalFaces<symmTensor>(mpm.faceMap());
    mapNewInternalFaces<tensor>(mpm.faceMap());
}

bool Foam::adaptiveFvMesh::firstUpdate()
{
    if (currentTimeIndex_ >= time().timeIndex())
    {
        return false;
    }
    currentTimeIndex_ = time().timeIndex();
    return true;
}

bool Foam::adaptiveFvMesh::update()
{
    if (!firstUpdate()) return false;

    bool changed =
        (amrCore_.refiner().canRefine(true) || amrCore_.refiner().canUnrefine(true))
     && refine();

    // Check for load balancing independently of refinement
    // This allows balanceInterval to differ from refineInterval
    if (amrCore_.lbInitialized() && Pstream::parRun())
    {
        if (amrCore_.balance())
        {
            changed = true;
        }
    }

    reduce(changed, orOp<bool>());

    return changed;
}


bool Foam::adaptiveFvMesh::refine()
{
    // Re-read dictionary for on-the-fly modifications
    readDict();

    // Delegate to amrCore
    return amrCore_.refine();
}


bool Foam::adaptiveFvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    return dynamicFvMesh::writeObject(streamOpt, valid);
}


// ************************************************************************* //
