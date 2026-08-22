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

#include "adaptiveOversetFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "cloudSupport.H"
#include "sampledSurfaceWorkaround.H"
#include "dynMeshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adaptiveOversetFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        adaptiveOversetFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::adaptiveOversetFvMesh::dynamicMeshDictIOobject
(
    const IOobject& io
)
{
    // Not registered: dynamicOversetFvMesh reads the same dictionary
    return IOobject
    (
        "dynamicMeshDict",
        io.time().constant(),
        (io.name() == polyMesh::defaultRegion ? "" : io.name()),
        io.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::adaptiveOversetFvMesh::readDict()
{
    amrCore_.readDict
    (
        dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );
}


void Foam::adaptiveOversetFvMesh::updateMesh(const mapPolyMesh& mpm)
{
    // zoneID has to be right before anything rebuilds the stencil, and
    // correcting the overset boundary conditions inside fvMesh::updateMesh
    // does exactly that
    if (oversetHandler* handler = amrCore_.oversetHandlerPtr())
    {
        handler->updateMesh(mpm);
    }

    // Motion solvers, mesh objects (including the overset stencil) and all
    // registered fields. Runs before amrCore so that flux correction sees
    // fields already resized to the new topology
    dynamicOversetFvMesh::updateMesh(mpm);

    // Flux correction and refiner update
    amrCore_.updateMesh(mpm);

    // Clouds have to be remapped after fvMesh::updateMesh, see adaptiveFvMesh
    cloudSupport::autoMapClouds(*this, mpm);

    expireSampledSurfaces(time(), amrCore_.expireSampledSurfacesOnLB());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adaptiveOversetFvMesh::adaptiveOversetFvMesh(const IOobject& io)
:
    dynamicOversetFvMesh(io),
    dynamicMeshDict_(dynamicMeshDictIOobject(io)),
    amrCore_(*this),
    currentTimeIndex_(-1)
{
    const dictionary refineDict
    (
        dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );

    amrCore_.initializeAMR(refineDict);

    if (refineDict.getOrDefault("balance", false))
    {
        amrCore_.initializeLB(refineDict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adaptiveOversetFvMesh::mapFields(const mapPolyMesh& mpm)
{
    dynamicOversetFvMesh::mapFields(mpm);

    // Correct surface fields on introduced internal faces. These get
    // created out-of-nothing so get an interpolated value.
    meshTools::mapNewInternalFaces<scalar>(*this, mpm.faceMap());
    meshTools::mapNewInternalFaces<vector>(*this, mpm.faceMap());
    meshTools::mapNewInternalFaces<sphericalTensor>(*this, mpm.faceMap());
    meshTools::mapNewInternalFaces<symmTensor>(*this, mpm.faceMap());
    meshTools::mapNewInternalFaces<tensor>(*this, mpm.faceMap());
}


bool Foam::adaptiveOversetFvMesh::refine()
{
    readDict();

    return amrCore_.refine();
}


bool Foam::adaptiveOversetFvMesh::update()
{
    bool adapted = false;

    if (currentTimeIndex_ < time().timeIndex())
    {
        currentTimeIndex_ = time().timeIndex();

        adapted =
            (
                amrCore_.refiner().canRefine(true)
             || amrCore_.refiner().canUnrefine(true)
            )
         && refine();

        if (amrCore_.lbInitialized() && Pstream::parRun() && amrCore_.balance())
        {
            adapted = true;
        }

        reduce(adapted, orOp<bool>());
    }

    // Mesh motion runs after adaptation, not before: balancing resets the
    // motion state, and only a subsequent movePoints re-establishes the mesh
    // flux the solver asks for. This also refreshes the overset addressing
    const bool moved = dynamicOversetFvMesh::update();

    if (adapted && !moved)
    {
        // Topology changed but nothing moved, so rebuild the extended
        // addressing and re-interpolate before the solver runs
        oversetFvMeshBase::update();
    }

    return moved || adapted;
}


bool Foam::adaptiveOversetFvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    // zoneID matching the adapted mesh, so that restarts do not fall back
    // to the zoneID of the original mesh
    if (amrCore_.oversetHandlerPtr())
    {
        amrCore_.oversetHandlerPtr()->writeZoneID();
    }

    return dynamicOversetFvMesh::writeObject(streamOpt, valid);
}


// ************************************************************************* //
