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

#include "oversetHandler.H"
#include "cellCellStencilObject.H"
#include "cellZoneMesh.H"
#include "mapDistributePolyMesh.H"
#include "mapPolyMesh.H"
#include "oversetPolyPatch.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetHandler, 0);
    defineRunTimeSelectionTable(oversetHandler, mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetHandler::oversetHandler(fvMesh& mesh)
:
    regIOobject
    (
        IOobject
        (
            oversetHandler::typeName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    oversetMesh_(mesh),
    oldZoneID_(),
    protectZones_()
{
    // Pull zoneID onto the registry while the mesh still matches the zoneID
    // field on disk. From here on it is only ever remapped
    zoneID();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

bool Foam::oversetHandler::isOverset(const fvMesh& mesh)
{
    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        if (isA<oversetPolyPatch>(pp))
        {
            return true;
        }
    }

    return false;
}


Foam::oversetHandler* Foam::oversetHandler::lookup(const fvMesh& mesh)
{
    return mesh.getObjectPtr<oversetHandler>(oversetHandler::typeName);
}


Foam::oversetHandler* Foam::oversetHandler::New
(
    fvMesh& mesh,
    const dictionary& dict
)
{
    const word handlerType
    (
        dict.getOrDefault<word>("oversetHandler", oversetHandler::typeName)
    );

    if (handlerType == "none" || !isOverset(mesh))
    {
        return nullptr;
    }

    if (oversetHandler* existing = lookup(mesh))
    {
        return existing;
    }

    autoPtr<oversetHandler> handlerPtr;

    if (handlerType == oversetHandler::typeName)
    {
        handlerPtr.reset(new oversetHandler(mesh));
    }
    else
    {
        const auto* tablePtr = meshConstructorTablePtr_;

        if (!tablePtr || !tablePtr->found(handlerType))
        {
            FatalIOErrorInFunction(dict)
                << "Unknown oversetHandler type " << handlerType << nl << nl
                << "Valid oversetHandler types are :" << nl
                << (tablePtr ? tablePtr->sortedToc() : wordList())
                << exit(FatalIOError);
        }

        handlerPtr.reset(tablePtr->cfind(handlerType).val()(mesh).ptr());
    }

    handlerPtr->read(dict);

    Info<< "oversetHandler: registering " << handlerPtr->type()
        << " for overset mesh " << mesh.name() << endl;

    oversetHandler* handler = handlerPtr.ptr();
    handler->store();

    return handler;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::labelIOList& Foam::oversetHandler::zoneID()
{
    // Creates and registers the list on first call
    cellCellStencil::zoneID(oversetMesh_);

    return *oversetMesh_.getObjectPtr<labelIOList>("zoneID");
}


void Foam::oversetHandler::expireStencil()
{
    auto* stencilPtr =
        oversetMesh_.getObjectPtr<cellCellStencilObject>
        (
            cellCellStencilObject::typeName
        );

    if (stencilPtr)
    {
        stencilPtr->checkOut();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::oversetHandler::read(const dictionary& dict)
{
    protectZones_ = dict.getOrDefault<wordRes>("protectZones", wordRes());
}


Foam::label Foam::oversetHandler::protectCells(volScalarField& error) const
{
    label nProtected = 0;

    // Cells carried by a motion solver cannot be refined: OpenFOAM places
    // the points introduced by refinement in the reference configuration by
    // assuming the motion is a pure scaling (points0MotionSolver::updateMesh),
    // which is wrong under rotation and tangles the mesh as the body turns
    for (const cellZone& cz : oversetMesh_.cellZones())
    {
        if (protectZones_.match(cz.name()))
        {
            for (const label celli : cz)
            {
                error[celli] = 0.0;
                nProtected++;
            }
        }
    }

    const auto* stencilPtr =
        oversetMesh_.findObject<cellCellStencilObject>
        (
            cellCellStencilObject::typeName
        );

    if (!stencilPtr)
    {
        return nProtected;
    }

    const labelUList& cellTypes = stencilPtr->cellTypes();

    forAll(cellTypes, celli)
    {
        if (cellTypes[celli] == cellCellStencil::HOLE)
        {
            error[celli] = 0.0;
            nProtected++;
        }
    }

    return nProtected;
}


void Foam::oversetHandler::updateMesh(const mapPolyMesh& mpm)
{
    auto* zonesPtr = oversetMesh_.getObjectPtr<labelIOList>("zoneID");

    if (!zonesPtr)
    {
        return;
    }

    labelIOList& zones = *zonesPtr;
    const labelList& cellMap = mpm.cellMap();

    labelList newZones(cellMap.size(), 0);
    forAll(cellMap, celli)
    {
        const label oldCelli = cellMap[celli];

        // Refined cells inherit the zone of the cell they came from
        if (oldCelli >= 0 && oldCelli < zones.size())
        {
            newZones[celli] = zones[oldCelli];
        }
    }

    zones.transfer(newZones);

    // The stencil itself is discarded by meshObject::updateMesh; the
    // mirrored field is mapped by fvMesh::updateMesh
}


void Foam::oversetHandler::preDistribute()
{
    oldZoneID_ = zoneID();
}


void Foam::oversetHandler::distribute(const mapDistributePolyMesh& map)
{
    map.distributeCellData(oldZoneID_);

    zoneID().transfer(oldZoneID_);

    expireStencil();
}


void Foam::oversetHandler::writeZoneID() const
{
    const auto* zonesPtr = oversetMesh_.findObject<labelIOList>("zoneID");

    if (!zonesPtr)
    {
        return;
    }

    // Written under the current time so that findInstance() picks this up
    // instead of the zoneID belonging to the original mesh
    volScalarField volZoneID
    (
        IOobject
        (
            "zoneID",
            oversetMesh_.time().timeName(),
            oversetMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        oversetMesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll(*zonesPtr, celli)
    {
        volZoneID[celli] = scalar((*zonesPtr)[celli]);
    }

    volZoneID.write();
}


// ************************************************************************* //
