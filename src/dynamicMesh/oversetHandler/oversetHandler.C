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
#include "oversetPolyPatch.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetHandler, 0);
    defineRunTimeSelectionTable(oversetHandler, mesh);
}

const Foam::word Foam::oversetHandler::zoneIDFieldName("oversetZoneID");


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
    oversetMesh_(mesh)
{
    // Mirror zoneID while the mesh still matches the zoneID field on disk
    zoneIDField();
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

    Info<< "oversetHandler: registering " << handlerPtr->type()
        << " for overset mesh " << mesh.name() << endl;

    oversetHandler* handler = handlerPtr.ptr();
    handler->store();

    return handler;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::volScalarField& Foam::oversetHandler::zoneIDField()
{
    auto* fieldPtr =
        oversetMesh_.getObjectPtr<volScalarField>(zoneIDFieldName);

    if (!fieldPtr)
    {
        fieldPtr = new volScalarField
        (
            IOobject
            (
                zoneIDFieldName,
                oversetMesh_.time().timeName(),
                oversetMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            oversetMesh_,
            dimensionedScalar(dimless, Zero),
            zeroGradientFvPatchScalarField::typeName
        );
        fieldPtr->store();

        // Seed from the stencil's own zoneID, read off the mesh as it is
        // now. Everything after this point rides on the field
        const labelIOList& zones = cellCellStencil::zoneID(oversetMesh_);

        forAll(zones, celli)
        {
            (*fieldPtr)[celli] = scalar(zones[celli]);
        }
        fieldPtr->correctBoundaryConditions();
    }

    return *fieldPtr;
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

Foam::label Foam::oversetHandler::protectCells(volScalarField& error) const
{
    const auto* stencilPtr =
        oversetMesh_.findObject<cellCellStencilObject>
        (
            cellCellStencilObject::typeName
        );

    if (!stencilPtr)
    {
        return 0;
    }

    const labelUList& cellTypes = stencilPtr->cellTypes();

    label nProtected = 0;
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


void Foam::oversetHandler::sync()
{
    const volScalarField& zoneIDf = zoneIDField();

    auto* zonesPtr = oversetMesh_.getObjectPtr<labelIOList>("zoneID");

    if (zonesPtr)
    {
        labelIOList& zones = *zonesPtr;
        zones.setSize(zoneIDf.size());

        forAll(zones, celli)
        {
            zones[celli] = label(zoneIDf[celli] + 0.5);
        }
    }

    expireStencil();
}


void Foam::oversetHandler::writeZoneID() const
{
    const auto* fieldPtr =
        oversetMesh_.findObject<volScalarField>(zoneIDFieldName);

    if (!fieldPtr)
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
        *fieldPtr
    );

    volZoneID.write();
}


// ************************************************************************* //
