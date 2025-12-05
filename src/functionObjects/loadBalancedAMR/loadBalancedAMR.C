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

#include "loadBalancedAMR.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(loadBalancedAMR, 0);
    addToRunTimeSelectionTable(functionObject, loadBalancedAMR, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::dynamicFvMesh& Foam::functionObjects::loadBalancedAMR::validateDynamicMesh
(
    const fvMesh& mesh
)
{
    // Try to cast to dynamicFvMesh
    dynamicFvMesh* dynMeshPtr =
        const_cast<dynamicFvMesh*>(dynamic_cast<const dynamicFvMesh*>(&mesh));

    if (!dynMeshPtr)
    {
        FatalErrorInFunction
            << "The loadBalancedAMR functionObject requires a dynamicFvMesh type."
            << nl << nl
            << "Your current mesh type is: " << mesh.polyMesh::type() << nl << nl
            << "To use the loadBalancedAMR functionObject, you must run "
            << "a solver that has a dynamicFvMesh (not a fvMesh one)" << nl
            << exit(FatalError);
    }

    return *dynMeshPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::loadBalancedAMR::loadBalancedAMR
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    dynMesh_(validateDynamicMesh(mesh_)),
    amrCore_(dynMesh_),
    currentTimeIndex_(-1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::loadBalancedAMR::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    bool doRefine = dict.getOrDefault("refine", true);
    bool doUnrefine = dict.getOrDefault("unrefine", true);
    bool doBalance = dict.getOrDefault("balance", false);

    // Only initialize AMR if refinement/unrefinement is enabled
    if (doRefine || doUnrefine)
    {
        amrCore_.initializeAMR(dict);
    }

    // Initialize load balancing if balance is enabled
    if (doBalance)
    {
        amrCore_.initializeLB(dict);
    }

    return true;
}


bool Foam::functionObjects::loadBalancedAMR::execute()
{
    // Prevent double execution per timestep
    if (currentTimeIndex_ >= time_.timeIndex())
    {
        return true;
    }
    currentTimeIndex_ = time_.timeIndex();

    // Check and perform refinement
    if (amrCore_.amrInitialized())
    {
        // Re-read dictionary for on-the-fly modifications
        // (same behavior as adaptiveFvMesh::refine)
        // Note: Could add rereading from functionObject dict if needed

        if (amrCore_.refine())
        {
            Info<< name()
                << ": Mesh topology changed (refinement)" << endl;
        }
    }

    // Check and perform load balancing (only in parallel)
    if (amrCore_.lbInitialized() && Pstream::parRun())
    {
        if (amrCore_.balance())
        {
            Info<< name()
                << ": Mesh redistributed across processors" << endl;
        }
    }

    return true;
}


bool Foam::functionObjects::loadBalancedAMR::write()
{
    // Write error field if desired
    if (amrCore_.amrInitialized())
    {
        const volScalarField& errorField = amrCore_.error().error();
        if (errorField.writeOpt() == IOobject::AUTO_WRITE)
        {
            errorField.write();
        }
    }

    // Report load balance statistics in debug mode
    if (Pstream::parRun() && debug)
    {
        label nLocalCells = mesh_.nCells();
        label nTotalCells = returnReduce(nLocalCells, sumOp<label>());
        scalar avgCells = scalar(nTotalCells) / Pstream::nProcs();
        scalar imbalance = mag(nLocalCells - avgCells) / max(avgCells, SMALL);

        Pout<< name()
            << ": Local cells = " << nLocalCells
            << ", avg = " << avgCells
            << ", imbalance = " << 100*imbalance << "%" << endl;
    }

    return true;
}


// ************************************************************************* //
