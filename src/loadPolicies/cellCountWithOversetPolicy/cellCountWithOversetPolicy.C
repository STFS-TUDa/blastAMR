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

#include "cellCountWithOversetPolicy.H"
#include "addToRunTimeSelectionTable.H"
#include "cellCellStencilObject.H"
#include "oversetPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellCountWithOversetPolicy, 0);
    addToRunTimeSelectionTable
    (
        loadPolicy,
        cellCountWithOversetPolicy,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCountWithOversetPolicy::cellCountWithOversetPolicy
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    cellCountPolicy(mesh, dict),
    holeWeight_(dict.getOrDefault<scalar>("holeWeight", 0.01)),
    interpolatedWeight_(dict.getOrDefault<scalar>("interpolatedWeight", 1.0)),
    porousWeight_(dict.getOrDefault<scalar>("porousWeight", 1.0)),
    specialWeight_(dict.getOrDefault<scalar>("specialWeight", 1.0))
{
    // A negative weight would corrupt the decomposition rather than just
    // bias it, so refuse it outright
    const scalar minWeight
    (
        min
        (
            min(holeWeight_, interpolatedWeight_),
            min(porousWeight_, specialWeight_)
        )
    );

    if (minWeight < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Overset cell weights must not be negative, got"
            << " holeWeight " << holeWeight_
            << ", interpolatedWeight " << interpolatedWeight_
            << ", porousWeight " << porousWeight_
            << ", specialWeight " << specialWeight_
            << exit(FatalIOError);
    }

    Info<< "    holeWeight: " << holeWeight_
        << ", interpolatedWeight: " << interpolatedWeight_
        << ", porousWeight: " << porousWeight_
        << ", specialWeight: " << specialWeight_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField Foam::cellCountWithOversetPolicy::cellWeights()
{
    scalarField weights(cellCountPolicy::cellWeights());

    bool isOverset = false;
    for (const polyPatch& pp : mesh_.boundaryMesh())
    {
        if (isA<oversetPolyPatch>(pp))
        {
            isOverset = true;
            break;
        }
    }

    if (!isOverset)
    {
        // Nothing overset about this mesh, so plain cellCount it is
        return weights;
    }

    // Construct rather than look up. The stencil is discarded on every
    // topology change and rebuilt lazily on first use, which happens after
    // the balancing decision - so looking it up finds nothing precisely when
    // the weights are wanted. Building it here is not wasted work: the solver
    // needs it in the same timestep and picks up the cached one, except on
    // the steps that do redistribute, which discard it either way
    const labelUList& cellTypes = Stencil::New(mesh_).cellTypes();

    if (cellTypes.size() != weights.size())
    {
        // Sized for another topology, so of no use here
        return weights;
    }

    forAll(weights, celli)
    {
        scalar factor = 1.0;

        switch (cellTypes[celli])
        {
            case cellCellStencil::HOLE:
                factor = holeWeight_;
                break;

            case cellCellStencil::INTERPOLATED:
                factor = interpolatedWeight_;
                break;

            case cellCellStencil::POROUS:
                factor = porousWeight_;
                break;

            case cellCellStencil::SPECIAL:
                factor = specialWeight_;
                break;

            default:
                // CALCULATED keeps the cellCount weight
                break;
        }

        // Scale the cell contribution only. The cellCount weight is
        // 1 + particleCoeff*nParticles, and a parcel costs the same to track
        // wherever it sits, so the particle part is carried over untouched
        weights[celli] = factor + (weights[celli] - 1.0);
    }

    return weights;
}


Foam::scalar Foam::cellCountWithOversetPolicy::localLoad()
{
    // The load has to be the same measure the decomposition is given,
    // otherwise the decision to balance and the placement of cells disagree
    return sum(cellWeights());
}


// ************************************************************************* //
