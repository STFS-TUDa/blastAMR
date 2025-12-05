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

#include "cellCountPolicy.H"
#include "addToRunTimeSelectionTable.H"
#include "cloudSupport.H"
#include "messageStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellCountPolicy, 0);
    addToRunTimeSelectionTable(loadPolicy, cellCountPolicy, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCountPolicy::cellCountPolicy
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    loadPolicy(mesh, dict),
    particleCoeff_(dict.getOrDefault<scalar>("particleCoeff", 1.0)),
    minCellsPerProc_(dict.getOrDefault<label>("minCellsPerProc", 5))
{
    Info<< "    particleCoeff: " << particleCoeff_ << endl;
    Info<< "    minCellsPerProc: " << minCellsPerProc_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCountPolicy::~cellCountPolicy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cellCountPolicy::canBalance()
{
    label nParticles = cloudSupport::countParticles(mesh_);
    myLoad_ = mesh_.nCells() + particleCoeff_ * nParticles;

    DebugPout<< "    cells: " << mesh_.nCells()
        << ", particles: " << nParticles
        << ", load: " << myLoad_ << endl;

    myLoadHistory_.set(mesh_.time().timeIndex(), myLoad_);

    scalar globalLoad = returnReduce(myLoad_, sumOp<scalar>());
    scalar idealLoad = globalLoad / scalar(Pstream::nProcs());
    scalar maxImbalance = returnReduce(mag(myLoad_ - idealLoad) / idealLoad, maxOp<scalar>());

    Info<< "Maximum imbalance found = " << 100*maxImbalance << " %" << endl;
    if (maxImbalance < allowedImbalance_)
    {
        return false;
    }
    return true;
}

Foam::scalarField Foam::cellCountPolicy::cellWeights()
{
    // Base weight of 1 per cell + particle contribution
    scalarField weights(mesh_.nCells(), 1.0);

    if (particleCoeff_ > SMALL)
    {
        tmp<labelField> tParticlesPerCell = cloudSupport::particlesPerCell(mesh_);
        const labelField& ppc = tParticlesPerCell();

        forAll(weights, celli)
        {
            weights[celli] += particleCoeff_ * ppc[celli];
        }
    }

    return weights;
}

bool Foam::cellCountPolicy::willBeBeneficial
(
    const labelList& distribution
) {
    // Get cell weights including particle contributions
    scalarField weights = cellWeights();

    // Calculate new load and cell count per processor
    scalarList procLoadNew(Pstream::nProcs(), 0.0);
    labelList procCellsNew(Pstream::nProcs(), 0);
    forAll(distribution, celli)
    {
        procLoadNew[distribution[celli]] += weights[celli];
        procCellsNew[distribution[celli]]++;
    }
    reduce(procLoadNew, sumOp<scalarList>());
    reduce(procCellsNew, sumOp<labelList>());

    if (min(procLoadNew) < SMALL)
    {
        DebugInfo
            << "New distribution results in a load of ~0. Skipping" << endl;
        return false;
    }

    label minCells = min(procCellsNew);
    if (minCells < minCellsPerProc_)
    {
        Info<< "    Not balancing because distribution would leave"
            << " processor with only " << minCells << " cells"
            << " (minimum: " << minCellsPerProc_ << ")" << nl
            << "    Cells per proc: " << procCellsNew << endl;
        return false;
    }

    scalar averageLoadNew = sum(procLoadNew) / scalar(Pstream::nProcs());
    scalar maxDevNew = 0;
    forAll(procLoadNew, proci)
    {
        maxDevNew = max(maxDevNew, mag(procLoadNew[proci] - averageLoadNew) / averageLoadNew);
    }

    // Calculate current imbalance using myLoad_ (already includes particles)
    scalar nGlobalLoad = returnReduce(myLoad_, sumOp<scalar>());
    scalar idealLoad = nGlobalLoad / Pstream::nProcs();
    scalar maxImbalance = returnReduce(mag(myLoad_ - idealLoad) / idealLoad, maxOp<scalar>());

    if (maxDevNew > maxImbalance * 0.99)
    {
        Info
            << "    Not balancing because the new distribution does" << nl
            << "    not improve the load. Skipping" << nl
            << "    old imbalance: " << maxImbalance << nl
            << "    new imbalance: " << maxDevNew << nl
            << endl;
        return false;
    }
    return true;
}

// ************************************************************************* //
