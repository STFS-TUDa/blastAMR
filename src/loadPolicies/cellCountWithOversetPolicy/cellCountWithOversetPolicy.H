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

Class
    Foam::cellCountWithOversetPolicy

Description
    cellCount load policy that prices overset cells by what they cost.

    cellCount weights every cell the same. On an overset mesh that is wrong:
    a hole cell is masked out of the solution (its matrix row is forced, it
    carries no physics) yet counts as a full cell, so a processor owning a
    large part of the hole region reports as fully loaded while doing almost
    no work. The balancer then equalises cell counts into a state that is
    still imbalanced in compute. An acceptor cell has the opposite problem:
    the extended overset addressing gives it extra matrix faces and, for
    remote donors, communication, so it costs more than a plain cell.

    This policy multiplies the cellCount weight of each cell by a factor
    chosen per overset cell type, and reports the load as the local sum of
    those weights so that the decision to balance and the placement of cells
    use the same measure. Particle contributions are left at full weight:
    tracking a parcel costs the same wherever it sits.

    All factors are user-configurable, and default to changing nothing except
    discounting holes - the one case that is unambiguous:

    \verbatim
        loadPolicy          cellCountWithOverset;

        holeWeight          0.01;   // masked out, no physics solved
        interpolatedWeight  1.0;    // raise above 1: extra faces + comms
        porousWeight        1.0;    // no donors found, special treatment
        specialWeight       1.0;    // hole that acts as a donor
    \endverbatim

    Calculated cells always keep their cellCount weight. A mesh with no
    stencil registered, or a stencil sized for another topology, falls back
    to plain cellCount behaviour.

SourceFiles
    cellCountWithOversetPolicy.C

\*---------------------------------------------------------------------------*/

#ifndef cellCountWithOversetPolicy_H
#define cellCountWithOversetPolicy_H

#include "cellCountPolicy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cellCountWithOversetPolicy Declaration
\*---------------------------------------------------------------------------*/

class cellCountWithOversetPolicy
:
    public cellCountPolicy
{
protected:

    // Protected data

        //- Weight factor for hole cells
        scalar holeWeight_;

        //- Weight factor for interpolated (acceptor) cells
        scalar interpolatedWeight_;

        //- Weight factor for porous cells (no donors found)
        scalar porousWeight_;

        //- Weight factor for special cells (hole acting as a donor)
        scalar specialWeight_;


    // Protected Member Functions

        //- Local sum of the weights, so that the balancing decision and the
        //  decomposition see the same load
        virtual scalar localLoad();


public:

    //- Runtime type information
    TypeName("cellCountWithOverset");


    // Constructors

        //- Construct from mesh and dictionary
        cellCountWithOversetPolicy
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- No copy construct
        cellCountWithOversetPolicy(const cellCountWithOversetPolicy&) = delete;


    //- Destructor
    virtual ~cellCountWithOversetPolicy() = default;


    // Member Functions

        //- Cell weights, the cellCount weights scaled per overset cell type
        virtual scalarField cellWeights();


    // Member Operators

        //- No copy assignment
        void operator=(const cellCountWithOversetPolicy&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
