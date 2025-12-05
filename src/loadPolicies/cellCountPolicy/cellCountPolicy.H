/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020 Synthetik Applied Technologies: |   Modified original
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


Class
    Foam::cellCountPolicy

SourceFiles
    cellCountPolicy.C

Description
    A load tracked based on local cell count. Imbalance is calculated irt.
    ideal average cell count per processor.

\*---------------------------------------------------------------------------*/

#ifndef cellCountPolicy_H
#define cellCountPolicy_H

#include "loadPolicy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class cellCountPolicy Declaration
\*---------------------------------------------------------------------------*/

class cellCountPolicy
:   public loadPolicy
{
    // Private data

        //- Coefficient for particle contribution to load
        //  load = cellCount + particleCoeff * particleCount
        scalar particleCoeff_;

        //- Minimum number of cells per processor to avoid degenerate meshes
        //  If a proposed distribution would leave any processor with fewer
        //  cells than this, the distribution is rejected
        label minCellsPerProc_;

public:

    //- Runtime type information
    TypeName("cellCount");

    // Constructors

        //- Construct from IOobject
        cellCountPolicy
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cellCountPolicy(const cellCountPolicy&) = delete;

    //- Destructor
    virtual ~cellCountPolicy();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool canBalance();

        //- Compute and return cell weights for decomposition
        virtual scalarField cellWeights();

        //- Predict if LB would be viable
        virtual bool willBeBeneficial
        (
            const labelList& distribution
        );

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellCountPolicy&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
