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
    Foam::loadPolicy

SourceFiles
    loadPolicy.C

Description
    A base class for tracking MPI-processor-local CPU load.

\*---------------------------------------------------------------------------*/

#ifndef loadPolicy_H
#define loadPolicy_H

#include "fvMesh.H"
#include "autoPtr.H"
#include "HashTable.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class loadPolicy Declaration
\*---------------------------------------------------------------------------*/


class loadPolicy
{
protected:

    //- Mesh refence
    const fvMesh& mesh_;

    //- configuration
    const dictionary& dict_;

    //- allowed imbalance as a fraction of mean of processor loads
    mutable scalar allowedImbalance_;

    //- Current load value for this processor
    mutable scalar myLoad_;

    //- Load history for this processor
    //  indexed with time index
    mutable HashTable<scalar, label> myLoadHistory_;

public:

    //- Runtime type information
    TypeName("loadPolicy");

    declareRunTimeSelectionTable
    (
        autoPtr,
        loadPolicy,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (mesh, dict)
    );

    // Constructors

        //- Construct from IOobject
        loadPolicy
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        loadPolicy(const loadPolicy&) = delete;

        static autoPtr<loadPolicy> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~loadPolicy();


    // Member Functions

        //- Return the dynamicMeshDict
        const dictionary& dynamicMeshDict() const
        {
            return dict_;
        }

        //- Return allowed imbalance
        scalar allowedImbalance() const {
            return allowedImbalance_;
        }

        //- Return current load value
        scalar myLoad() const {
            if (myLoad_ == -1)
            {
                FatalErrorInFunction
                    << "The load-balancer is requesting a load before it is computed... This is not allowed."
                    << abort(FatalError);
            }
            return myLoad_;
        }

        const HashTable<scalar, label>& myLoadHistory() const{
            return myLoadHistory_;
        }

        //- Compute loads and see if should balance
        virtual bool canBalance() = 0;

        //- Compute and return cell weights for decomposition
        virtual scalarField cellWeights() = 0;

        //- Predict if LB would be viable
        virtual bool willBeBeneficial
        (
            const labelList& distribution
        ) = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const loadPolicy&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
