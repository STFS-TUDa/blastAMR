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
    Foam::cpuLoadPolicy

SourceFiles
    cpuLoadPolicy.C

Description
    A load tracked based on local CPU time. Imbalance is calculated irt.
    ideal average (computation-related) CPU execution time per processor.

Note
    Relies on a static global for MPI profiling 

\*---------------------------------------------------------------------------*/

#ifndef cpuLoadPolicy_H
#define cpuLoadPolicy_H

#include "loadPolicy.H"
#include "mpiOverrides.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class cpuLoadPolicy Declaration
\*---------------------------------------------------------------------------*/

class cpuLoadPolicy
:   public loadPolicy
{
protected:
    // Maximum cycle length (load is tracked over this many timesteps)
    label maxCycleLength_;

    // Time index when stats were last reset
    label lastResetTimeIndex_;

    // Is MPI call interception active
    bool isActive_;

    // Time unit for reporting (nano, micro, milli)
    TimeUnit timeUnit_;

public:

    //- Runtime type information
    TypeName("cpuLoad");

    // Constructors

        //- Construct from IOobject
        cpuLoadPolicy
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cpuLoadPolicy(const cpuLoadPolicy&) = delete;

    //- Destructor
    virtual ~cpuLoadPolicy();


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

        inline bool stillSameCycle() const {
            if (myLoadHistory_.empty()) return false;
            return mesh_.time().timeIndex() == myLoadHistory_.toc().last();
        }
    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cpuLoadPolicy&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
