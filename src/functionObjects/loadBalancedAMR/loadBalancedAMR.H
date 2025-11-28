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
    Foam::functionObjects::loadBalancedAMR

Group
    grpMeshFunctionObjects

Description
    FunctionObject that performs adaptive mesh refinement (AMR) and dynamic
    load balancing during runtime. This allows using AMR/LB with any
    dynamicFvMesh type, not just adaptiveFvMesh.

    The dictionary settings are identical to those used in dynamicMeshDict
    for adaptiveFvMesh, enabling seamless switching between the two approaches.

Usage
    Example of function object specification:
    \verbatim
    functions
    {
        loadBalancedAMR
        {
            type            loadBalancedAMR;
            libs            (libamrFunctionObjects);

            // Execution control
            executeControl  timeStep;
            executeInterval 1;

            // --- AMR Settings (same as adaptiveFvMesh) ---

            // Refiner type: hexRefiner, polyRefiner, directionalRefiner
            refiner         hexRefiner;

            // Error estimator: delta, gradient, densityGradient, coded, etc.
            errorEstimator  delta;

            // Field for error estimation (for delta/gradient estimators)
            field           alpha.water;

            // Refinement on/off
            refine          true;
            unrefine        true;

            // Timing control
            refineInterval  1;
            beginRefine     0;
            endRefine       1e10;
            beginUnrefine   0;
            endUnrefine     1e10;

            // Refinement levels
            maxRefinement   3;
            nRefinementBufferLayers 1;
            nUnrefinementBufferLayers 1;

            // Error thresholds
            lowerRefineLevel  0.01;
            upperRefineLevel  0.99;
            unrefineLevel     10;

            // Cell limits
            maxCells        1000000;

            // Protected patches (no refinement near these)
            protectedPatches ();

            // Flux correction
            correctFluxes   ((phi U) (phi_0 U_0));

            // --- Load Balancing Settings ---

            balance             true;
            allowableImbalance  0.1;
            balanceInterval     1;
            beginBalance        0;
            endBalance          1e10;
        }
    }
    \endverbatim

Note
    - Requires a dynamicFvMesh type (not staticFvMesh)
    - Settings are identical to adaptiveFvMesh dynamicMeshDict format
    - Load balancing only active in parallel runs

See also
    Foam::adaptiveFvMesh
    Foam::amrCore

SourceFiles
    loadBalancedAMR.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_loadBalancedAMR_H
#define functionObjects_loadBalancedAMR_H

#include "fvMeshFunctionObject.H"
#include "amrCore.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class loadBalancedAMR Declaration
\*---------------------------------------------------------------------------*/

class loadBalancedAMR
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Reference to dynamicFvMesh (validated at construction)
        dynamicFvMesh& dynMesh_;

        //- Core AMR/LB functionality
        amrCore amrCore_;

        //- Current time index to prevent double operations
        label currentTimeIndex_;


    // Private Member Functions

        //- Check and return dynamicFvMesh reference
        //  Throws FatalError if mesh is not a dynamicFvMesh
        static dynamicFvMesh& validateDynamicMesh(const fvMesh& mesh);

        //- No copy construct
        loadBalancedAMR(const loadBalancedAMR&) = delete;

        //- No copy assignment
        void operator=(const loadBalancedAMR&) = delete;


public:

    //- Runtime type information
    TypeName("loadBalancedAMR");


    // Constructors

        //- Construct from Time and dictionary
        loadBalancedAMR
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~loadBalancedAMR() = default;


    // Member Functions

        //- Read the dictionary settings
        virtual bool read(const dictionary& dict);

        //- Execute AMR and load balancing
        virtual bool execute();

        //- Write diagnostic fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
