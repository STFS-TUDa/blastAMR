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
    Foam::adaptiveOversetFvMesh

Description
    dynamicOversetFvMesh with adaptive mesh refinement and load balancing.

    Inherits the whole overset machinery (extended ldu addressing, implicit
    and explicit interpolation, mesh motion through the motion solver list)
    from dynamicOversetFvMesh and adds blastAMR's amrCore on top, so a
    single dynamicMeshDict drives motion, refinement and balancing:

    \verbatim
        dynamicFvMesh   adaptiveOversetFvMesh;

        solver          multiSolidBodyMotionSolver;
        multiSolidBodyMotionSolverCoeffs { ... }

        adaptiveOversetFvMeshCoeffs
        {
            errorEstimator  ...;
            balance         yes;
            ...
        }
    \endverbatim

    Overset zone book-keeping across refinement and balancing is done by the
    oversetHandler registered by amrCore.

SourceFiles
    adaptiveOversetFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef adaptiveOversetFvMesh_H
#define adaptiveOversetFvMesh_H

#include "dynamicOversetFvMesh.H"
#include "amrCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class adaptiveOversetFvMesh Declaration
\*---------------------------------------------------------------------------*/

class adaptiveOversetFvMesh
:
    public dynamicOversetFvMesh
{
        //- Return the dynamicMeshDict IOobject
        static IOobject dynamicMeshDictIOobject(const IOobject& io);


protected:

        //- The dynamicMeshDict
        IOdictionary dynamicMeshDict_;

        //- Core AMR/LB functionality
        amrCore amrCore_;

        //- Do mesh updates only once per timestep
        label currentTimeIndex_;


    // Protected Member Functions

        //- Read the AMR settings
        void readDict();

        //- Overload update mesh to include other methods
        virtual void updateMesh(const mapPolyMesh& mpm);


public:

    //- Runtime type information
    TypeName("adaptiveOversetFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit adaptiveOversetFvMesh(const IOobject& io);

        //- No copy construct
        adaptiveOversetFvMesh(const adaptiveOversetFvMesh&) = delete;


    //- Destructor
    virtual ~adaptiveOversetFvMesh() = default;


    // Member Functions

        //- Return the dynamicMeshDict
        const dictionary& dynamicMeshDict() const
        {
            return dynamicMeshDict_;
        }

        //- Direct access to the refinement engine
        const fvMeshRefiner& refiner() const
        {
            return amrCore_.refiner();
        }

        //- Map all fields in time using given map
        virtual void mapFields(const mapPolyMesh& mpm);

        //- Update the mesh for motion, topology change and balancing
        virtual bool update();

        //- Refine the mesh
        bool refine();


    // Writing

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool valid
        ) const;


    // Member Operators

        //- No copy assignment
        void operator=(const adaptiveOversetFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
