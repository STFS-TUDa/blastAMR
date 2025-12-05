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
    Foam::amrCore

Description
    Core implementation of Adaptive Mesh Refinement (AMR) and Load Balancing
    functionality. This class provides the shared logic used by both:
    - adaptiveFvMesh (dynamicFvMesh-based approach)
    - meshRefinement/loadBalance functionObjects

    The class manages:
    - Error estimation for refinement criteria
    - Mesh refinement/unrefinement via fvMeshRefiner
    - Flux correction after topology changes
    - Cloud position storage and remapping
    - Load balancing via fvMeshBalance

SourceFiles
    amrCore.C

\*---------------------------------------------------------------------------*/

#ifndef amrCore_H
#define amrCore_H

#include "fvMesh.H"
#include "fvMeshRefiner.H"
#include "errorEstimator.H"
#include "HashTable.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class amrCore Declaration
\*---------------------------------------------------------------------------*/

class amrCore
{
protected:

    // Protected Data

        //- Reference to the mesh
        fvMesh& mesh_;

        //- Error estimator for refinement criteria
        autoPtr<errorEstimator> error_;

        //- Mesh refiner
        autoPtr<fvMeshRefiner> refiner_;

        //- Mapping between flux field names and velocity field names
        //  for flux correction after topology changes
        HashTable<word> correctFluxes_;

        //- Enable expiration of sampledSurface caches after mesh topology changes
        //  This is opt-in (default false) because it uses private member access
        bool expireSampledSurfacesOnLB_;


    // Protected Member Functions

        //- No copy construct
        amrCore(const amrCore&) = delete;

        //- No copy assignment
        void operator=(const amrCore&) = delete;


public:

    //- Runtime type information
    TypeName("amrCore");


    // Constructors

        //- Construct from mesh and dictionary
        amrCore(fvMesh& mesh, const dictionary& dict);

        //- Construct from mesh only (deferred initialization)
        explicit amrCore(fvMesh& mesh);


    //- Destructor
    virtual ~amrCore() = default;


    // Member Functions

        // Initialization

            //- Initialize AMR components from dictionary
            void initializeAMR(const dictionary& dict);

            //- Initialize load balancing from dictionary
            void initializeLB(const dictionary& dict);

            //- Read/update settings from dictionary
            void readDict(const dictionary& dict);

            //- Read flux correction settings from dictionary
            void readCorrectFluxes(const dictionary& dict);


        // Access

            //- Return const access to mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return non-const access to mesh
            fvMesh& mesh()
            {
                return mesh_;
            }

            //- Return const access to error estimator
            const errorEstimator& error() const
            {
                return error_();
            }

            //- Return non-const access to error estimator
            errorEstimator& error()
            {
                return error_();
            }

            //- Return const access to refiner
            const fvMeshRefiner& refiner() const
            {
                return refiner_();
            }

            //- Return non-const access to refiner
            fvMeshRefiner& refiner()
            {
                return refiner_();
            }

            //- Check if AMR is initialized
            bool amrInitialized() const
            {
                return error_.valid() && refiner_.valid();
            }

            //- Check if load balancing is initialized
            bool lbInitialized() const
            {
                return refiner_.valid();
            }

            //- Return if sampledSurface cache expiration is enabled
            bool expireSampledSurfacesOnLB() const
            {
                return expireSampledSurfacesOnLB_;
            }


        // Operations

            //- Perform mesh refinement
            //  Returns true if mesh topology changed
            bool refine();

            //- Perform load balancing
            //  Returns true if mesh was redistributed
            bool balance();

            //- Correct fluxes on new/modified faces after topology change
            void correctFluxes(const mapPolyMesh& map);

            //- Update after mesh topology change
            void updateMesh(const mapPolyMesh& map);

            //- Distribute data after mesh redistribution
            void distribute(const mapDistributePolyMesh& map);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
