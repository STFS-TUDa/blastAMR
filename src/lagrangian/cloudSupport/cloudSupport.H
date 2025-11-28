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
    Foam::cloudSupport

Description
    Utility functions for Lagrangian cloud handling during mesh refinement
    and load balancing operations.

    Built-in support includes:
    - passiveParticleCloud (basic tracking particles)
    - basicKinematicCloud (kinematic particles with forces)
    - basicKinematicCollidingCloud (with DEM collisions)
    - basicKinematicMPPICCloud (MPPIC method)
    - basicThermoCloud (heat transfer particles)
    - basicReactingCloud (evaporating/reacting particles)
    - basicReactingMultiphaseCloud (multiphase reactions)

    Supported operations:
    - storeGlobalPositions() - Store particle coordinates before topology changes
    - autoMapClouds() - Remap particles after refinement (mapPolyMesh)
    - distributeClouds() - Redistribute particles during load balancing
    - relocateClouds() - Recreate particles after mesh redistribution

See Also
    Foam::cloudHandler

SourceFiles
    cloudSupport.C

\*---------------------------------------------------------------------------*/

#ifndef cloudSupport_H
#define cloudSupport_H

#include "fvMesh.H"
#include "Cloud.H"
#include "passiveParticle.H"
#include "mapPolyMesh.H"
#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cloudSupport
{

// NOTE: There is an "expected" order of operations:
// storeGlobalPositions -> distributClouds -> mesh.distribute -> relocateClouds
// so, might be better to be called in preDistribute/postDistribute routines

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Store global positions for all clouds in the mesh
void storeGlobalPositions(const fvMesh& mesh);

//- AutoMap all clouds after mesh topology changes
void autoMapClouds
(
    const fvMesh& mesh,
    const mapPolyMesh& map
);

//- Distribute all clouds to new processors
//  Transfers particles between processors during domain redistribution
//  Must be called BEFORE mesh distribution with the cell->processor mapping
void distributeClouds
(
    const fvMesh& mesh,
    const labelList& distribution
);

//- Relocate particles after distribution
//  Must be called AFTER mesh distribution to find new cells
void relocateClouds(const fvMesh& mesh);

//- Count total number of particles across all clouds
label countParticles(const fvMesh& mesh);

//- Count particles per cell
tmp<labelField> particlesPerCell(const fvMesh& mesh);

//- Return list of cloud names present in the mesh
wordList cloudNames(const fvMesh& mesh);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cloudSupport
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
