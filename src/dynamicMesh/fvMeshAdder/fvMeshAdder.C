/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "fvMesh.H"
#include "fvMeshAdder.H"
#include "faceCoupleInfo.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(fvMeshAdder, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::fvMeshAdder::calcPatchMap
(
    const label oldStart,
    const label oldSize,
    const labelList& oldToNew,
    const polyPatch& newPatch,
    const label unmappedValue
)
{
    labelList newToOld(newPatch.size(), unmappedValue);

    label newStart = newPatch.start();
    label newSize = newPatch.size();

    for (label i = 0; i < oldSize; i++)
    {
        label newFacei = oldToNew[oldStart+i];

        if (newFacei >= newStart && newFacei < newStart+newSize)
        {
            newToOld[newFacei-newStart] = i;
        }
    }
    return newToOld;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapAddedPolyMesh> Foam::fvMeshAdder::add
(
    fvMesh& mesh0,
    const fvMesh& mesh1,
    const faceCoupleInfo& coupleInfo,
    const bool validBoundary,
    const bool fullyMapped
)
{
    mesh0.clearOut();

    // Resulting merged mesh (polyMesh only!)
    autoPtr<mapAddedPolyMesh> mapPtr
    (
        polyMeshAdder::add
        (
            mesh0,
            mesh1,
            coupleInfo,
            validBoundary
        )
    );
    mapAddedPolyMesh& map = *mapPtr;

    // Adjust the fvMesh part.
    const polyBoundaryMesh& patches = mesh0.boundaryMesh();

    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh0.boundary());
    fvPatches.setSize(patches.size());
    forAll(patches, patchi)
    {
        fvPatches.set(patchi, fvPatch::New(patches[patchi], fvPatches));
    }

    // Do the mapping of the stored fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fvMeshAdder::MapVolFields<scalar>(map, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapVolFields<vector>(map, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapVolFields<sphericalTensor>(map, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapVolFields<symmTensor>(map, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapVolFields<tensor>(map, mesh0, mesh1, fullyMapped);

    fvMeshAdder::MapSurfaceFields<scalar>(map, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapSurfaceFields<vector>(map, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapSurfaceFields<sphericalTensor>
    (
        map, mesh0, mesh1, fullyMapped
    );
    fvMeshAdder::MapSurfaceFields<symmTensor>(map, mesh0, mesh1, fullyMapped);
    fvMeshAdder::MapSurfaceFields<tensor>(map, mesh0, mesh1, fullyMapped);

    fvMeshAdder::MapDimFields<scalar>(map, mesh0, mesh1);
    fvMeshAdder::MapDimFields<vector>(map, mesh0, mesh1);
    fvMeshAdder::MapDimFields<sphericalTensor>(map, mesh0, mesh1);
    fvMeshAdder::MapDimFields<symmTensor>(map, mesh0, mesh1);
    fvMeshAdder::MapDimFields<tensor>(map, mesh0, mesh1);

    return mapPtr;
}


// ************************************************************************* //
