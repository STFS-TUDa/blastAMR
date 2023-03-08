/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
#include "polyTopoChange.H"

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

    const label newStart = newPatch.start();
    const label newSize = newPatch.size();

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


Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshAdder::add
(
    const label myProci,            // index of mesh to modify
    UPtrList<fvMesh>& fvMeshes,
    const labelList& oldFaceOwner,  // face owner for myProci mesh

    // Coupling info
    const labelListList& localBoundaryFace,
    const labelListList& remoteFaceProc,
    const labelListList& remoteBoundaryFace,

    labelListList& constructPatchMap,
    labelListList& constructCellMap,
    labelListList& constructFaceMap,
    labelListList& constructPointMap
)
{
    // Do in-place addition. Modifies fvMeshes[myProci]

    UPtrList<polyMesh> meshes(fvMeshes.size());
    forAll(fvMeshes, proci)
    {
        if (fvMeshes.set(proci))
        {
            meshes.set(proci, &fvMeshes[proci]);
        }
    }


    // All matched faces assumed to have vertex0 matched
    labelListList remoteFaceStart(meshes.size());
    forAll(localBoundaryFace, proci)
    {
        const labelList& procFaces = localBoundaryFace[proci];
        remoteFaceStart[proci].setSize(procFaces.size(), 0);
    }

    // Assume all meshes have same global patches
    labelListList patchMap(meshes.size());
    labelListList pointZoneMap(meshes.size());
    labelListList faceZoneMap(meshes.size());
    labelListList cellZoneMap(meshes.size());
    forAll(meshes, proci)
    {
        if (meshes.set(proci))
        {
            const polyMesh& mesh = meshes[proci];
            patchMap[proci] = identity(mesh.boundaryMesh().nNonProcessor());
            pointZoneMap[proci] = identity(mesh.pointZones().size());
            faceZoneMap[proci] = identity(mesh.faceZones().size());
            cellZoneMap[proci] = identity(mesh.cellZones().size());
        }
    }


    // Swap myProci to 0th element
    if (myProci != 0)
    {
        polyMesh* pm0 = meshes.get(0);
        polyMesh* pmi = meshes.get(myProci);
        meshes.set(0, pmi);
        meshes.set(myProci, pm0);

        fvMesh* fvm0 = fvMeshes.get(0);
        fvMesh* fvmi = fvMeshes.get(myProci);
        fvMeshes.set(0, fvmi);
        fvMeshes.set(myProci, fvm0);

        //Pout<< "swapped from 0 to " << myProci << endl;
        //forAll(meshes, meshi)
        //{
        //    Pout<< "meshi:" << meshi << endl;
        //    if (meshes.set(meshi))
        //    {
        //        Pout<< "    nCells:" << meshes[meshi].nCells() << endl;
        //    }
        //}


        // Swap (and renumber) patch face information
        labelListList& lbf = const_cast<labelListList&>(localBoundaryFace);
        std::swap(lbf[0], lbf[myProci]);

        labelListList& rfp = const_cast<labelListList&>(remoteFaceProc);
        std::swap(rfp[0], rfp[myProci]);
        forAll(rfp, proci)
        {
            for (label& proc : rfp[proci])
            {
                if (proc == 0) proc = myProci;
                else if (proc == myProci) proc = 0;
            }
        }
        labelListList& rbf = const_cast<labelListList&>(remoteBoundaryFace);
        std::swap(rbf[0], rbf[myProci]);

        labelListList& rfs = const_cast<labelListList&>(remoteFaceStart);
        std::swap(rfs[0], rfs[myProci]);

        // Swap optional renumbering maps
        std::swap(patchMap[0], patchMap[myProci]);
        std::swap(pointZoneMap[0], pointZoneMap[myProci]);
        std::swap(faceZoneMap[0], faceZoneMap[myProci]);
        std::swap(cellZoneMap[0], cellZoneMap[myProci]);
    }

    polyTopoChange meshMod(meshes[0].boundaryMesh().size(), true);
    // Collect statistics for sizing
    label nCells = 0;
    label nFaces = 0;
    label nPoints = 0;
    forAll(meshes, proci)
    {
        if (meshes.set(proci))
        {
            const polyMesh& mesh = meshes[proci];
            nCells += mesh.nCells();
            nFaces += mesh.nFaces();
            nPoints += mesh.nPoints();
        }
    }
    meshMod.setCapacity(nPoints, nFaces, nCells);

    // Add all cells in meshes' order
    polyMeshAdder::add
    (
        meshes,
        patchMap,

        // Information on (one-to-one) boundary faces to stitch
        localBoundaryFace,
        remoteFaceProc,
        remoteBoundaryFace,
        remoteFaceStart,

        pointZoneMap,
        faceZoneMap,
        cellZoneMap,

        meshMod,

        constructCellMap,
        constructFaceMap,
        constructPointMap
    );

    // Replace mesh
    autoPtr<mapPolyMesh> mapPtr
    (
        meshMod.changeMesh
        (
            fvMeshes[0],        // note: still swapped to position 0
            false               // no inflation
        )
    );

    // Update fields. Note that this tries to interpolate from the mesh
    // before adding all the remote meshes so is quite wrong.
    //fvMeshes[0.updateMesh(mapPtr());
    // Update polyMesh but not fvMesh to avoid mapping all the fields
    fvMeshes[0].polyMesh::updateMesh(mapPtr());

    // Now reverseFaceMap contains from order in which face was added
    // to mesh face (after e.g. upper-triangular ordering)

    // Renumber output of any mesh ordering done by changeMesh
    for (labelList& cellMap : constructCellMap)
    {
        inplaceRenumber(mapPtr().reverseCellMap(), cellMap);
    }
    for (labelList& faceMap : constructFaceMap)
    {
        inplaceRenumber(mapPtr().reverseFaceMap(), faceMap);
    }
    for (labelList& pointMap : constructPointMap)
    {
        inplaceRenumber(mapPtr().reversePointMap(), pointMap);
    }

    // constructPatchMap is patchMap with -1 for the removed
    // patches
    forAll(meshes, meshi)
    {
        if (meshes.set(meshi))
        {
            constructPatchMap[meshi] = patchMap[meshi];
            constructPatchMap[meshi].setSize
            (
                meshes[meshi].boundaryMesh().size(),
                -1
            );
        }
    }


    // Map all fields
    fvMeshAdder::MapVolFields<scalar>
    (
        fvMeshes,
        mapPtr().oldPatchStarts(),
        mapPtr().oldPatchSizes(),
        patchMap,
        constructCellMap,
        constructFaceMap,
        constructPointMap
    );
    fvMeshAdder::MapVolFields<vector>
    (
        fvMeshes,
        mapPtr().oldPatchStarts(),
        mapPtr().oldPatchSizes(),
        patchMap,
        constructCellMap,
        constructFaceMap,
        constructPointMap
    );
    fvMeshAdder::MapVolFields<sphericalTensor>
    (
        fvMeshes,
        mapPtr().oldPatchStarts(),
        mapPtr().oldPatchSizes(),
        patchMap,
        constructCellMap,
        constructFaceMap,
        constructPointMap
    );
    fvMeshAdder::MapVolFields<symmTensor>
    (
        fvMeshes,
        mapPtr().oldPatchStarts(),
        mapPtr().oldPatchSizes(),
        patchMap,
        constructCellMap,
        constructFaceMap,
        constructPointMap
    );
    fvMeshAdder::MapVolFields<tensor>
    (
        fvMeshes,
        mapPtr().oldPatchStarts(),
        mapPtr().oldPatchSizes(),
        patchMap,
        constructCellMap,
        constructFaceMap,
        constructPointMap
    );
    fvMeshAdder::MapSurfaceFields<scalar>
    (
        fvMeshes,
        oldFaceOwner,
        mapPtr().oldPatchStarts(),
        mapPtr().oldPatchSizes(),
        patchMap,
        constructCellMap,
        constructFaceMap,
        constructPointMap
    );
    fvMeshAdder::MapSurfaceFields<vector>
    (
        fvMeshes,
        oldFaceOwner,
        mapPtr().oldPatchStarts(),
        mapPtr().oldPatchSizes(),
        patchMap,
        constructCellMap,
        constructFaceMap,
        constructPointMap
    );
    fvMeshAdder::MapSurfaceFields<sphericalTensor>
    (
        fvMeshes,
        oldFaceOwner,
        mapPtr().oldPatchStarts(),
        mapPtr().oldPatchSizes(),
        patchMap,
        constructCellMap,
        constructFaceMap,
        constructPointMap
    );
    fvMeshAdder::MapSurfaceFields<symmTensor>
    (
        fvMeshes,
        oldFaceOwner,
        mapPtr().oldPatchStarts(),
        mapPtr().oldPatchSizes(),
        patchMap,
        constructCellMap,
        constructFaceMap,
        constructPointMap
    );
    fvMeshAdder::MapSurfaceFields<tensor>
    (
        fvMeshes,
        oldFaceOwner,
        mapPtr().oldPatchStarts(),
        mapPtr().oldPatchSizes(),
        patchMap,
        constructCellMap,
        constructFaceMap,
        constructPointMap
    );
    fvMeshAdder::MapDimFields<scalar>(fvMeshes, constructCellMap);
    fvMeshAdder::MapDimFields<vector>(fvMeshes, constructCellMap);
    fvMeshAdder::MapDimFields<sphericalTensor>(fvMeshes, constructCellMap);
    fvMeshAdder::MapDimFields<symmTensor>(fvMeshes, constructCellMap);
    fvMeshAdder::MapDimFields<tensor>(fvMeshes, constructCellMap);

    // Swap returned data back to processor order
    if (myProci != 0)
    {
        fvMesh* fvm0 = fvMeshes.get(0);
        fvMesh* fvmi = fvMeshes.get(myProci);
        fvMeshes.set(0, fvmi);
        fvMeshes.set(myProci, fvm0);

        // Swap (and renumber) patch face information
        labelListList& lbf = const_cast<labelListList&>(localBoundaryFace);
        std::swap(lbf[0], lbf[myProci]);
        labelListList& rfp = const_cast<labelListList&>(remoteFaceProc);
        std::swap(rfp[0], rfp[myProci]);
        forAll(rfp, proci)
        {
            for (label& proc : rfp[proci])
            {
                if (proc == 0) proc = myProci;
                else if (proc == myProci) proc = 0;
            }
        }
        labelListList& rbf = const_cast<labelListList&>(remoteBoundaryFace);
        std::swap(rbf[0], rbf[myProci]);

        std::swap(constructPatchMap[0], constructPatchMap[myProci]);
        std::swap(constructCellMap[0], constructCellMap[myProci]);
        std::swap(constructFaceMap[0], constructFaceMap[myProci]);
        std::swap(constructPointMap[0], constructPointMap[myProci]);
    }

    return mapPtr;
}


// ************************************************************************* //
