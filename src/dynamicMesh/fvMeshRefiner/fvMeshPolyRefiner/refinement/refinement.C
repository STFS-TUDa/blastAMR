/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "refinement.H"
#include "polyTopoChanger.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "polyModifyFace.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "meshTools.H"
#include "hexRef.H"
#include "dynMeshTools.H"
#include "mapPolyMesh.H"
#include "mapDistributePolyMesh.H"
#include "globalIndex.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(refinement, 0);

    // Note: do not add to run-time selection table since this is abstract base
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::refinement::setInstance(const fileName& inst) const
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << "Resetting file instance of refinement data to " << inst
            << endl;
    }

    cellLevel_.instance() = inst;
    pointLevel_.instance() = inst;
    parentCells_.instance() = inst;
}


void Foam::refinement::walkFaceToMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label faceI,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    // Get the face and its edges
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges()[faceI];

    label fp = startFp;

    // Starting from fp store all (1 or 2) vertices until where the face
    // gets split
    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] > -1)
        {
            // Edge is split, append its mid point
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (pointLevel_[f[fp]] <= cLevel)
        {
            // Found next anchor. Already appended the split point just above
            return;
        }
        else if (pointLevel_[f[fp]] == cLevel + 1)
        {
            // Mid level, append and return
            faceVerts.append(f[fp]);

            return;
        }
        else if (pointLevel_[f[fp]] == cLevel + 2)
        {
            // Store and continue to cLevel + 1
            faceVerts.append(f[fp]);
        }
    }
}


void Foam::refinement::walkFaceFromMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label faceI,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    // Get the face and its edges
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges()[faceI];

    label fp = f.rcIndex(startFp);

    while (true)
    {
        if (pointLevel_[f[fp]] <= cLevel)
        {
            // Anchor point, break
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel + 1)
        {
            // Mid level, append and break
            faceVerts.append(f[fp]);
            break;
        }
        // Continue to next point
        fp = f.rcIndex(fp);
    }

    // Store
    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] > -1)
        {
            // Edge is split, append its mid point
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (fp == startFp)
        {
            break;
        }
        faceVerts.append(f[fp]);
    }
}


Foam::label Foam::refinement::findMinLevel(const labelList& f) const
{
    return findMin(labelField(pointLevel_, f));
}


Foam::label Foam::refinement::findMaxLevel(const labelList& f) const
{
    return findMax(labelField(pointLevel_, f));
}


Foam::label Foam::refinement::countAnchors
(
    const labelList& f,
    const label anchorLevel
) const
{
    label nAnchors = 0;

    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= anchorLevel)
        {
            nAnchors++;
        }
    }
    return nAnchors;
}


Foam::label Foam::refinement::faceConsistentRefinement
(
    boolList& refineCell,
    const bool maxSet
) const
{
    // Count number of cells that will be changed
    label nChanged = 0;

    // Get necessary mesh data
    const label nFaces = mesh_.nFaces();
    const label nInternalFaces = mesh_.nInternalFaces();

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    // Loop through internal faces and check consistency
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        // Get owner and neighbour labels
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Get owner and neighbour cell levels
        // Note: If the cell is marked for refinement, the level is current
        // level + 1, otherwise it is equal to the current level
        const label ownLevel = cellLevel_[own] + refineCell.get(own);
        const label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        if (ownLevel > (neiLevel + 1))
        {
            if (maxSet)
            {
                nChanged += refineCell.set(nei);
            }
            else
            {
                nChanged += refineCell.unset(own);
            }
        }
        else if (neiLevel > (ownLevel + 1))
        {
            if (maxSet)
            {
                nChanged += refineCell.set(own);
            }
            else
            {
                nChanged += refineCell.unset(nei);
            }
        }
    }

    // Create owner level for boundary faces to prepare for swapping on coupled
    // boundaries
    labelList neiLevel(nFaces - nInternalFaces);
    forAll (neiLevel, i)
    {
        // Get owner of the face and update owner cell levels
        const label& own = owner[i + nInternalFaces];
        neiLevel[i] = cellLevel_[own] + refineCell.get(own);
    }

    // Swap boundary face lists (coupled boundary update)
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Loop through boundary faces
    forAll (neiLevel, i)
    {
        // Get owner of the face and owner level
        const label& own = owner[i + nInternalFaces];
        const label ownLevel = cellLevel_[own] + refineCell.get(own);

        // Note: we are using more stringent 1:1 consistency across coupled
        // boundaries in order to simplify handling of edge based consistency
        // checks for parallel runs
        // Bugfix related to PLB: Check whether owner is already marked for
        // refinement. Will allow 2:1 consistency across certain processor
        // faces where we have a new processor boundary. VV, 23/Jan/2019.
        if (ownLevel > (neiLevel[i] + 1))
        {
            if (!maxSet)
            {
                nChanged += refineCell.unset(own);
            }
        }
        else if (neiLevel[i] > (ownLevel + 1))
        {
            if (maxSet)
            {
                nChanged += refineCell.set(own);
            }
        }

        // Note: other possibility (that owner level is higher than neighbour
        // level) is taken into account on the other side automatically
    }

    // Return number of changed cells
    return nChanged;
}


Foam::label Foam::refinement::edgeConsistentRefinement
(
    boolList& refineCell,
    const bool maxSet
) const
{
    // Count number of cells that will be added
    label nChanged = 0;

    // Algorithm: loop over all edges and visit all unique cell pairs sharing
    // this particular edge. Then, ensure 2:1 edge consistency by marking
    // cell with lower level for refinement

    // Get edge cells
    const labelListList& meshEdgeCells = mesh_.edgeCells();

    // Loop through all mesh edges
    forAll (meshEdgeCells, edgeI)
    {
        // Get current edge cells
        const labelList& curEdgeCells = meshEdgeCells[edgeI];

        // Loop through all edge cells
        forAll (curEdgeCells, i)
        {
            // Get first cell index
            const label& cellI = curEdgeCells[i];

            // Loop through remaining edge cells
            for (label j = i + 1; j < curEdgeCells.size(); ++j)
            {
                // Get second cell index
                const label& cellJ = curEdgeCells[j];

                // Get levels of the two cells. If the cell is marked for
                // refinement, the level is current level + 1, otherwise it is
                // equal to the current level

                // Note: refineCell flag for both cellI and cellJ might
                // change, this is why we need to recalculate cellI level here
                const label cellILevel =
                    cellLevel_[cellI] + int(refineCell.get(cellI));
                const label cellJLevel =
                    cellLevel_[cellJ] + int(refineCell.get(cellJ));

                if (cellILevel > (cellJLevel + 1))
                {
                    if (maxSet)
                    {
                        nChanged += refineCell.set(cellJ);
                    }
                    else
                    {
                        nChanged += refineCell.unset(cellI);
                    }
                }
                else if (cellJLevel > (cellILevel + 1))
                {
                    if (maxSet)
                    {
                        nChanged += refineCell.set(cellI);
                    }
                    else
                    {
                        nChanged += refineCell.unset(cellJ);
                    }
                }
            }
        }
    }

    // Note: in order to avoid very difficult and time-consuming parallelisation
    // of edge cell connectivity and edge cell values, we enforce a more
    // stringent face-based consistency across processor boundaries. Basically,
    // if a face-based consistency of 1:1 (not 2:1 as for ordinary faces) is
    // ensured, the edge-based consistency becomes a local operation (I'm not
    // 100% sure to be honest since there are countless variants when dealing
    // with arbitrary polyhedral cells).
    // See faceConsistentRefinement for details. VV, 17/Apr/2018.

    // Return number of changed cells
    return nChanged;
}


Foam::label Foam::refinement::faceConsistentUnrefinement
(
    boolList& unrefineCell,
    const bool maxSet
) const
{
    // Count number of removed cells from unrefinement
    label nChanged = 0;

    // Get necessary mesh data
    const label nFaces = mesh_.nFaces();
    const label nInternalFaces = mesh_.nInternalFaces();

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    // Loop through internal faces and check consistency
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        // Get owner and neighbour labels
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Get owner and neighbour cell levels
        // Note: If the cell is marked for unrefinement, the level is current
        // level - 1, otherwise it is equal to the current level
        const label ownLevel = cellLevel_[own] - unrefineCell.get(own);
        const label neiLevel = cellLevel_[nei] - unrefineCell.get(nei);

        if (ownLevel < (neiLevel - 1))
        {
            // Since was 2:1 this can only occur if own is marked for
            // unrefinement.

            if (maxSet)
            {
                nChanged += unrefineCell.set(nei);
            }
            else
            {
                if (!unrefineCell.get(own))
                {
                    FatalErrorInFunction
                        << "problem" << abort(FatalError);
                }

                nChanged += unrefineCell.unset(own);
            }
        }
        else if (neiLevel < (ownLevel - 1))
        {
            if (maxSet)
            {
                nChanged += unrefineCell.set(own);
            }
            else
            {
                if (!unrefineCell.get(nei))
                {
                    FatalErrorInFunction
                        << "problem" << abort(FatalError);
                }

                nChanged += unrefineCell.unset(nei);
            }
        }
    }

    // Create owner level for boundary faces to prepare for swapping on coupled
    // boundaries
    labelList neiLevel(nFaces - nInternalFaces);
    forAll (neiLevel, i)
    {
        // Get owner of the face and update owner cell levels
        const label& own = owner[i + nInternalFaces];
        neiLevel[i] = cellLevel_[own] - unrefineCell.get(own);
    }

    // Swap boundary face lists (coupled boundary update)
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Loop through boundary faces
    forAll (neiLevel, i)
    {
        // Get owner of the face and owner level
        const label& own = owner[i + nInternalFaces];
        const label ownLevel = cellLevel_[own] - unrefineCell.get(own);

        if (ownLevel < (neiLevel[i] - 1))
        {
            if (!maxSet)
            {
                if (!unrefineCell.get(own))
                {
                    FatalErrorInFunction
                        << "problem" << abort(FatalError);
                }

                nChanged += unrefineCell.unset(own);
            }
        }
        else if (neiLevel[i] < (ownLevel - 1))
        {
            if (maxSet)
            {
                if (!unrefineCell.get(own))
                {
                    FatalErrorInFunction
                        << "problem" << abort(FatalError);
                }

                nChanged += unrefineCell.set(own);
            }
        }
    }

    // Return number of local changed
    return nChanged;
}


Foam::label Foam::refinement::edgeConsistentUnrefinement
(
    boolList& unrefineCell,
    const bool maxSet
) const
{
    // Count number of cells that will be removed
    label nChanged = 0;

    // Algorithm: loop over all edges and visit all unique cell pairs sharing
    // this particular edge. Then, ensure 2:1 edge consistency by protecting the
    // cell with lower level from unrefinement

    // Get edge cells
    const labelListList& meshEdgeCells = mesh_.edgeCells();

    // Loop through all mesh edges
    forAll (meshEdgeCells, edgeI)
    {
        // Get current edge cells
        const labelList& curEdgeCells = meshEdgeCells[edgeI];

        // Loop through all edge cells
        forAll (curEdgeCells, i)
        {
            // Get first cell index
            const label& cellI = curEdgeCells[i];

            // Loop through remaining edge cells
            for (label j = i + 1; j < curEdgeCells.size(); ++j)
            {
                // Get second cell index
                const label& cellJ = curEdgeCells[j];

                // Get levels of the two cells. If the cell is marked for
                // unrefinement, the level is current level - 1, otherwise it is
                // equal to the current level

                // Note: unrefineCell flag for both cellI and cellJ might
                // change, this is why we need to recalculate cellI level here
                const label cellILevel =
                    cellLevel_[cellI] - unrefineCell.get(cellI);
                const label cellJLevel =
                    cellLevel_[cellJ] - unrefineCell.get(cellJ);

                if (cellILevel < (cellJLevel - 1))
                {
                    if (maxSet)
                    {
                        nChanged += unrefineCell.set(cellJ);
                    }
                    else
                    {
                        if (!unrefineCell.get(cellI))
                        {
                            FatalErrorInFunction
                                << "problem" << abort(FatalError);
                        }

                        nChanged += unrefineCell.unset(cellI);
                    }
                }
                else if (cellJLevel < (cellILevel - 1))
                {
                    if (maxSet)
                    {
                        nChanged += unrefineCell.set(cellI);
                    }
                    else
                    {
                        if (!unrefineCell.get(cellJ))
                        {
                            FatalErrorInFunction
                                << "problem" << abort(FatalError);
                        }
                        nChanged += unrefineCell.unset(cellJ);
                    }
                }
            }
        }
    }

    // Note: in order to avoid very difficult and time-consuming parallelisation
    // of edge cell connectivity and edge cell values, we enforce a more
    // stringent face-based consistency across processor boundaries. Basically,
    // if a face-based consistency of 1:1 (not 2:1 as for ordinary faces) is
    // ensured, the edge-based consistency becomes a local operation (I'm not
    // 100% sure to be honest whether this is true all the time since there are
    // countless variants when dealing with arbitrary polyhedral cells).
    // See faceConsistentRefinement for details. VV, 3/Apr/2018.

    // Return number of changed cells
    return nChanged;
}


Foam::label Foam::refinement::getCellClusters(labelList& clusters) const
{
    // Count up the number of clusters that have been refined atleast once
    Map<label> count;
    forAll(parentCells_, ci)
    {
        const label celli = parentCells_[ci];
        if (parentCells_[ci] == -1)
        {}
        else if (!count.found(celli))
        {
            count.insert(celli, 1);
        }
        else
        {
            count[celli]++;
        }
    }

    // Create a map from the parent index to the cluster number
    Map<label> map;
    label i = 0;
    forAllConstIter(Map<label>, count, iter)
    {
        if (iter() > 1)
        {
            map.insert(iter.key(), i++);
        }
    }

    // Return the list of clusters and their children cells
    clusters.setSize(mesh_.nCells(), -1);
    forAll(parentCells_, celli)
    {
        if (map.found(parentCells_[celli]))
        {
            clusters[celli] = map[parentCells_[celli]];
        }
    }
    return i;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinement::refinement
(
    const polyMesh& mesh,
    const dictionary& dict,
    const bool read
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            read ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh_.nCells(), 0)
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            read ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh_.nPoints(), 0)
    ),
    parentCells_
    (
        IOobject
        (
            "parentCells",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            read ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        identity(mesh_.nCells())
    ),
    faceRemover_(mesh_, GREAT),   // Merge boundary faces wherever possible
    edgeBasedConsistency_
    (
        dict.lookupOrDefault<Switch>("edgeBasedConsistency", true)
    ),
    locationMapper_(mesh_)
{
    if (!parentCells_.headerOk())
    {
        globalIndex gI(mesh_.nCells());
        forAll(parentCells_, celli)
        {
            parentCells_[celli] = gI.toGlobal(parentCells_[celli]);
        }
    }
    DebugInfo<< "Created pointLevel and cellLevel" << endl;

    // Check consistency between cellLevel and number of cells and pointLevel
    // and number of points in the mesh
    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "Restarted from inconsistent cellLevel or pointLevel files."
            << endl
            << "Number of cells in mesh: " << mesh_.nCells()
            << " does not equal size of cellLevel: " << cellLevel_.size() << nl
            << "Number of points in mesh: " << mesh_.nPoints()
            << " does not equal size of pointLevel: " << pointLevel_.size()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::refinement::~refinement()
{}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

Foam::labelList Foam::refinement::consistentRefinement
(
    const labelList& refinementCellCandidates,
    const bool maxSet
) const
{
    if (debug)
    {
        InfoInFunction<< "Setting cells to refine" << endl;
    }

    // Create a mark-up field for cells to refine
    boolList refineCell(mesh_.nCells(), false);
    forAll(refinementCellCandidates, i)
    {
        refineCell.set(refinementCellCandidates[i], true);
    }

    // Make sure that the refinement is face consistent (2:1 consistency) and
    // point consistent (4:1 consistency) if necessary

    // Counter for additional cells to refine due to consistency in each
    // iteration and number of iterations
    label nChanged = 0;

    while (true)
    {
        // Check for 2:1 face based consistent refinement. Updates cellsToRefine
        // and returns number of cells added in this iteration
        nChanged = 0;

        if (edgeBasedConsistency_)
        {
            // Check for 4:1 edge based consistent refinement. Updates
            // cellsToRefine and returns number of cells added in this iteration
            nChanged += edgeConsistentRefinement(refineCell, maxSet);
        }
        nChanged += faceConsistentRefinement(refineCell, maxSet);;

        // Global reduction
        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }

    }

    // Convert back to labelList.
    label nRefined = std::count(refineCell.begin(), refineCell.end(), true);

    // Collect all cells to refine in a dynamic list
    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll (refineCell, cellI)
    {
        if (refineCell.get(cellI))
        {
            // Cell marked for refinement, append it
            newCellsToRefine[nRefined++] = cellI;
        }
    }

    return newCellsToRefine;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::refinement::refine
(
    polyMesh& mesh,
    const labelList& cellsToRefine
)
{
    polyTopoChange meshMod(mesh);
    this->setRefinement(meshMod, cellsToRefine);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
    mesh.updateMesh(map());

    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << mesh_.globalData().nTotalCells() << " cells." << endl;
    return map;
}


bool Foam::refinement::unrefine
(
    polyMesh& mesh,
    const labelList& splitPointsToUnrefine
)
{
    polyTopoChange meshMod(mesh);
    this->setUnrefinement(meshMod, splitPointsToUnrefine);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
    mesh.updateMesh(map());

    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << mesh_.globalData().nTotalCells() << " cells."
        << endl;

    return true;
}

void Foam::refinement::updateMesh(const mapPolyMesh& map)
{
    if (debug)
    {
        InfoInFunction
            << "Updating cell and point levels."
            << endl;

        Pout<< "hexRef::updateMesh :"
            << " reverseCellMap:" << map.reverseCellMap().size()
            << " cellMap:" << map.cellMap().size()
            << " nCells:" << mesh_.nCells()
            << " nOldCells:" << map.nOldCells()
            << " cellLevel_:" << cellLevel_.size()
            << " reversePointMap:" << map.reversePointMap().size()
            << " pointMap:" << map.pointMap().size()
            << " nPoints:" << mesh_.nPoints()
            << " nOldPoints:" << map.nOldPoints()
            << " pointLevel_:" << pointLevel_.size()
            << endl;
    }

    // Mesh has changed topologically, we need to update cell and point levels
    // and optionally face removal object
    {
        const labelList& reverseCellMap = map.reverseCellMap();
        if (reverseCellMap.size() == cellLevel_.size())
        {
            // Assume it is after hexRef that this routine is called.
            // Just account for reordering. We cannot use cellMap since
            // then cells created from cells would get cellLevel_ of
            // cell they were created from.
            hexRef::reorder(reverseCellMap, mesh_.nCells(), -1, cellLevel_);
        }
        else
        {
            // Map data
            const labelList& cellMap = map.cellMap();

            label newParentIndex = gMax(parentCells_)+1;

            labelList newCellLevel(cellMap.size());
            labelList newParentCells(cellMap.size(), -1);
            forAll(cellMap, newCelli)
            {
                label oldCelli = cellMap[newCelli];

                if (oldCelli == -1)
                {
                    newCellLevel[newCelli] = -1;
                    newParentCells[newCelli] = newParentIndex++;
                }
                else
                {
                    newCellLevel[newCelli] = cellLevel_[oldCelli];
                    newParentCells[newCelli] = parentCells_[oldCelli];
                }
            }
            cellLevel_.transfer(newCellLevel);
            parentCells_.transfer(newParentCells);
        }

        // Update the parent cells
        {
            // Map data
            const labelList& cellMap = map.cellMap();

            label newParentIndex = gMax(parentCells_)+1;

            labelList newParentCells(cellMap.size(), -1);
            forAll(cellMap, newCelli)
            {
                label oldCelli = cellMap[newCelli];

                if (oldCelli == -1)
                {
                    newParentCells[newCelli] = newParentIndex++;
                }
                else
                {
                    newParentCells[newCelli] = parentCells_[oldCelli];
                }
            }
            parentCells_.transfer(newParentCells);
        }

        const labelList& reversePointMap = map.reversePointMap();
        if (reversePointMap.size() == pointLevel_.size())
        {
            // Assume it is after hexRef that this routine is called.
            hexRef::reorder(reversePointMap, mesh_.nPoints(), -1,  pointLevel_);
        }
        else
        {
            // Map data
            const labelList& pointMap = map.pointMap();

            labelList newPointLevel(pointMap.size());

            forAll(pointMap, newPointi)
            {
                label oldPointi = pointMap[newPointi];

                if (oldPointi == -1)
                {
                    FatalErrorInFunction
                       << "Problem : point " << newPointi
                       << " at " << mesh_.points()[newPointi]
                       << " does not originate from another point"
                       << " (i.e. is inflated)." << nl
                       << "Hence we cannot determine the new pointLevel"
                       << " for it." << abort(FatalError);
                    newPointLevel[newPointi] = -1;
                }
                else
                {
                    newPointLevel[newPointi] = pointLevel_[oldPointi];
                }
            }
            pointLevel_.transfer(newPointLevel);
        }
    }

    // Update face remover
    faceRemover_.updateMesh(map);

    // Mark files as changed
    setInstance(mesh_.facesInstance());
}


void Foam::refinement::distribute(const mapDistributePolyMesh& map)
{
    if (debug)
    {
        InfoInFunction
            << "Distributing cell and point levels."
            << endl;
    }

    // Update cell level
    map.distributeCellData(cellLevel_);

    // Update pointlevel
    map.distributePointData(pointLevel_);

    // Distribute the parent cells;
    map.distributeCellData(parentCells_);

    // Update face removal engine
    faceRemover_.distribute(map);

    // Mark files as changed
    setInstance(mesh_.facesInstance());
}


void Foam::refinement::add
(
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    blockedFace.setSize(mesh_.nFaces(), true);

    // Find common parent for all cells
    labelList cellToCluster;
    getCellClusters(cellToCluster);

    // Unblock all faces inbetween same cluster
    label nUnblocked = 0;

    forAll(mesh_.faceNeighbour(), faceI)
    {
        label ownCluster = cellToCluster[mesh_.faceOwner()[faceI]];
        label neiCluster = cellToCluster[mesh_.faceNeighbour()[faceI]];

        if (ownCluster != -1 && ownCluster == neiCluster)
        {
            if (blockedFace[faceI])
            {
                blockedFace[faceI] = false;
                nUnblocked++;
            }
        }
    }

    if (debug)
    {
        reduce(nUnblocked, sumOp<label>());
        Info<< type() << " : unblocked " << nUnblocked << " faces" << endl;
    }

    syncTools::syncFaceList(mesh_, blockedFace, andEqOp<bool>());
}


void Foam::refinement::apply
(
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    // Find common parent for all cells
    labelList cellToCluster;
    label nClusters = getCellClusters(cellToCluster);

    // Unblock all faces inbetween same cluster
    labelList clusterToProc(nClusters, -1);

    label nChanged = 0;

    forAll(mesh_.faceNeighbour(), faceI)
    {
        label own = mesh_.faceOwner()[faceI];
        label nei = mesh_.faceNeighbour()[faceI];

        label ownCluster = cellToCluster[own];
        label neiCluster = cellToCluster[nei];

        if (ownCluster != -1 && ownCluster == neiCluster)
        {
            if (clusterToProc[ownCluster] == -1)
            {
                clusterToProc[ownCluster] = decomposition[own];
            }

            if (decomposition[own] != clusterToProc[ownCluster])
            {
                decomposition[own] = clusterToProc[ownCluster];
                nChanged++;
            }
            if (decomposition[nei] != clusterToProc[ownCluster])
            {
                decomposition[nei] = clusterToProc[ownCluster];
                nChanged++;
            }
        }
    }

    if (debug)
    {
        reduce(nChanged, sumOp<label>());
        Info<< typeName << ": changed decomposition on " << nChanged
            << " cells" << endl;
    }
}


bool Foam::refinement::write() const
{
    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Write necessary data before writing dictionary
    return cellLevel_.write() && pointLevel_.write();
}


// ************************************************************************* //
