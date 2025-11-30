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

\*---------------------------------------------------------------------------*/

#include "passiveCloudHandler.H"
#include "addToRunTimeSelectionTable.H"
#include "UOPstream.H"
#include "UIPstream.H"
#include "Cloud.H"
#include "passiveParticle.H"

// Debug switch
extern int cloudSupportDebug;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(passiveCloudHandler, 0);
    addToRunTimeSelectionTable(cloudHandler, passiveCloudHandler, cloud);

    // Static storage for pending positions
    HashTable<DynamicList<point>> passiveCloudHandler::pendingPositions_;

    // Static storage for autoMap positions
    HashTable<DynamicList<point>> passiveCloudHandler::autoMapPositions_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::passiveCloudHandler::canHandle(cloud& c) const
{
    // Check if cloud is a passiveParticleCloud
    return dynamic_cast<Cloud<passiveParticle>*>(&c) != nullptr;
}


void Foam::passiveCloudHandler::storePositions(cloud& c)
{
    if (auto* ptr = dynamic_cast<Cloud<passiveParticle>*>(&c))
    {
        const word& cloudName = c.name();

        // Check if we already have stored positions for this cloud
        // (this prevents overwriting when storePositions is called multiple times
        // before autoMap has a chance to run)
        if (autoMapPositions_.found(cloudName) && autoMapPositions_[cloudName].size() > 0)
        {
            if (cloudSupportDebug)
            {
                Info<< "passiveCloudHandler::storePositions: already have "
                    << autoMapPositions_[cloudName].size()
                    << " positions stored for cloud '" << cloudName
                    << "', skipping" << endl;
            }
            // Still need to ensure OpenFOAM's internal autoMap doesn't fail
            ptr->storeGlobalPositions();
            ptr->clear();
            ptr->storeGlobalPositions();
            return;
        }

        // Store positions for our graceful autoMap handling
        DynamicList<point> positions;
        for (const passiveParticle& p : *ptr)
        {
            positions.append(p.position());
        }

        // Always store entry (even if empty) to ensure all processors take
        // the same code path in autoMap and hit returnReduce collectively
        autoMapPositions_.set(cloudName, positions);

        if (cloudSupportDebug)
        {
            Info<< "passiveCloudHandler::storePositions: stored "
                << positions.size() << " positions for cloud '" << cloudName
                << "', autoMapPositions_ now has " << autoMapPositions_.size()
                << " entries" << endl;
        }

        Info<< "    Cloud '" << c.name() << "' (passiveParticleCloud): "
            << ptr->size() << " particles" << endl;

        // Call storeGlobalPositions to set the internal globalPositionsPtr_
        // This is required by OpenFOAM's internal Cloud::autoMap
        ptr->storeGlobalPositions();

        // Clear the cloud so OpenFOAM's internal Cloud::autoMap (called via
        // regIOobject callback during mesh topology changes) has no particles
        // to process. This works around "Particle mapped to a location outside
        // of the mesh" errors from OpenFOAM's particle.autoMap().
        // Our handler's autoMap will recreate particles from stored positions
        // using mesh.findCell().
        ptr->clear();

        // Also reset globalPositionsPtr_ to work around OpenFOAM's autoMap from
        // failing with "size mismatch" since the cloud is now empty but
        // globalPositionsPtr_ still has the old positions.
        ptr->storeGlobalPositions();
    }
}


void Foam::passiveCloudHandler::distribute
(
    cloud& c,
    const fvMesh& mesh,
    const labelList& distribution,
    PstreamBuffers& pBufs
)
{
    const word& cloudName = c.name();

    auto* ptr = dynamic_cast<Cloud<passiveParticle>*>(&c);
    if (!ptr)
    {
        return;
    }

    // Lists of positions to be transferred to each processor
    List<DynamicList<point>> posTransferLists(UPstream::nProcs());

    // Check if cloud was cleared by storePositions() - if so, use autoMapPositions_
    if (ptr->size() == 0 && autoMapPositions_.found(cloudName))
    {
        // Use stored positions (cloud was cleared for autoMap safety)
        const DynamicList<point>& storedPositions = autoMapPositions_[cloudName];

        if (cloudSupportDebug)
        {
            Pout<< "passiveCloudHandler::distribute: using "
                << storedPositions.size() << " stored positions for cloud '"
                << cloudName << "'" << endl;
        }

        // For each stored position, find its cell and determine destination processor
        for (const point& pos : storedPositions)
        {
            const label celli = mesh.findCell(pos);
            if (celli >= 0 && celli < distribution.size())
            {
                posTransferLists[distribution[celli]].append(pos);
            }
        }

        // Clear autoMapPositions_ since we're distributing now
        autoMapPositions_.erase(cloudName);
    }
    else
    {
        for (const auto& p : *ptr)
        {
            const label celli = p.cell();
            if (celli >= 0 && celli < distribution.size())
            {
                posTransferLists[distribution[celli]].append(p.position());
            }
        }
    }

    // Clear the cloud
    ptr->clear();

    if (cloudSupportDebug)
    {
        forAll(posTransferLists, procI)
        {
            Pout<< "      -> proc " << procI << ": "
                << posTransferLists[procI].size() << " positions" << nl;
        }
    }

    // Stream positions into send buffers
    forAll(posTransferLists, procI)
    {
        if (posTransferLists[procI].size())
        {
            UOPstream os(procI, pBufs);
            os << posTransferLists[procI];
        }
    }

    pBufs.finishedSends();

    // Receive positions
    DynamicList<point> receivedPositions;
    for (const int proci : pBufs.allProcs())
    {
        if (pBufs.recvDataCount(proci))
        {
            UIPstream is(proci, pBufs);
            List<point> positions(is);
            receivedPositions.append(positions);

            if (cloudSupportDebug)
            {
                Pout<< "      <- proc " << proci << ": "
                    << positions.size() << " positions" << nl;
            }
        }
    }

    // Store positions for later relocation
    pendingPositions_.set(cloudName, receivedPositions);
}


void Foam::passiveCloudHandler::relocate(cloud& c, const fvMesh& mesh)
{
    const word& cloudName = c.name();

    if (!pendingPositions_.found(cloudName))
    {
        if (cloudSupportDebug)
        {
            Pout<< "    No pending positions for cloud '" << cloudName << "'" << endl;
        }
        return;
    }

    auto* ptr = dynamic_cast<Cloud<passiveParticle>*>(&c);
    if (!ptr)
    {
        return;
    }

    const DynamicList<point>& positions = pendingPositions_[cloudName];

    label nRelocated = 0;
    label nLost = 0;

    for (const point& pos : positions)
    {
        const label celli = mesh.findCell(pos);
        if (celli >= 0)
        {
            ptr->addParticle(new passiveParticle(mesh, pos, celli));
            nRelocated++;
        }
        else
        {
            nLost++;
        }
    }

    const label globalRelocated = returnReduce(nRelocated, sumOp<label>());
    const label globalLost = returnReduce(nLost, sumOp<label>());

    Info<< "    Cloud '" << cloudName
        << "': relocated " << globalRelocated << " parcels";

    if (globalLost > 0)
    {
        Info<< " (" << globalLost << " lost)";
    }
    Info<< endl;

    // Clear pending data for this cloud
    pendingPositions_.erase(cloudName);
}


Foam::label Foam::passiveCloudHandler::countPerCell
(
    cloud& c,
    labelField& ppCell
) const
{
    label nTotal = 0;

    auto* ptr = dynamic_cast<Cloud<passiveParticle>*>(&c);
    if (ptr)
    {
        for (const auto& p : *ptr)
        {
            const label celli = p.cell();
            if (celli >= 0 && celli < ppCell.size())
            {
                ppCell[celli]++;
                nTotal++;
            }
        }
    }

    return nTotal;
}


void Foam::passiveCloudHandler::updateMesh(cloud& c)
{
    // Passive particle clouds have no injection models or other mesh-dependent
    // data that needs updating after redistribution. This is a no-op.
    if (cloudSupportDebug)
    {
        Pout<< "    Cloud '" << c.name()
            << "': no mesh-dependent data to update (passive)" << endl;
    }
}


void Foam::passiveCloudHandler::autoMap(cloud& c, const mapPolyMesh& map)
{
    const word& cloudName = c.name();

    if (cloudSupportDebug)
    {
        Info<< "passiveCloudHandler::autoMap: cloud '" << cloudName << "'"
            << ", autoMapPositions_ has " << autoMapPositions_.size() << " entries"
            << ", found=" << autoMapPositions_.found(cloudName) << endl;
    }

    if (!autoMapPositions_.found(cloudName))
    {
        if (cloudSupportDebug)
        {
            Info<< "passiveCloudHandler::autoMap: no stored positions for '"
                << cloudName << "', falling back to default" << endl;
        }
        cloudHandler::autoMap(c, map);
        return;
    }

    auto* ptr = dynamic_cast<Cloud<passiveParticle>*>(&c);
    if (!ptr)
    {
        cloudHandler::autoMap(c, map);
        autoMapPositions_.erase(cloudName);
        return;
    }

    const fvMesh& mesh = dynamic_cast<const fvMesh&>(ptr->pMesh());

    // Trigger tet base point calculation on all processors
    (void)mesh.tetBasePtIs();
    (void)mesh.oldCellCentres();

    // Get our stored positions (from before topology change)
    const DynamicList<point>& positions = autoMapPositions_[cloudName];

    // Clear the cloud
    ptr->clear();

    label nRelocated = 0;
    label nLost = 0;

    // Try to relocate each particle using stored positions
    for (const point& pos : positions)
    {
        const label celli = mesh.findCell(pos);
        if (celli >= 0)
        {
            ptr->addParticle(new passiveParticle(mesh, pos, celli));
            nRelocated++;
        }
        else
        {
            nLost++;
        }
    }

    const label globalRelocated = returnReduce(nRelocated, sumOp<label>());
    const label globalLost = returnReduce(nLost, sumOp<label>());

    if (globalLost > 0 || cloudSupportDebug)
    {
        Info<< "    Cloud '" << cloudName
            << "': autoMap relocated " << globalRelocated << " particles";
        if (globalLost > 0)
        {
            Info<< " (" << globalLost << " lost)";
        }
        Info<< endl;
    }

    // Clear stored positions for this cloud
    autoMapPositions_.erase(cloudName);
}


// ************************************************************************* //
