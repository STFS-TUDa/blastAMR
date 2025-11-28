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
        ptr->storeGlobalPositions();
        Info<< "    Cloud '" << c.name() << "' (passiveParticleCloud): "
            << ptr->size() << " particles" << endl;
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

    // Extract positions and sort by destination processor
    for (const auto& p : *ptr)
    {
        const label celli = p.cell();
        if (celli >= 0 && celli < distribution.size())
        {
            posTransferLists[distribution[celli]].append(p.position());
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


// ************************************************************************* //
