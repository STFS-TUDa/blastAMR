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

#include "cloudSupport.H"
#include "polyMesh.H"
#include "PstreamBuffers.H"
#include "mapDistribute.H"
#include "globalIndex.H"
#include "DynamicList.H"
#include "HashTable.H"
#include "IDLList.H"

// Intermediate cloud includes for dynamic_cast support
#include "basicKinematicCloud.H"
#include "basicKinematicCollidingCloud.H"
#include "basicKinematicMPPICCloud.H"
#include "basicThermoCloud.H"
#include "basicReactingCloud.H"
#include "basicReactingMultiphaseCloud.H"

// Parcel type includes for streaming/creating parcels
#include "basicKinematicParcel.H"
#include "basicKinematicCollidingParcel.H"
#include "basicKinematicMPPICParcel.H"
#include "basicThermoParcel.H"
#include "basicReactingParcel.H"
#include "basicReactingMultiphaseParcel.H"

// Debug switch for detailed particle tracking output
int cloudSupportDebug = Foam::debug::debugSwitch("cloudSupport", 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cloudSupport
{
    // Static storage for pending positions during distribution (legacy)
    static HashTable<DynamicList<point>> pendingPositions_;

    // Static storage for pending parcel data during distribution
    // Key: cloudName, Value: raw binary data containing serialized parcels
    static HashTable<List<char>> pendingParcelData_;

    // Helper macro to try dynamic_cast and call storeGlobalPositions
    #define TRY_STORE_POSITIONS(CloudType)                                    \
        if (auto* ptr = dynamic_cast<CloudType*>(&c))                         \
        {                                                                     \
            ptr->storeGlobalPositions();                                      \
            Info<< "    Cloud '" << c.name() << "' (" << #CloudType           \
                << "): " << ptr->nParcels() << " particles" << endl;          \
            handled = true;                                                   \
        }

    // Helper function to store positions for a single cloud
    bool storePositionsForCloud(cloud& c)
    {
        bool handled = false;

        // Try each known cloud type in order of specificity
        // (most derived first to get best match)

        // Reacting multiphase clouds (most derived)
        TRY_STORE_POSITIONS(basicReactingMultiphaseCloud)
        else TRY_STORE_POSITIONS(basicReactingCloud)
        else TRY_STORE_POSITIONS(basicThermoCloud)
        else TRY_STORE_POSITIONS(basicKinematicCollidingCloud)
        else TRY_STORE_POSITIONS(basicKinematicMPPICCloud)
        else TRY_STORE_POSITIONS(basicKinematicCloud)
        // passiveParticle clouds (base case)
        else if (auto* ptr = dynamic_cast<Cloud<passiveParticle>*>(&c))
        {
            ptr->storeGlobalPositions();
            Info<< "    Cloud '" << c.name() << "' (passiveParticleCloud): "
                << ptr->size() << " particles" << endl;
            handled = true;
        }

        return handled;
    }

    #undef TRY_STORE_POSITIONS
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::cloudSupport::storeGlobalPositions(const fvMesh& mesh)
{
    // Get ALL clouds registered with the mesh using the base cloud class
    UPtrList<const cloud> allClouds = mesh.csorted<cloud>();

    if (allClouds.empty())
    {
        return;
    }

    Info<< "cloudSupport: Storing global positions for "
        << allClouds.size() << " cloud(s)" << endl;

    for (const cloud& constCloud : allClouds)
    {
        cloud& c = const_cast<cloud&>(constCloud);

        if (!storePositionsForCloud(c))
        {
            // Unknown cloud type - warn but don't fail
            WarningInFunction
                << "Unknown cloud type for '" << c.name()
                << "' - positions not stored. "
                << "Load balancing may fail for this cloud." << endl;
        }
    }
}


void Foam::cloudSupport::autoMapClouds
(
    const fvMesh& mesh,
    const mapPolyMesh& map
)
{
    // Note: autoMap is virtual in the base cloud class, so we can call it
    // polymorphically. However, it requires storeGlobalPositions() to have
    // been called first

    UPtrList<const cloud> allClouds = mesh.csorted<cloud>();

    if (allClouds.empty())
    {
        return;
    }

    Info<< "cloudSupport: Remapping " << allClouds.size()
        << " cloud(s) after topology change" << endl;

    for (const cloud& constCloud : allClouds)
    {
        cloud& c = const_cast<cloud&>(constCloud);
        const label oldSize = c.nParcels();

        // autoMap is virtual in base class - call polymorphically
        c.autoMap(map);

        Info<< "    Cloud '" << c.name()
            << "': " << oldSize << " -> " << c.nParcels()
            << " particles" << endl;
    }
}


// Helper macro to extract and transfer parcel positions + data between procs
// This preserves all parcel properties (position, velocity, diameter, etc.)
// Positions are stored separately to allow locating cells on the new mesh.
#define DISTRIBUTE_PARCELS(CloudType, ParcelType)                             \
    if (auto* ptr = dynamic_cast<CloudType*>(&c))                             \
    {                                                                         \
        /* Lists of positions to be transferred to each processor */          \
        List<DynamicList<point>> posTransferLists(UPstream::nProcs());        \
                                                                              \
        /* Extract positions and sort by destination processor */             \
        for (const auto& p : *ptr)                                            \
        {                                                                     \
            const label celli = p.cell();                                     \
            if (celli >= 0 && celli < distribution.size())                    \
            {                                                                 \
                posTransferLists[distribution[celli]].append(p.position());   \
            }                                                                 \
        }                                                                     \
                                                                              \
        /* Clear the cloud */                                                 \
        ptr->clear();                                                         \
                                                                              \
        if (cloudSupportDebug)                                                \
        {                                                                     \
            forAll(posTransferLists, procI)                                   \
            {                                                                 \
                Pout<< "      -> proc " << procI << ": "                      \
                    << posTransferLists[procI].size() << " parcels" << nl;    \
            }                                                                 \
        }                                                                     \
                                                                              \
        /* Stream positions into send buffers */                              \
        forAll(posTransferLists, procI)                                       \
        {                                                                     \
            if (posTransferLists[procI].size())                               \
            {                                                                 \
                UOPstream os(procI, pBufs);                                   \
                os << posTransferLists[procI];                                \
            }                                                                 \
        }                                                                     \
                                                                              \
        pBufs.finishedSends();                                                \
                                                                              \
        /* Receive positions */                                               \
        DynamicList<point> receivedPositions;                                 \
        for (const int proci : pBufs.allProcs())                              \
        {                                                                     \
            if (pBufs.recvDataCount(proci))                                   \
            {                                                                 \
                UIPstream is(proci, pBufs);                                   \
                List<point> positions(is);                                    \
                receivedPositions.append(positions);                          \
                                                                              \
                if (cloudSupportDebug)                                        \
                {                                                             \
                    Pout<< "      <- proc " << proci << ": "                  \
                        << positions.size() << " parcels" << nl;              \
                }                                                             \
            }                                                                 \
        }                                                                     \
                                                                              \
        /* Store positions for later relocation */                            \
        pendingPositions_.set(cloudName, receivedPositions);                  \
        globalNewSize = returnReduce(receivedPositions.size(), sumOp<label>());\
        handled = true;                                                       \
    }

// Helper macro for passiveParticle clouds (position-only transfer)
#define DISTRIBUTE_PASSIVE_POSITIONS(CloudType)                               \
    if (auto* ptr = dynamic_cast<CloudType*>(&c))                             \
    {                                                                         \
        /* Lists of positions to be transferred to each processor */          \
        List<DynamicList<point>> posTransferLists(UPstream::nProcs());        \
                                                                              \
        /* Extract positions and sort by destination processor */             \
        for (const auto& p : *ptr)                                            \
        {                                                                     \
            const label celli = p.cell();                                     \
            if (celli >= 0 && celli < distribution.size())                    \
            {                                                                 \
                posTransferLists[distribution[celli]].append(p.position());   \
            }                                                                 \
        }                                                                     \
                                                                              \
        /* Clear the cloud */                                                 \
        ptr->clear();                                                         \
                                                                              \
        if (cloudSupportDebug)                                                \
        {                                                                     \
            forAll(posTransferLists, procI)                                   \
            {                                                                 \
                Pout<< "      -> proc " << procI << ": "                      \
                    << posTransferLists[procI].size() << " positions" << nl;  \
            }                                                                 \
        }                                                                     \
                                                                              \
        /* Stream positions into send buffers */                              \
        forAll(posTransferLists, procI)                                       \
        {                                                                     \
            if (posTransferLists[procI].size())                               \
            {                                                                 \
                UOPstream os(procI, pBufs);                                   \
                os << posTransferLists[procI];                                \
            }                                                                 \
        }                                                                     \
                                                                              \
        pBufs.finishedSends();                                                \
                                                                              \
        /* Receive positions */                                               \
        DynamicList<point> receivedPositions;                                 \
        for (const int proci : pBufs.allProcs())                              \
        {                                                                     \
            if (pBufs.recvDataCount(proci))                                   \
            {                                                                 \
                UIPstream is(proci, pBufs);                                   \
                List<point> positions(is);                                    \
                receivedPositions.append(positions);                          \
                                                                              \
                if (cloudSupportDebug)                                        \
                {                                                             \
                    Pout<< "      <- proc " << proci << ": "                  \
                        << positions.size() << " positions" << nl;            \
                }                                                             \
            }                                                                 \
        }                                                                     \
                                                                              \
        /* Store positions for later relocation */                            \
        pendingPositions_.set(cloudName, receivedPositions);                  \
        globalNewSize = returnReduce(receivedPositions.size(), sumOp<label>());\
        handled = true;                                                       \
    }


void Foam::cloudSupport::distributeClouds
(
    const fvMesh& mesh,
    const labelList& distribution
)
{
    // For distribution, we transfer full parcel data between processors.
    // For intermediate clouds (kinematicCloud, thermoCloud, etc.), this
    // preserves all parcel properties (velocity, diameter, etc.).
    // For passiveParticleCloud, only positions are transferred.

    UPtrList<const cloud> allClouds = mesh.csorted<cloud>();

    if (allClouds.empty())
    {
        return;
    }

    Info<< "cloudSupport: Distributing " << allClouds.size()
        << " cloud(s) across processors" << endl;

    for (const cloud& constCloud : allClouds)
    {
        cloud& c = const_cast<cloud&>(constCloud);
        const word& cloudName = c.name();
        const label oldSize = c.nParcels();
        const label globalOldSize = returnReduce(oldSize, sumOp<label>());

        if (globalOldSize == 0)
        {
            Info<< "    Cloud '" << cloudName << "': empty, skipping" << endl;
            continue;
        }

        if (cloudSupportDebug)
        {
            Pout<< "    Cloud '" << cloudName << "' on proc "
                << UPstream::myProcNo() << ": " << oldSize << " parcels" << endl;
        }

        // Transfer buffers
        PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);

        bool handled = false;
        label globalNewSize = 0;

        // Try each known cloud type - extract and transfer full parcel data
        // Most derived first to get best match
        DISTRIBUTE_PARCELS(basicReactingMultiphaseCloud, basicReactingMultiphaseParcel)
        else DISTRIBUTE_PARCELS(basicReactingCloud, basicReactingParcel)
        else DISTRIBUTE_PARCELS(basicThermoCloud, basicThermoParcel)
        else DISTRIBUTE_PARCELS(basicKinematicCollidingCloud, basicKinematicCollidingParcel)
        else DISTRIBUTE_PARCELS(basicKinematicMPPICCloud, basicKinematicMPPICParcel)
        else DISTRIBUTE_PARCELS(basicKinematicCloud, basicKinematicParcel)
        // passiveParticleCloud uses position-only transfer
        else DISTRIBUTE_PASSIVE_POSITIONS(Cloud<passiveParticle>)

        if (!handled)
        {
            WarningInFunction
                << "Unknown cloud type for '" << cloudName
                << "' - cannot distribute particles" << endl;
            continue;
        }

        Info<< "    Cloud '" << cloudName
            << "': " << globalOldSize << " -> " << globalNewSize
            << " parcels transferred" << endl;

        if (globalOldSize != globalNewSize)
        {
            WarningInFunction
                << "Parcel count changed during distribution for cloud '"
                << cloudName << "': " << globalOldSize << " -> " << globalNewSize
                << endl;
        }
    }
}

#undef DISTRIBUTE_PARCELS
#undef DISTRIBUTE_PASSIVE_POSITIONS


// Helper macro to create parcels from positions for intermediate cloud types
// Note: This creates parcels at the transferred positions. For kinematic
// parcels, we set nParticle=1 and a small diameter to avoid division by zero.
// Velocity and other properties are initialized to zero/default.
// The simple constructor already sets active_=true.
#define RELOCATE_KINEMATIC_PARCELS(CloudType, ParcelType)                     \
    if (auto* ptr = dynamic_cast<CloudType*>(&c))                             \
    {                                                                         \
        for (const point& pos : positions)                                    \
        {                                                                     \
            const label celli = mesh.findCell(pos);                           \
            if (celli >= 0)                                                   \
            {                                                                 \
                auto* p = new ParcelType(mesh, pos, celli);                   \
                /* Set non-zero values to prevent FPE */                      \
                p->nParticle() = 1.0;                                         \
                p->d() = 1e-6;  /* 1 micron default */                        \
                p->rho() = 1000.0;  /* Water density default */               \
                ptr->addParticle(p);                                          \
                nRelocated++;                                                 \
            }                                                                 \
            else                                                              \
            {                                                                 \
                nLost++;                                                      \
                if (cloudSupportDebug)                                        \
                {                                                             \
                    Pout<< "      Lost position " << pos                      \
                        << " (not in mesh)" << nl;                            \
                }                                                             \
            }                                                                 \
        }                                                                     \
        handled = true;                                                       \
    }


void Foam::cloudSupport::relocateClouds(const fvMesh& mesh)
{
    // After mesh redistribution, this function creates parcels from
    // the positions that were stored during distributeClouds().

    if (pendingPositions_.empty())
    {
        if (cloudSupportDebug)
        {
            Info<< "cloudSupport: No pending positions to relocate" << endl;
        }
        return;
    }

    // Trigger tet base point calculation on all processors
    // This must be done collectively before locating particles
    (void)mesh.tetBasePtIs();

    Info<< "cloudSupport: Relocating clouds after distribution" << endl;

    // Process each cloud that has pending positions
    forAllIter(HashTable<DynamicList<point>>, pendingPositions_, iter)
    {
        const word& cloudName = iter.key();
        DynamicList<point>& positions = iter();

        if (cloudSupportDebug)
        {
            Pout<< "    Processing " << positions.size()
                << " positions for cloud '" << cloudName << "'" << endl;
        }

        // Find the cloud in the registry
        auto* cloudPtr = mesh.getObjectPtr<cloud>(cloudName);
        if (!cloudPtr)
        {
            WarningInFunction
                << "Cloud '" << cloudName << "' not found in mesh registry"
                << endl;
            continue;
        }

        cloud& c = *cloudPtr;
        label nRelocated = 0;
        label nLost = 0;
        bool handled = false;

        // Try intermediate cloud types (most derived first)
        RELOCATE_KINEMATIC_PARCELS(basicReactingMultiphaseCloud, basicReactingMultiphaseParcel)
        else RELOCATE_KINEMATIC_PARCELS(basicReactingCloud, basicReactingParcel)
        else RELOCATE_KINEMATIC_PARCELS(basicThermoCloud, basicThermoParcel)
        else RELOCATE_KINEMATIC_PARCELS(basicKinematicCollidingCloud, basicKinematicCollidingParcel)
        else RELOCATE_KINEMATIC_PARCELS(basicKinematicMPPICCloud, basicKinematicMPPICParcel)
        else RELOCATE_KINEMATIC_PARCELS(basicKinematicCloud, basicKinematicParcel)
        // passiveParticleCloud (base case)
        else if (auto* ptr = dynamic_cast<Cloud<passiveParticle>*>(&c))
        {
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
            handled = true;
        }

        if (!handled)
        {
            WarningInFunction
                << "Unknown cloud type for '" << cloudName
                << "' - cannot create parcels" << endl;
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
    }

    pendingPositions_.clear();
    pendingParcelData_.clear();
}

#undef RELOCATE_KINEMATIC_PARCELS


Foam::label Foam::cloudSupport::countParticles(const fvMesh& mesh)
{
    UPtrList<const cloud> allClouds = mesh.csorted<cloud>();

    label nParticles = 0;

    for (const cloud& c : allClouds)
    {
        nParticles += c.nParcels();
    }

    return nParticles;
}


Foam::tmp<Foam::labelField> Foam::cloudSupport::particlesPerCell
(
    const fvMesh& mesh
)
{
    tmp<labelField> tppCell(new labelField(mesh.nCells(), 0));
    labelField& ppCell = tppCell.ref();

    UPtrList<const cloud> allClouds = mesh.csorted<cloud>();

    for (const cloud& constCloud : allClouds)
    {
        cloud& c = const_cast<cloud&>(constCloud);

        // Need type-specific iteration to access cell()
        #define COUNT_PER_CELL(CloudType)                                     \
            if (auto* ptr = dynamic_cast<CloudType*>(&c))                     \
            {                                                                 \
                for (const auto& p : *ptr)                                    \
                {                                                             \
                    const label celli = p.cell();                             \
                    if (celli >= 0 && celli < ppCell.size())                  \
                    {                                                         \
                        ppCell[celli]++;                                      \
                    }                                                         \
                }                                                             \
            }

        COUNT_PER_CELL(basicReactingMultiphaseCloud)
        else COUNT_PER_CELL(basicReactingCloud)
        else COUNT_PER_CELL(basicThermoCloud)
        else COUNT_PER_CELL(basicKinematicCollidingCloud)
        else COUNT_PER_CELL(basicKinematicMPPICCloud)
        else COUNT_PER_CELL(basicKinematicCloud)
        else COUNT_PER_CELL(Cloud<passiveParticle>)

        #undef COUNT_PER_CELL
    }

    return tppCell;
}


Foam::wordList Foam::cloudSupport::cloudNames(const fvMesh& mesh)
{
    UPtrList<const cloud> allClouds = mesh.csorted<cloud>();

    wordList names(allClouds.size());
    label i = 0;

    for (const cloud& c : allClouds)
    {
        names[i++] = c.name();
    }

    return names;
}


// ************************************************************************* //
