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

    // Track last timestep when autoMapClouds was called to prevent double-mapping
    static label lastAutoMapTimeIndex_ = -1;

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

    // Check if we've already done autoMap this timestep (prevents double-mapping
    // when both refinement and load balancing occur in the same timestep)
    const label currentTimeIndex = mesh.time().timeIndex();
    if (lastAutoMapTimeIndex_ == currentTimeIndex)
    {
        if (cloudSupportDebug)
        {
            Info<< "cloudSupport: Skipping duplicate autoMapClouds call "
                << "at timeIndex " << currentTimeIndex << endl;
        }
        return;
    }
    lastAutoMapTimeIndex_ = currentTimeIndex;

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


// Helper struct to store parcel data during transfer
struct ParcelData
{
    point position;
    vector U;           // velocity
    scalar d;           // diameter
    scalar rho;         // density
    scalar nParticle;   // number of particles
    scalar age;         // parcel age
    scalar dTarget;     // target diameter
    label typeId;       // parcel type ID
    bool active;        // active flag
};

// Static storage for pending parcel data during distribution
static HashTable<DynamicList<ParcelData>> pendingParcelDataFull_;

// Helper macro to extract and transfer full parcel data between procs
// This preserves all parcel properties (position, velocity, diameter, etc.)
#define DISTRIBUTE_PARCELS(CloudType, ParcelType)                             \
    if (auto* ptr = dynamic_cast<CloudType*>(&c))                             \
    {                                                                         \
        /* Lists of parcel data to be transferred to each processor */        \
        List<DynamicList<ParcelData>> dataTransferLists(UPstream::nProcs());  \
                                                                              \
        /* Extract parcel data and sort by destination processor */           \
        for (const auto& p : *ptr)                                            \
        {                                                                     \
            const label celli = p.cell();                                     \
            if (celli >= 0 && celli < distribution.size())                    \
            {                                                                 \
                ParcelData pd;                                                \
                pd.position = p.position();                                   \
                pd.U = p.U();                                                 \
                pd.d = p.d();                                                 \
                pd.rho = p.rho();                                             \
                pd.nParticle = p.nParticle();                                 \
                pd.age = p.age();                                             \
                pd.dTarget = p.dTarget();                                     \
                pd.typeId = p.typeId();                                       \
                pd.active = p.active();                                       \
                dataTransferLists[distribution[celli]].append(pd);            \
            }                                                                 \
        }                                                                     \
                                                                              \
        /* Clear the cloud */                                                 \
        ptr->clear();                                                         \
                                                                              \
        if (cloudSupportDebug)                                                \
        {                                                                     \
            forAll(dataTransferLists, procI)                                  \
            {                                                                 \
                Pout<< "      -> proc " << procI << ": "                      \
                    << dataTransferLists[procI].size() << " parcels" << nl;   \
            }                                                                 \
        }                                                                     \
                                                                              \
        /* Stream parcel data into send buffers */                            \
        forAll(dataTransferLists, procI)                                      \
        {                                                                     \
            if (dataTransferLists[procI].size())                              \
            {                                                                 \
                UOPstream os(procI, pBufs);                                   \
                const DynamicList<ParcelData>& list = dataTransferLists[procI]; \
                os << label(list.size());                                     \
                for (const ParcelData& pd : list)                             \
                {                                                             \
                    os << pd.position << pd.U << pd.d << pd.rho               \
                       << pd.nParticle << pd.age << pd.dTarget                \
                       << pd.typeId << pd.active;                             \
                }                                                             \
            }                                                                 \
        }                                                                     \
                                                                              \
        pBufs.finishedSends();                                                \
                                                                              \
        /* Receive parcel data */                                             \
        DynamicList<ParcelData> receivedData;                                 \
        for (const int proci : pBufs.allProcs())                              \
        {                                                                     \
            if (pBufs.recvDataCount(proci))                                   \
            {                                                                 \
                UIPstream is(proci, pBufs);                                   \
                label nParcels;                                               \
                is >> nParcels;                                               \
                for (label i = 0; i < nParcels; i++)                          \
                {                                                             \
                    ParcelData pd;                                            \
                    is >> pd.position >> pd.U >> pd.d >> pd.rho               \
                       >> pd.nParticle >> pd.age >> pd.dTarget                \
                       >> pd.typeId >> pd.active;                             \
                    receivedData.append(pd);                                  \
                }                                                             \
                                                                              \
                if (cloudSupportDebug)                                        \
                {                                                             \
                    Pout<< "      <- proc " << proci << ": "                  \
                        << nParcels << " parcels" << nl;                      \
                }                                                             \
            }                                                                 \
        }                                                                     \
                                                                              \
        /* Store parcel data for later relocation */                          \
        pendingParcelDataFull_.set(cloudName, receivedData);                  \
        globalNewSize = returnReduce(receivedData.size(), sumOp<label>());    \
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


// Helper macro to create parcels from full parcel data for intermediate cloud types
// This restores all parcel properties (position, velocity, diameter, etc.)
#define RELOCATE_KINEMATIC_PARCELS(CloudType, ParcelType)                     \
    if (auto* ptr = dynamic_cast<CloudType*>(&c))                             \
    {                                                                         \
        if (pendingParcelDataFull_.found(cloudName))                          \
        {                                                                     \
            const DynamicList<ParcelData>& dataList =                         \
                pendingParcelDataFull_[cloudName];                            \
            for (const ParcelData& pd : dataList)                             \
            {                                                                 \
                const label celli = mesh.findCell(pd.position);               \
                if (celli >= 0)                                               \
                {                                                             \
                    auto* p = new ParcelType(mesh, pd.position, celli);       \
                    /* Restore all parcel properties */                       \
                    p->U() = pd.U;                                            \
                    p->d() = pd.d;                                            \
                    p->rho() = pd.rho;                                        \
                    p->nParticle() = pd.nParticle;                            \
                    p->age() = pd.age;                                        \
                    p->dTarget() = pd.dTarget;                                \
                    p->typeId() = pd.typeId;                                  \
                    /* Note: active flag handled by constructor (true) */     \
                    ptr->addParticle(p);                                      \
                    nRelocated++;                                             \
                }                                                             \
                else                                                          \
                {                                                             \
                    nLost++;                                                  \
                    if (cloudSupportDebug)                                    \
                    {                                                         \
                        Pout<< "      Lost parcel at " << pd.position         \
                            << " with U=" << pd.U << " (not in mesh)" << nl;  \
                    }                                                         \
                }                                                             \
            }                                                                 \
        }                                                                     \
        handled = true;                                                       \
    }


void Foam::cloudSupport::relocateClouds(const fvMesh& mesh)
{
    // After mesh redistribution, this function creates parcels from
    // the full parcel data that was stored during distributeClouds().

    // Check if we have full parcel data (kinematic clouds) or positions only (passive)
    if (pendingParcelDataFull_.empty() && pendingPositions_.empty())
    {
        if (cloudSupportDebug)
        {
            Info<< "cloudSupport: No pending parcel data to relocate" << endl;
        }
        return;
    }

    // Trigger tet base point calculation on all processors
    // This must be done collectively before locating particles
    (void)mesh.tetBasePtIs();

    Info<< "cloudSupport: Relocating clouds after distribution" << endl;

    // Process each cloud that has pending full parcel data (kinematic clouds)
    forAllIter(HashTable<DynamicList<ParcelData>>, pendingParcelDataFull_, iter)
    {
        const word& cloudName = iter.key();
        const DynamicList<ParcelData>& dataList = iter();

        if (cloudSupportDebug)
        {
            Pout<< "    Processing " << dataList.size()
                << " parcels for cloud '" << cloudName << "'" << endl;
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

    // Process passive particle clouds (positions only)
    forAllIter(HashTable<DynamicList<point>>, pendingPositions_, iter)
    {
        const word& cloudName = iter.key();
        DynamicList<point>& positions = iter();

        // Skip if already handled as kinematic cloud
        if (pendingParcelDataFull_.found(cloudName))
        {
            continue;
        }

        if (cloudSupportDebug)
        {
            Pout<< "    Processing " << positions.size()
                << " positions for passive cloud '" << cloudName << "'" << endl;
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

        // passiveParticleCloud (base case)
        if (auto* ptr = dynamic_cast<Cloud<passiveParticle>*>(&c))
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

    pendingParcelDataFull_.clear();
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
