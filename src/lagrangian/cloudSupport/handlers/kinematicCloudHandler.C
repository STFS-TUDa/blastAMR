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

#include "kinematicCloudHandler.H"
#include "addToRunTimeSelectionTable.H"
#include "UOPstream.H"
#include "UIPstream.H"

// Intermediate cloud includes
#include "basicKinematicCloud.H"
#include "basicKinematicCollidingCloud.H"
#include "basicKinematicMPPICCloud.H"
#include "basicThermoCloud.H"
#include "basicReactingCloud.H"
#include "basicReactingMultiphaseCloud.H"

// Parcel type includes
#include "basicKinematicParcel.H"
#include "basicKinematicCollidingParcel.H"
#include "basicKinematicMPPICParcel.H"
#include "basicThermoParcel.H"
#include "basicReactingParcel.H"
#include "basicReactingMultiphaseParcel.H"

// Debug switch
extern int cloudSupportDebug;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kinematicCloudHandler, 0);
    addToRunTimeSelectionTable(cloudHandler, kinematicCloudHandler, cloud);

    // Static storage for pending parcel data
    HashTable<DynamicList<kinematicCloudHandler::ParcelData>>
        kinematicCloudHandler::pendingParcelData_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::kinematicCloudHandler::canHandle(cloud& c) const
{
    // Check if cloud is any of the supported kinematic cloud types
    // Note: Each cloud type is a distinct template instantiation, so we must
    // check each one explicitly. They don't share a common base class that
    // we can cast to.
    return
        dynamic_cast<basicReactingMultiphaseCloud*>(&c) != nullptr
     || dynamic_cast<basicReactingCloud*>(&c) != nullptr
     || dynamic_cast<basicThermoCloud*>(&c) != nullptr
     || dynamic_cast<basicKinematicCollidingCloud*>(&c) != nullptr
     || dynamic_cast<basicKinematicMPPICCloud*>(&c) != nullptr
     || dynamic_cast<basicKinematicCloud*>(&c) != nullptr;
}


void Foam::kinematicCloudHandler::storePositions(cloud& c)
{
    // Use dynamic_cast chain from most derived to least derived
    // to get accurate particle count reporting

    #define TRY_STORE(CloudType)                                              \
        if (auto* ptr = dynamic_cast<CloudType*>(&c))                         \
        {                                                                     \
            ptr->storeGlobalPositions();                                      \
            Info<< "    Cloud '" << c.name() << "' (" << #CloudType           \
                << "): " << ptr->nParcels() << " particles" << endl;          \
            return;                                                           \
        }

    TRY_STORE(basicReactingMultiphaseCloud)
    else TRY_STORE(basicReactingCloud)
    else TRY_STORE(basicThermoCloud)
    else TRY_STORE(basicKinematicCollidingCloud)
    else TRY_STORE(basicKinematicMPPICCloud)
    else TRY_STORE(basicKinematicCloud)

    #undef TRY_STORE
}


void Foam::kinematicCloudHandler::distribute
(
    cloud& c,
    const fvMesh& mesh,
    const labelList& distribution,
    PstreamBuffers& pBufs
)
{
    const word& cloudName = c.name();

    // Helper macro to extract and transfer full parcel data between procs
    #define DISTRIBUTE_CLOUD(CloudType)                                       \
        if (auto* ptr = dynamic_cast<CloudType*>(&c))                         \
        {                                                                     \
            /* Lists of parcel data to be transferred to each processor */    \
            List<DynamicList<ParcelData>> dataTransferLists(UPstream::nProcs()); \
                                                                              \
            /* Extract parcel data and sort by destination processor */       \
            for (const auto& p : *ptr)                                        \
            {                                                                 \
                const label celli = p.cell();                                 \
                if (celli >= 0 && celli < distribution.size())                \
                {                                                             \
                    ParcelData pd;                                            \
                    pd.position = p.position();                               \
                    pd.U = p.U();                                             \
                    pd.d = p.d();                                             \
                    pd.rho = p.rho();                                         \
                    pd.nParticle = p.nParticle();                             \
                    pd.age = p.age();                                         \
                    pd.dTarget = p.dTarget();                                 \
                    pd.typeId = p.typeId();                                   \
                    pd.active = p.active();                                   \
                    dataTransferLists[distribution[celli]].append(pd);        \
                }                                                             \
            }                                                                 \
                                                                              \
            /* Clear the cloud */                                             \
            ptr->clear();                                                     \
                                                                              \
            if (cloudSupportDebug)                                            \
            {                                                                 \
                forAll(dataTransferLists, procI)                              \
                {                                                             \
                    Pout<< "      -> proc " << procI << ": "                  \
                        << dataTransferLists[procI].size() << " parcels" << nl; \
                }                                                             \
            }                                                                 \
                                                                              \
            /* Stream parcel data into send buffers */                        \
            forAll(dataTransferLists, procI)                                  \
            {                                                                 \
                if (dataTransferLists[procI].size())                          \
                {                                                             \
                    UOPstream os(procI, pBufs);                               \
                    const DynamicList<ParcelData>& list = dataTransferLists[procI]; \
                    os << label(list.size());                                 \
                    for (const ParcelData& pd : list)                         \
                    {                                                         \
                        os << pd.position << pd.U << pd.d << pd.rho           \
                           << pd.nParticle << pd.age << pd.dTarget            \
                           << pd.typeId << pd.active;                         \
                    }                                                         \
                }                                                             \
            }                                                                 \
                                                                              \
            pBufs.finishedSends();                                            \
                                                                              \
            /* Receive parcel data */                                         \
            DynamicList<ParcelData> receivedData;                             \
            for (const int proci : pBufs.allProcs())                          \
            {                                                                 \
                if (pBufs.recvDataCount(proci))                               \
                {                                                             \
                    UIPstream is(proci, pBufs);                               \
                    label nParcels;                                           \
                    is >> nParcels;                                           \
                    for (label i = 0; i < nParcels; i++)                      \
                    {                                                         \
                        ParcelData pd;                                        \
                        is >> pd.position >> pd.U >> pd.d >> pd.rho           \
                           >> pd.nParticle >> pd.age >> pd.dTarget            \
                           >> pd.typeId >> pd.active;                         \
                        receivedData.append(pd);                              \
                    }                                                         \
                                                                              \
                    if (cloudSupportDebug)                                    \
                    {                                                         \
                        Pout<< "      <- proc " << proci << ": "              \
                            << nParcels << " parcels" << nl;                  \
                    }                                                         \
                }                                                             \
            }                                                                 \
                                                                              \
            /* Store parcel data for later relocation */                      \
            pendingParcelData_.set(cloudName, receivedData);                  \
            return;                                                           \
        }

    // Try each cloud type (most derived first)
    DISTRIBUTE_CLOUD(basicReactingMultiphaseCloud)
    else DISTRIBUTE_CLOUD(basicReactingCloud)
    else DISTRIBUTE_CLOUD(basicThermoCloud)
    else DISTRIBUTE_CLOUD(basicKinematicCollidingCloud)
    else DISTRIBUTE_CLOUD(basicKinematicMPPICCloud)
    else DISTRIBUTE_CLOUD(basicKinematicCloud)

    #undef DISTRIBUTE_CLOUD
}


void Foam::kinematicCloudHandler::relocate(cloud& c, const fvMesh& mesh)
{
    const word& cloudName = c.name();

    if (!pendingParcelData_.found(cloudName))
    {
        if (cloudSupportDebug)
        {
            Pout<< "    No pending data for cloud '" << cloudName << "'" << endl;
        }
        return;
    }

    const DynamicList<ParcelData>& dataList = pendingParcelData_[cloudName];

    // Helper macro to create parcels from stored data
    #define RELOCATE_CLOUD(CloudType, ParcelType)                             \
        if (auto* ptr = dynamic_cast<CloudType*>(&c))                         \
        {                                                                     \
            label nRelocated = 0;                                             \
            label nLost = 0;                                                  \
                                                                              \
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
                                                                              \
            const label globalRelocated = returnReduce(nRelocated, sumOp<label>()); \
            const label globalLost = returnReduce(nLost, sumOp<label>());     \
                                                                              \
            Info<< "    Cloud '" << cloudName                                 \
                << "': relocated " << globalRelocated << " parcels";          \
                                                                              \
            if (globalLost > 0)                                               \
            {                                                                 \
                Info<< " (" << globalLost << " lost)";                        \
            }                                                                 \
            Info<< endl;                                                      \
                                                                              \
            /* Clear pending data for this cloud */                           \
            pendingParcelData_.erase(cloudName);                              \
            return;                                                           \
        }

    // Try each cloud type (most derived first)
    RELOCATE_CLOUD(basicReactingMultiphaseCloud, basicReactingMultiphaseParcel)
    else RELOCATE_CLOUD(basicReactingCloud, basicReactingParcel)
    else RELOCATE_CLOUD(basicThermoCloud, basicThermoParcel)
    else RELOCATE_CLOUD(basicKinematicCollidingCloud, basicKinematicCollidingParcel)
    else RELOCATE_CLOUD(basicKinematicMPPICCloud, basicKinematicMPPICParcel)
    else RELOCATE_CLOUD(basicKinematicCloud, basicKinematicParcel)

    #undef RELOCATE_CLOUD
}


Foam::label Foam::kinematicCloudHandler::countPerCell
(
    cloud& c,
    labelField& ppCell
) const
{
    label nTotal = 0;

    // Helper macro to count particles per cell
    #define COUNT_CLOUD(CloudType)                                            \
        if (auto* ptr = dynamic_cast<CloudType*>(&c))                         \
        {                                                                     \
            for (const auto& p : *ptr)                                        \
            {                                                                 \
                const label celli = p.cell();                                 \
                if (celli >= 0 && celli < ppCell.size())                      \
                {                                                             \
                    ppCell[celli]++;                                          \
                    nTotal++;                                                 \
                }                                                             \
            }                                                                 \
            return nTotal;                                                    \
        }

    COUNT_CLOUD(basicReactingMultiphaseCloud)
    else COUNT_CLOUD(basicReactingCloud)
    else COUNT_CLOUD(basicThermoCloud)
    else COUNT_CLOUD(basicKinematicCollidingCloud)
    else COUNT_CLOUD(basicKinematicMPPICCloud)
    else COUNT_CLOUD(basicKinematicCloud)

    #undef COUNT_CLOUD

    return nTotal;
}


void Foam::kinematicCloudHandler::updateMesh(cloud& c)
{
    // Helper macro to call updateMesh on the correct cloud type
    #define UPDATE_CLOUD_MESH(CloudType)                                      \
        if (auto* ptr = dynamic_cast<CloudType*>(&c))                         \
        {                                                                     \
            ptr->updateMesh();                                                \
            Info<< "    Cloud '" << c.name()                                  \
                << "': mesh-dependent data updated" << endl;                  \
            return;                                                           \
        }

    UPDATE_CLOUD_MESH(basicReactingMultiphaseCloud)
    else UPDATE_CLOUD_MESH(basicReactingCloud)
    else UPDATE_CLOUD_MESH(basicThermoCloud)
    else UPDATE_CLOUD_MESH(basicKinematicCollidingCloud)
    else UPDATE_CLOUD_MESH(basicKinematicMPPICCloud)
    else UPDATE_CLOUD_MESH(basicKinematicCloud)

    #undef UPDATE_CLOUD_MESH
}


// ************************************************************************* //
