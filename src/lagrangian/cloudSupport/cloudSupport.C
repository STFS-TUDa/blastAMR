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
#include "globalIndex.H"
#include "IOdictionary.H"

// Handler system
#include "cloudHandler.H"
#include "kinematicCloudHandler.H"
#include "passiveCloudHandler.H"

// Debug switch for detailed particle tracking output
int cloudSupportDebug = Foam::debug::debugSwitch("cloudSupport", 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cloudSupport
{
    // Track last timestep when autoMapClouds was called to prevent double-mapping
    static label lastAutoMapTimeIndex_ = -1;

    // Cache handlers for each cloud to avoid repeated lookups
    static HashTable<autoPtr<cloudHandler>> handlerCache_;
}
}

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
namespace cloudSupport
{

// Get the cloudHandler name from dynamicMeshDict (if specified)
static word getConfiguredHandlerName(const fvMesh& mesh)
{
    // Try to read dynamicMeshDict
    IOobject dictIO
    (
        "dynamicMeshDict",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (dictIO.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary dict(dictIO);
        return dict.getOrDefault<word>("cloudHandler", word::null);
    }

    return word::null;
}


// Get or create handler for a cloud (returns reference to cached handler)
// Returns nullptr equivalent via empty autoPtr if no handler found
static autoPtr<cloudHandler>& getHandler(cloud& c, const fvMesh& mesh)
{
    static autoPtr<cloudHandler> nullHandler;  // Empty handler for "not found"
    const word& cloudName = c.name();

    // Check cache first
    if (handlerCache_.found(cloudName) && handlerCache_[cloudName])
    {
        // Verify handler still valid for this cloud
        if (handlerCache_[cloudName]->canHandle(c))
        {
            return handlerCache_[cloudName];
        }
        // Handler no longer valid, remove from cache
        handlerCache_.erase(cloudName);
    }

    // Check if user specified an explicit handler in dynamicMeshDict
    const word handlerName = getConfiguredHandlerName(mesh);

    // Try to create new handler
    autoPtr<cloudHandler> handler = cloudHandler::tryNew(c, mesh, handlerName);
    if (handler)
    {
        handlerCache_.set(cloudName, std::move(handler));
        return handlerCache_[cloudName];
    }

    nullHandler.reset();
    return nullHandler;
}

} // End namespace cloudSupport
} // End namespace Foam


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

        autoPtr<cloudHandler>& handler = getHandler(c, mesh);
        if (handler)
        {
            handler->storePositions(c);
        }
        else
        {
            // Unknown cloud type - warn but don't fail
            WarningInFunction
                << "No handler found for cloud '" << c.name()
                << "' - positions not stored. "
                << "Load balancing may fail for this cloud." << nl
                << "Consider loading a library with a handler for this cloud type "
                << "via system/controlDict: libs (libmyCloudHandler);" << endl;
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

        autoPtr<cloudHandler>& handler = getHandler(c, mesh);
        if (handler)
        {
            handler->autoMap(c, map);
        }
        else
        {
            c.autoMap(map);
        }

        Info<< "    Cloud '" << c.name()
            << "': " << oldSize << " -> " << c.nParcels()
            << " particles" << endl;
    }
}


void Foam::cloudSupport::distributeClouds
(
    const fvMesh& mesh,
    const labelList& distribution
)
{
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

        if (cloudSupportDebug)
        {
            const label oldSize = c.nParcels();
            Pout<< "    Cloud '" << cloudName << "' on proc "
                << UPstream::myProcNo() << ": " << oldSize << " parcels" << endl;
        }

        autoPtr<cloudHandler>& handler = getHandler(c, mesh);
        if (!handler)
        {
            WarningInFunction
                << "No handler found for cloud '" << cloudName
                << "' - cannot distribute particles" << endl;
            continue;
        }

        // Transfer buffers
        PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);

        // Use handler to distribute
        // NOTE: cloud may appear empty if storePositions() cleared it,
        // but handler has stored positions to distribute
        handler->distribute(c, mesh, distribution, pBufs);

        Info<< "    Cloud '" << cloudName << "': distributed" << endl;
    }
}


void Foam::cloudSupport::relocateClouds(const fvMesh& mesh)
{
    // After mesh redistribution, this method uses handlers to create parcels
    // from the data that was stored during distributeClouds().

    // Check if any handler has pending data
    bool hasPendingData =
        !kinematicCloudHandler::pendingData().empty()
     || !passiveCloudHandler::pendingData().empty();

    if (!hasPendingData)
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

    // Process each cloud that might have pending data
    UPtrList<const cloud> allClouds = mesh.csorted<cloud>();

    for (const cloud& constCloud : allClouds)
    {
        cloud& c = const_cast<cloud&>(constCloud);

        autoPtr<cloudHandler>& handler = getHandler(c, mesh);
        if (handler)
        {
            handler->relocate(c, mesh);
        }
    }

    // Update mesh-dependent data (injection cell indices, cell occupancy, etc.)
    // This is essential for injection models like ManualInjection that cache
    // cell indices - those indices become invalid after mesh redistribution.
    Info<< "cloudSupport: Updating mesh-dependent data after redistribution" << endl;

    for (const cloud& constCloud : allClouds)
    {
        cloud& c = const_cast<cloud&>(constCloud);

        autoPtr<cloudHandler>& handler = getHandler(c, mesh);
        if (handler)
        {
            handler->updateMesh(c);
        }
    }

    // Clear any remaining handler caches (mesh topology changed)
    handlerCache_.clear();
}


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
    auto tppCell = tmp<labelField>::New(mesh.nCells(), Zero);
    labelField& ppCell = tppCell.ref();

    UPtrList<const cloud> allClouds = mesh.csorted<cloud>();

    for (const cloud& constCloud : allClouds)
    {
        cloud& c = const_cast<cloud&>(constCloud);

        autoPtr<cloudHandler>& handler = getHandler(c, mesh);
        if (handler)
        {
            handler->countPerCell(c, ppCell);
        }
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
