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
    Foam::cloudHandler

Description
    Abstract base class for cloud handling during mesh redistribution and
    topology changes. Uses OpenFOAM's Run-Time Selection (RTS) mechanism
    to allow users to add support for custom cloud types via external libraries.

    To add support for a custom cloud type:
    1. Derive from cloudHandler
    2. Implement canHandle() to return true for your cloud type
    3. Implement the four virtual methods for your cloud operations
    4. Add static registration: addToRunTimeSelectionTable(cloudHandler, myHandler, cloud)
    5. Load your library via system/controlDict: libs (libmyCloudHandler);

SourceFiles
    cloudHandler.C

\*---------------------------------------------------------------------------*/

#ifndef cloudHandler_H
#define cloudHandler_H

#include "cloud.H"
#include "fvMesh.H"
#include "PstreamBuffers.H"
#include "runTimeSelectionTables.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cloudHandler Declaration
\*---------------------------------------------------------------------------*/

class cloudHandler
{
public:

    //- Runtime type information
    TypeName("cloudHandler");


    // Declare run-time selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        cloudHandler,
        cloud,
        (),
        ()
    );


    // Constructors

        //- Default construct
        cloudHandler() = default;


    // Selectors

        //- Try to find a handler that can process the given cloud
        //  Returns nullptr if no handler can handle this cloud type
        //  If handlerName is provided, uses that handler directly (must exist in RTS)
        //  If handlerName is empty, searches registered handlers using canHandle()
        static autoPtr<cloudHandler> tryNew
        (
            cloud& c,
            const fvMesh& mesh,
            const word& handlerName = word::null
        );

        //- Create handler by explicit name (fails if not found)
        static autoPtr<cloudHandler> New(const word& handlerName);


    //- Destructor
    virtual ~cloudHandler() = default;


    // Member Functions

        //- Check if this handler can process the given cloud type
        //  Derived classes should use dynamic_cast to test the cloud type
        virtual bool canHandle(cloud& c) const = 0;

        //- Priority for handler selection (higher = tried first)
        //  Use higher priority for more derived cloud types to ensure
        //  correct handler selection in inheritance hierarchies
        virtual label priority() const
        {
            return 0;
        }

        //- Store global positions before mesh changes
        virtual void storePositions(cloud& c) = 0;

        //- Distribute parcels to target processors
        //  Called BEFORE mesh distribution with cell->processor mapping
        //  Should extract parcel data from the cloud and send to target procs
        virtual void distribute
        (
            cloud& c,
            const fvMesh& mesh,
            const labelList& distribution,
            PstreamBuffers& pBufs
        ) = 0;

        //- Relocate parcels after mesh redistribution
        //  Called AFTER mesh distribution to recreate parcels in new cells
        virtual void relocate(cloud& c, const fvMesh& mesh) = 0;

        //- Count particles per cell for this cloud
        //  Returns total count for this cloud
        virtual label countPerCell(cloud& c, labelField& ppCell) const = 0;

        //- Update mesh-dependent data (e.g., injection cell indices)
        //  Called after mesh redistribution to refresh cached cell lookups.
        //  This is essential for injection models that cache cell indices,
        //  as those indices become invalid after mesh topology changes.
        virtual void updateMesh(cloud& c) = 0;

        //- Remap particles after mesh refinement/unrefinement
        //  Called after topology changes (via mapPolyMesh).
        //  Default implementation calls the cloud's virtual autoMap method.
        //  Override for custom error handling (e.g., to gracefully handle
        //  particles that cannot be relocated after topology changes).
        //  I hate that this is a virtual with default implmentation
        virtual void autoMap(cloud& c, const mapPolyMesh& map);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
