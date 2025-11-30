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
    Foam::passiveCloudHandler

Description
    Cloud handler for passive particle clouds (Cloud<passiveParticle>).

    Passive particles only store position, so this handler transfers
    positions only (not velocity/diameter/etc. like kinematic clouds).

SourceFiles
    passiveCloudHandler.C

\*---------------------------------------------------------------------------*/

#ifndef passiveCloudHandler_H
#define passiveCloudHandler_H

#include "cloudHandler.H"
#include "DynamicList.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class passiveCloudHandler Declaration
\*---------------------------------------------------------------------------*/

class passiveCloudHandler
:
    public cloudHandler
{
private:

    // Private Static Data

        //- Storage for pending positions during distribution
        //  Key: cloudName, Value: list of positions
        static HashTable<DynamicList<point>> pendingPositions_;

        //- Storage for global positions during autoMap
        //  Key: cloudName, Value: list of positions (stored before topology change)
        static HashTable<DynamicList<point>> autoMapPositions_;


public:

    //- Runtime type information
    TypeName("passiveCloudHandler");


    // Constructors

        //- Default construct
        passiveCloudHandler() = default;


    //- Destructor
    virtual ~passiveCloudHandler() = default;


    // Member Functions

        //- Check if this handler can process the given cloud type
        virtual bool canHandle(cloud& c) const override;

        //- Priority (use 5 - lower than kinematic as it's more generic)
        virtual label priority() const override
        {
            return 5;
        }

        //- Store global positions before mesh changes
        virtual void storePositions(cloud& c) override;

        //- Distribute parcels to target processors
        virtual void distribute
        (
            cloud& c,
            const fvMesh& mesh,
            const labelList& distribution,
            PstreamBuffers& pBufs
        ) override;

        //- Relocate parcels after mesh redistribution
        virtual void relocate(cloud& c, const fvMesh& mesh) override;

        //- Count particles per cell
        virtual label countPerCell(cloud& c, labelField& ppCell) const override;

        //- Update mesh-dependent data (no-op for passive clouds)
        virtual void updateMesh(cloud& c) override;

        //- Remap particles after mesh refinement/unrefinement
        //  Override to gracefully handle particles that cannot be located
        virtual void autoMap(cloud& c, const mapPolyMesh& map) override;


    // Static Member Functions

        //- Access pending positions
        static HashTable<DynamicList<point>>& pendingData()
        {
            return pendingPositions_;
        }

        //- Clear pending data for a specific cloud
        static void clearPendingData(const word& cloudName)
        {
            pendingPositions_.erase(cloudName);
        }

        //- Clear all pending data
        static void clearAllPendingData()
        {
            pendingPositions_.clear();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
