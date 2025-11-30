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
    Foam::kinematicCloudHandler

Description
    Cloud handler for kinematic-derived cloud types:
    - basicKinematicCloud
    - basicKinematicCollidingCloud
    - basicKinematicMPPICCloud
    - basicThermoCloud
    - basicReactingCloud
    - basicReactingMultiphaseCloud

    All these cloud types share the same kinematic parcel interface and can
    be handled uniformly by checking for basicKinematicCloud (the base).

    Since the handler uses dynamic_cast to basicKinematicCloud, it will match
    all derived types. The priority mechanism ensures more specific handlers
    (if any) take precedence.

SourceFiles
    kinematicCloudHandler.C

\*---------------------------------------------------------------------------*/

#ifndef kinematicCloudHandler_H
#define kinematicCloudHandler_H

#include "cloudHandler.H"
#include "DynamicList.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class kinematicCloudHandler Declaration
\*---------------------------------------------------------------------------*/

class kinematicCloudHandler
:
    public cloudHandler
{
public:

    //- Data structure for transferring parcel properties between processors
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


private:

    // Private Static Data

        //- Storage for pending parcel data during distribution
        //  Key: cloudName, Value: list of parcel data
        static HashTable<DynamicList<ParcelData>> pendingParcelData_;

        //- Storage for parcel data during autoMap (mesh topology changes)
        //  Key: cloudName, Value: list of parcel data
        static HashTable<DynamicList<ParcelData>> autoMapParcelData_;


public:

    //- Runtime type information
    TypeName("kinematicCloudHandler");


    // Constructors

        //- Default construct
        kinematicCloudHandler() = default;


    //- Destructor
    virtual ~kinematicCloudHandler() = default;


    // Member Functions

        //- Check if this handler can process the given cloud type
        virtual bool canHandle(cloud& c) const override;

        //- Priority (use 10 for kinematic base - derived handlers can use higher)
        virtual label priority() const override
        {
            return 10;
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

        //- Update mesh-dependent data (injection cell indices, cell occupancy)
        virtual void updateMesh(cloud& c) override;

        //- Remap parcels after mesh topology change (refinement/unrefinement)
        virtual void autoMap(cloud& c, const mapPolyMesh& map) override;


    // Static Member Functions

        //- Access pending parcel data (for relocate after mesh redistribution)
        static HashTable<DynamicList<ParcelData>>& pendingData()
        {
            return pendingParcelData_;
        }

        //- Clear pending data for a specific cloud
        static void clearPendingData(const word& cloudName)
        {
            pendingParcelData_.erase(cloudName);
        }

        //- Clear all pending data
        static void clearAllPendingData()
        {
            pendingParcelData_.clear();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
