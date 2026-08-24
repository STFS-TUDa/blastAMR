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
    Foam::oversetHandler

Description
    Keeps overset book-keeping consistent across AMR and load balancing.

    OpenFOAM discards the cellCellStencilObject on every topology change (it
    is a MoveableMeshObject, so meshObject::updateMesh destroys it) and
    rebuilds it on the next lookup. What it does not handle is the zoneID
    labelIOList the stencil is built from: it is a plain registered list, so
    it is neither mapped on refinement nor distributed on balancing, and on
    restart it is re-read from the original, undecomposed mesh.

    The handler remaps zoneID itself: by cellMap on a topology change, and
    by the redistribution map on balancing. Mirroring it into a registered
    volScalarField is not an option — the overset machinery interpolates
    every registered field into the acceptor cells except the ones the
    stencil names in nonInterpolatedFields_, so a mirror would be smeared
    across the fringe.

    Selection is automatic for meshes carrying an overset patch. Users with
    custom stencil code derive from this class, add themselves to the
    run-time selection table and select their implementation with

    \verbatim
        oversetHandler  myHandler;   // or "none" to disable entirely
    \endverbatim

    in the AMR dictionary.

SourceFiles
    oversetHandler.C

\*---------------------------------------------------------------------------*/

#ifndef oversetHandler_H
#define oversetHandler_H

#include "fvMesh.H"
#include "labelIOList.H"
#include "wordRes.H"
#include "regIOobject.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                       Class oversetHandler Declaration
\*---------------------------------------------------------------------------*/

class oversetHandler
:
    public regIOobject
{
protected:

    // Protected Data

        //- Non-const access to the mesh
        fvMesh& oversetMesh_;


    // Protected Member Functions

        //- The zone identifiers, sized nCells. Created on first call
        labelIOList& zoneID();

        //- Drop the cached stencil so the next lookup rebuilds it
        void expireStencil();


    // Protected Data

        //- zoneID as it was before redistribution
        labelList oldZoneID_;

        //- Cell zones excluded from adaptation, see protectCells()
        wordRes protectZones_;


public:

    //- Runtime type information
    TypeName("oversetHandler");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            oversetHandler,
            mesh,
            (fvMesh& mesh),
            (mesh)
        );


    // Constructors

        //- Construct from mesh
        explicit oversetHandler(fvMesh& mesh);

        //- No copy construct
        oversetHandler(const oversetHandler&) = delete;


    //- Destructor
    virtual ~oversetHandler() = default;


    // Selectors

        //- Does the mesh carry an overset patch?
        static bool isOverset(const fvMesh& mesh);

        //- The handler registered on the mesh, or nullptr
        static oversetHandler* lookup(const fvMesh& mesh);

        //- Register a handler on an overset mesh if there is none yet.
        //  Returns the registered handler, or nullptr if the mesh is not
        //  overset or handling was switched off with "oversetHandler none"
        static oversetHandler* New(fvMesh& mesh, const dictionary& dict);


    // Member Functions

        //- Read settings
        virtual void read(const dictionary& dict);

        //- Exclude cells from adaptation by zeroing their error estimate:
        //  hole cells, and any cell in a zone named by protectZones.
        //  Returns the number of cells protected on this processor
        virtual label protectCells(volScalarField& error) const;

        //- Remap zoneID onto the new cells. Must run *before*
        //  fvMesh::updateMesh: correcting the overset boundary conditions
        //  in there rebuilds the stencil, which reads zoneID
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Take a copy of zoneID on the pre-distribution mesh
        virtual void preDistribute();

        //- Redistribute zoneID onto the new decomposition
        virtual void distribute(const mapDistributePolyMesh& map);

        //- Write zoneID for the current mesh so that a restart reads a
        //  zoneID matching the adapted mesh rather than the original one
        virtual void writeZoneID() const;

        //- Nothing of our own to write
        virtual bool writeData(Ostream&) const
        {
            return true;
        }


    // Member Operators

        //- No copy assignment
        void operator=(const oversetHandler&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
