/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "volFields.H"
#include "surfaceFields.H"
#include "emptyFvPatchField.H"
#include "directFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMeshAdder::MapVolField
(
    const mapAddedPolyMesh& meshMap,

    GeometricField<Type, fvPatchField, volMesh>& fld,
    const GeometricField<Type, fvPatchField, volMesh>& fldToAdd,
    const bool fullyMapped
)
{
    const fvMesh& mesh = fld.mesh();

    // Internal field
    // ~~~~~~~~~~~~~~

    {
        // Store old internal field
        const Field<Type> oldInternalField(fld.primitiveField());

        // Modify internal field
        Field<Type>& intFld = fld.primitiveFieldRef();

        intFld.setSize(mesh.nCells());

        intFld.rmap(oldInternalField, meshMap.oldCellMap());
        intFld.rmap(fldToAdd.primitiveField(), meshMap.addedCellMap());
    }


    // Patch fields from old mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    auto& bfld = fld.boundaryFieldRef();

    {
        const labelList& oldPatchMap = meshMap.oldPatchMap();
        const labelList& oldPatchStarts = meshMap.oldPatchStarts();
        const labelList& oldPatchSizes = meshMap.oldPatchSizes();

        // Reorder old patches in order of new ones. Put removed patches at end.

        label unusedPatchi = 0;

        forAll(oldPatchMap, patchi)
        {
            label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                unusedPatchi++;
            }
        }

        label nUsedPatches = unusedPatchi;

        // Reorder list for patchFields
        labelList oldToNew(oldPatchMap.size());

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                oldToNew[patchi] = newPatchi;
            }
            else
            {
                oldToNew[patchi] = unusedPatchi++;
            }
        }


        // Sort deleted ones last so is now in newPatch ordering
        bfld.reorder(oldToNew);
        // Extend to covers all patches
        bfld.setSize(mesh.boundaryMesh().size());
        // Delete unused patches
        for
        (
            label newPatchi = nUsedPatches;
            newPatchi < bfld.size();
            newPatchi++
        )
        {
            bfld.set(newPatchi, nullptr);
        }


        // Map old values
        // ~~~~~~~~~~~~~~

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                labelList newToOld
                (
                    calcPatchMap
                    (
                        oldPatchStarts[patchi],
                        oldPatchSizes[patchi],
                        meshMap.oldFaceMap(),
                        mesh.boundaryMesh()[newPatchi],
                        -1              // unmapped value
                    )
                );

                directFvPatchFieldMapper patchMapper(newToOld);

                // Override mapping (for use in e.g. fvMeshDistribute where
                // it sorts mapping out itself)
                if (fullyMapped)
                {
                    patchMapper.hasUnmapped() = false;
                }

                // Create new patchField with same type as existing one.
                // Note:
                // - boundaryField already in new order so access with newPatchi
                // - fld.boundaryField()[newPatchi] both used for type and old
                //   value
                // - hope that field mapping allows aliasing since old and new
                //   are same memory!
                bfld.set
                (
                    newPatchi,
                    fvPatchField<Type>::New
                    (
                        bfld[newPatchi],                // old field
                        mesh.boundary()[newPatchi],     // new fvPatch
                        fld(), // new internal field
                        patchMapper                     // mapper (new to old)
                    )
                );
            }
        }
    }



    // Patch fields from added mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& addedPatchMap = meshMap.addedPatchMap();

        // Add addedMesh patches
        forAll(addedPatchMap, patchi)
        {
            const label newPatchi = addedPatchMap[patchi];

            if (newPatchi != -1)
            {
                const polyPatch& newPatch = mesh.boundaryMesh()[newPatchi];
                const polyPatch& oldPatch =
                    fldToAdd.mesh().boundaryMesh()[patchi];

                if (!bfld(newPatchi))
                {
                    // First occurrence of newPatchi. Map from existing
                    // patchField

                    // From new patch faces to patch faces on added mesh.
                    labelList newToAdded
                    (
                        calcPatchMap
                        (
                            oldPatch.start(),
                            oldPatch.size(),
                            meshMap.addedFaceMap(),
                            newPatch,
                            -1          // unmapped values
                        )
                    );

                    directFvPatchFieldMapper patchMapper(newToAdded);

                    // Override mapping (for use in e.g. fvMeshDistribute where
                    // it sorts mapping out itself)
                    if (fullyMapped)
                    {
                        patchMapper.hasUnmapped() = false;
                    }

                    bfld.set
                    (
                        newPatchi,
                        fvPatchField<Type>::New
                        (
                            fldToAdd.boundaryField()[patchi], // added field
                            mesh.boundary()[newPatchi],       // new fvPatch
                            fld(),   // new int. field
                            patchMapper                       // mapper
                        )
                    );
                }
                else
                {
                    // PatchField will have correct size already. Just slot in
                    // my elements.

                    labelList addedToNew(oldPatch.size(), -1);
                    forAll(addedToNew, i)
                    {
                        label addedFacei = oldPatch.start()+i;
                        label newFacei = meshMap.addedFaceMap()[addedFacei];
                        label patchFacei = newFacei-newPatch.start();
                        if (patchFacei >= 0 && patchFacei < newPatch.size())
                        {
                            addedToNew[i] = patchFacei;
                        }
                    }

                    bfld[newPatchi].rmap
                    (
                        fldToAdd.boundaryField()[patchi],
                        addedToNew
                    );
                }
            }
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapVolFields
(
    const mapAddedPolyMesh& meshMap,
    const fvMesh& mesh,
    const fvMesh& meshToAdd,
    const bool fullyMapped
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fldType;

    HashTable<const fldType*> fields
    (
        mesh.objectRegistry::lookupClass<fldType>()
    );

    HashTable<const fldType*> fieldsToAdd
    (
        meshToAdd.objectRegistry::lookupClass<fldType>()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    forAllIters(fields, fieldIter)
    {
        fldType& fld = const_cast<fldType&>(*fieldIter());

        DebugPout
            << "MapVolFields : Storing old time for " << fld.name() << endl;

        fld.storeOldTimes();
    }


    forAllIters(fields, fieldIter)
    {
        fldType& fld = const_cast<fldType&>(*fieldIter());

        if (fieldsToAdd.found(fld.name()))
        {
            const fldType& fldToAdd = *fieldsToAdd[fld.name()];

            DebugPout
                << "MapVolFields : mapping " << fld.name()
                << " and " << fldToAdd.name() << endl;

            MapVolField<Type>(meshMap, fld, fldToAdd, fullyMapped);
        }
        else
        {
            WarningInFunction
                << "Not mapping field " << fld.name()
                << " since not present on mesh to add" << endl;
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapSurfaceField
(
    const mapAddedPolyMesh& meshMap,

    GeometricField<Type, fvsPatchField, surfaceMesh>& fld,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fldToAdd,
    const bool fullyMapped
)
{
    const fvMesh& mesh = fld.mesh();
    const labelList& oldPatchStarts = meshMap.oldPatchStarts();

    auto& bfld = fld.boundaryFieldRef();

    // Internal field
    // ~~~~~~~~~~~~~~

    // Store old internal field
    {
        const Field<Type> oldField(fld);

        // Modify internal field
        Field<Type>& intFld = fld.primitiveFieldRef();

        intFld.setSize(mesh.nInternalFaces());

        intFld.rmap(oldField, meshMap.oldFaceMap());
        intFld.rmap(fldToAdd, meshMap.addedFaceMap());


        // Faces that were boundary faces but are not anymore.
        // Use owner value (so lowest numbered cell, i.e. from 'old' not 'added'
        // mesh)
        forAll(bfld, patchi)
        {
            const fvsPatchField<Type>& pf = bfld[patchi];

            label start = oldPatchStarts[patchi];

            forAll(pf, i)
            {
                label newFacei = meshMap.oldFaceMap()[start + i];

                if (newFacei >= 0 && newFacei < mesh.nInternalFaces())
                {
                    intFld[newFacei] = pf[i];
                }
            }
        }
    }


    // Patch fields from old mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& oldPatchMap = meshMap.oldPatchMap();
        const labelList& oldPatchSizes = meshMap.oldPatchSizes();

        // Reorder old patches in order of new ones. Put removed patches at end.

        label unusedPatchi = 0;

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                unusedPatchi++;
            }
        }

        label nUsedPatches = unusedPatchi;

        // Reorder list for patchFields
        labelList oldToNew(oldPatchMap.size());

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                oldToNew[patchi] = newPatchi;
            }
            else
            {
                oldToNew[patchi] = unusedPatchi++;
            }
        }


        // Sort deleted ones last so is now in newPatch ordering
        bfld.reorder(oldToNew);
        // Extend to covers all patches
        bfld.setSize(mesh.boundaryMesh().size());
        // Delete unused patches
        for
        (
            label newPatchi = nUsedPatches;
            newPatchi < bfld.size();
            newPatchi++
        )
        {
            bfld.set(newPatchi, nullptr);
        }


        // Map old values
        // ~~~~~~~~~~~~~~

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                labelList newToOld
                (
                    calcPatchMap
                    (
                        oldPatchStarts[patchi],
                        oldPatchSizes[patchi],
                        meshMap.oldFaceMap(),
                        mesh.boundaryMesh()[newPatchi],
                        -1      // unmapped value
                    )
                );

                directFvPatchFieldMapper patchMapper(newToOld);

                // Override mapping (for use in e.g. fvMeshDistribute where
                // it sorts mapping out itself)
                if (fullyMapped)
                {
                    patchMapper.hasUnmapped() = false;
                }

                // Create new patchField with same type as existing one.
                // Note:
                // - boundaryField already in new order so access with newPatchi
                // - bfld[newPatchi] both used for type and old
                //   value
                // - hope that field mapping allows aliasing since old and new
                //   are same memory!
                bfld.set
                (
                    newPatchi,
                    fvsPatchField<Type>::New
                    (
                        bfld[newPatchi],                // old field
                        mesh.boundary()[newPatchi],     // new fvPatch
                        fld(), // new internal field
                        patchMapper                     // mapper (new to old)
                    )
                );
            }
        }
    }



    // Patch fields from added mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& addedPatchMap = meshMap.addedPatchMap();

        // Add addedMesh patches
        forAll(addedPatchMap, patchi)
        {
            const label newPatchi = addedPatchMap[patchi];

            if (newPatchi != -1)
            {
                const polyPatch& newPatch = mesh.boundaryMesh()[newPatchi];
                const polyPatch& oldPatch =
                    fldToAdd.mesh().boundaryMesh()[patchi];

                if (!bfld(newPatchi))
                {
                    // First occurrence of newPatchi. Map from existing
                    // patchField

                    // From new patch faces to patch faces on added mesh.
                    labelList newToAdded
                    (
                        calcPatchMap
                        (
                            oldPatch.start(),
                            oldPatch.size(),
                            meshMap.addedFaceMap(),
                            newPatch,
                            -1                  // unmapped values
                        )
                    );

                    directFvPatchFieldMapper patchMapper(newToAdded);

                    // Override mapping (for use in e.g. fvMeshDistribute where
                    // it sorts mapping out itself)
                    if (fullyMapped)
                    {
                        patchMapper.hasUnmapped() = false;
                    }

                    bfld.set
                    (
                        newPatchi,
                        fvsPatchField<Type>::New
                        (
                            fldToAdd.boundaryField()[patchi],// added field
                            mesh.boundary()[newPatchi],      // new fvPatch
                            fld(),  // new int. field
                            patchMapper                      // mapper
                        )
                    );
                }
                else
                {
                    // PatchField will have correct size already. Just slot in
                    // my elements.

                    labelList addedToNew(oldPatch.size(), -1);
                    forAll(addedToNew, i)
                    {
                        label addedFacei = oldPatch.start()+i;
                        label newFacei = meshMap.addedFaceMap()[addedFacei];
                        label patchFacei = newFacei-newPatch.start();
                        if (patchFacei >= 0 && patchFacei < newPatch.size())
                        {
                            addedToNew[i] = patchFacei;
                        }
                    }

                    bfld[newPatchi].rmap
                    (
                        fldToAdd.boundaryField()[patchi],
                        addedToNew
                    );
                }
            }
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapSurfaceFields
(
    const mapAddedPolyMesh& meshMap,
    const fvMesh& mesh,
    const fvMesh& meshToAdd,
    const bool fullyMapped
)
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fldType;

    HashTable<const fldType*> fields
    (
        mesh.objectRegistry::lookupClass<fldType>()
    );

    HashTable<const fldType*> fieldsToAdd
    (
        meshToAdd.objectRegistry::lookupClass<fldType>()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    forAllIters(fields, fieldIter)
    {
        fldType& fld = const_cast<fldType&>(*fieldIter());

        DebugPout
            << "MapSurfaceFields : Storing old time for " << fld.name() << endl;

        fld.storeOldTimes();
    }


    forAllIters(fields, fieldIter)
    {
        fldType& fld = const_cast<fldType&>(*fieldIter());

        if (fieldsToAdd.found(fld.name()))
        {
            const fldType& fldToAdd = *fieldsToAdd[fld.name()];

            DebugPout
                << "MapSurfaceFields : mapping " << fld.name()
                << " and " << fldToAdd.name() << endl;

            MapSurfaceField<Type>(meshMap, fld, fldToAdd, fullyMapped);
        }
        else
        {
            WarningInFunction
                << "Not mapping field " << fld.name()
                << " since not present on mesh to add" << endl;
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapDimField
(
    const mapAddedPolyMesh& meshMap,

    DimensionedField<Type, volMesh>& fld,
    const DimensionedField<Type, volMesh>& fldToAdd
)
{
    const fvMesh& mesh = fld.mesh();

    // Store old field
    const Field<Type> oldField(fld);

    fld.setSize(mesh.nCells());

    fld.rmap(oldField, meshMap.oldCellMap());
    fld.rmap(fldToAdd, meshMap.addedCellMap());
}


template<class Type>
void Foam::fvMeshAdder::MapDimFields
(
    const mapAddedPolyMesh& meshMap,
    const fvMesh& mesh,
    const fvMesh& meshToAdd
)
{
    typedef DimensionedField<Type, volMesh> fldType;

    // Note: use strict flag on lookupClass to avoid picking up
    //       volFields
    HashTable<const fldType*> fields
    (
        mesh.objectRegistry::lookupClass<fldType>(true)
    );

    HashTable<const fldType*> fieldsToAdd
    (
        meshToAdd.objectRegistry::lookupClass<fldType>(true)
    );

    forAllIters(fields, fieldIter)
    {
        fldType& fld = const_cast<fldType&>(*fieldIter());

        if (fieldsToAdd.found(fld.name()))
        {
            const fldType& fldToAdd = *fieldsToAdd[fld.name()];

            DebugPout
                << "MapDimFields : mapping " << fld.name()
                << " and " << fldToAdd.name() << endl;

            MapDimField<Type>(meshMap, fld, fldToAdd);
        }
        else
        {
            WarningInFunction
                << "Not mapping field " << fld.name()
                << " since not present on mesh to add" << endl;
        }
    }
}


//- Multi-mesh mapping

template<class Type>
void Foam::fvMeshAdder::MapDimField
(
    UPtrList<DimensionedField<Type, volMesh>>& flds,
    const labelListList& cellProcAddressing,
    const bool fullyMapped
)
{
    // Add fields to fields[0] after adding the meshes to meshes[0].
    // Mesh[0] is the sum of all meshes. Fields are not yet mapped.

    if
    (
        flds.size() == 0
    || !flds.set(0)
    || cellProcAddressing.size() != flds.size()
    )
    {
        FatalErrorInFunction << "Not valid field at element 0"
            << " in field list of size " << flds.size() << exit(FatalError);
    }


    // Internal field
    // ~~~~~~~~~~~~~~

    {
        // Store old internal field
        const Field<Type> oldInternalField(flds[0]);

        // Modify internal field
        Field<Type>& intFld = flds[0];

        // Set to new mesh size
        intFld.setSize(flds[0].mesh().nCells());
        // Add fld0
        intFld.rmap(oldInternalField, cellProcAddressing[0]);

        for (label meshi = 1; meshi < flds.size(); meshi++)
        {
            if (flds.set(meshi))
            {
                const Field<Type>& addFld = flds[meshi];
                intFld.rmap(addFld, cellProcAddressing[meshi]);
            }
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapVolField
(
    UPtrList<GeometricField<Type, fvPatchField, volMesh>>& flds,
    const labelList& oldPatchStarts0,
    const labelList& oldPatchSizes0,
    const labelListList& patchProcAddressing,
    const labelListList& cellProcAddressing,
    const labelListList& faceProcAddressing,
    const bool fullyMapped
)
{
    // Add fields to fields[0] after adding the meshes to meshes[0].
    // Mesh[0] is the sum of all meshes. Fields are not yet mapped.

    if (flds.size() == 0 || !flds.set(0))
    {
        FatalErrorInFunction << "Not valid field at element 0"
            << " in field list of size " << flds.size() << exit(FatalError);
    }


    // Internal field
    // ~~~~~~~~~~~~~~

    {
        // Store old internal field
        const Field<Type> oldInternalField(flds[0].primitiveField());

        // Modify internal field
        Field<Type>& intFld = flds[0].primitiveFieldRef();

        // Set to new mesh size
        intFld.setSize(flds[0].mesh().nCells());
        // Add fld0
        intFld.rmap(oldInternalField, cellProcAddressing[0]);

        for (label meshi = 1; meshi < flds.size(); meshi++)
        {
            if (flds.set(meshi))
            {
                const Field<Type>& addFld = flds[meshi].primitiveFieldRef();
                intFld.rmap(addFld, cellProcAddressing[meshi]);
            }
        }
    }


    // Patch fields from old mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    auto& bfld0 = flds[0].boundaryFieldRef();
    //Pout<< "   Adding patchfields for field " << flds[0].name() << endl;
    forAll(bfld0, patchi)
    {
        //Pout<< "    patch:" << patchi
        //    << " patch:" << flds[0].mesh().boundaryMesh()[patchi].name()
        //    << " size:" << flds[0].mesh().boundaryMesh()[patchi].size()
        //    << endl;
        //Pout<< "    patchField:" << bfld0[patchi].size()
        //    << endl;

        labelList newToOld
        (
            calcPatchMap
            (
                oldPatchStarts0[patchi],
                oldPatchSizes0[patchi],
                faceProcAddressing[0],
                bfld0[patchi].patch().patch(),
                -1              // unmapped value
            )
        );

        directFvPatchFieldMapper patchMapper(newToOld);

        // Override mapping (for use in e.g. fvMeshDistribute where
        // it sorts mapping out itself)
        if (fullyMapped)
        {
            patchMapper.hasUnmapped() = false;
        }

        bfld0[patchi].autoMap(patchMapper);
    }

    for (label meshi = 1; meshi < flds.size(); meshi++)
    {
        if (flds.set(meshi))
        {
            const auto& bfld = flds[meshi].boundaryFieldRef();

            const labelList& patchMap = patchProcAddressing[meshi];

            forAll(patchMap, oldPatchi)
            {
                const auto& fvp = bfld[oldPatchi].patch();
                const label newPatchi = patchMap[oldPatchi];

                //Pout<< "    oldPatch:" << oldPatchi
                //    << " newPatch:" << newPatchi << endl;

                if (newPatchi >= 0 && newPatchi < bfld0.size())
                {
                    const auto& fvp0 = bfld0[newPatchi].patch();
                    labelList addedToNew(bfld[oldPatchi].size(), -1);
                    forAll(addedToNew, i)
                    {
                        const label newFacei =
                            faceProcAddressing[meshi][fvp.start()+i];
                        const label patchFacei = newFacei-fvp0.start();
                        if
                        (
                            patchFacei >= 0
                         && patchFacei < fvp0.size()
                        )
                        {
                            addedToNew[i] = patchFacei;
                        }
                    }

                    bfld0[newPatchi].rmap(bfld[oldPatchi], addedToNew);
                }
                else
                {
                    WarningInFunction << "Ignoring old patch "
                        << bfld[oldPatchi].patch().name() << " on field "
                        << flds[meshi].name() << endl;  //exit(FatalError);
                }
            }
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapSurfaceField
(
    UPtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>& flds,
    const labelList& oldFaceOwner0,
    const labelList& oldPatchStarts0,
    const labelList& oldPatchSizes0,
    const labelListList& patchProcAddressing,
    const labelListList& cellProcAddressing,
    const labelListList& faceProcAddressing,
    const bool fullyMapped
)
{
    // Add fields to fields[0] after adding the meshes to meshes[0].
    // Mesh[0] is the sum of all meshes. Fields are not yet mapped.

    if (flds.size() == 0 || !flds.set(0))
    {
        FatalErrorInFunction << "Not valid field at element 0"
            << " in field list of size " << flds.size() << exit(FatalError);
    }

    const fvMesh& mesh0 = flds[0].mesh();


    // Internal field
    // ~~~~~~~~~~~~~~

    {
        // Store old internal field
        const Field<Type> oldInternalField(flds[0].primitiveField());

        // Modify internal field
        Field<Type>& intFld = flds[0].primitiveFieldRef();

        // Set to new mesh size
        intFld.setSize(mesh0.nInternalFaces());

        // Map
        forAll(flds, meshi)
        {
            if (flds.set(meshi))
            {
                const labelList& faceMap = faceProcAddressing[meshi];
                const auto& fld = flds[meshi];

                // Map internal field
                if (meshi == 0)
                {
                    intFld.rmap(oldInternalField, faceMap);
                }
                else
                {
                    intFld.rmap(fld.primitiveField(), faceMap);
                }

                // Map faces that were boundary faces but are not anymore.
                // There will be two meshes that provide the same face. Use
                // owner one.
                const auto& bfld = flds[meshi].boundaryField();

                forAll(bfld, oldPatchi)
                {
                    const fvsPatchField<Type>& pf = bfld[oldPatchi];
                    //Pout<< "patch:" << pf.patch().name() << endl;
                    forAll(pf, patchFacei)
                    {
                        // Get new face, mapped face owner
                        label newFacei;
                        label newOwn;
                        if (meshi == 0)
                        {
                            // Do not access mesh information since in-place
                            // modified
                            const label oldFacei =
                                oldPatchStarts0[oldPatchi]+patchFacei;
                            newFacei = faceProcAddressing[meshi][oldFacei];
                            const label oldOwn = oldFaceOwner0[oldFacei];
                            newOwn = cellProcAddressing[meshi][oldOwn];

                            //Pout<< "MESH0: pfi:" << patchFacei
                            //    << " old face:" << oldFacei
                            //    << " new face:" << newFacei
                            //    << " at:" << mesh0.faceCentres()[newFacei]
                            //    << " oldOwn:" << oldOwn
                            //    << " newOwn:" << newOwn << endl;
                        }
                        else
                        {
                            const label oldFacei =
                                pf.patch().start()+patchFacei;
                            newFacei = faceProcAddressing[meshi][oldFacei];
                            const label oldOwn =
                                fld.mesh().faceOwner()[oldFacei];
                            newOwn = cellProcAddressing[meshi][oldOwn];

                            //Pout<< "MESH:" << meshi << " pfi:" << patchFacei
                            //    << " old face:" << oldFacei
                            //    << " new face:" << newFacei
                            //    << " at:" << mesh0.faceCentres()[newFacei]
                            //    << " oldOwn:" << oldOwn
                            //    << " newOwn:" << newOwn << endl;
                        }

                        if
                        (
                            newFacei >= 0
                         && newFacei < mesh0.nInternalFaces()
                         && (newOwn == mesh0.faceOwner()[newFacei])
                        )
                        {
                            intFld[newFacei] = pf[patchFacei];
                        }
                        //else
                        //{
                        //    Pout<< "    ignoring pfi:" << patchFacei
                        //        << " value:" << pf[patchFacei]
                        //        << " since newFacei:" << newFacei
                        //        << " since newOwn:" << newOwn
                        //        << " own:" << mesh0.faceOwner()[newFacei]
                        //        << endl;
                        //}
                    }
                }
            }
        }
    }


    // Patch fields from old mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    auto& bfld0 = flds[0].boundaryFieldRef();
    forAll(bfld0, patchi)
    {
        labelList newToOld
        (
            calcPatchMap
            (
                oldPatchStarts0[patchi],
                oldPatchSizes0[patchi],
                faceProcAddressing[0],
                bfld0[patchi].patch().patch(),
                -1              // unmapped value
            )
        );

        directFvPatchFieldMapper patchMapper(newToOld);

        // Override mapping (for use in e.g. fvMeshDistribute where
        // it sorts mapping out itself)
        if (fullyMapped)
        {
            patchMapper.hasUnmapped() = false;
        }

        bfld0[patchi].autoMap(patchMapper);
    }

    for (label meshi = 1; meshi < flds.size(); meshi++)
    {
        if (flds.set(meshi))
        {
            const auto& bfld = flds[meshi].boundaryFieldRef();

            const labelList& patchMap = patchProcAddressing[meshi];

            forAll(patchMap, oldPatchi)
            {
                const auto& fvp = bfld[oldPatchi].patch();
                const label newPatchi = patchMap[oldPatchi];
                if (newPatchi >= 0 && newPatchi < bfld0.size())
                {
                    const auto& fvp0 = bfld0[newPatchi].patch();
                    labelList addedToNew(bfld[oldPatchi].size(), -1);
                    forAll(addedToNew, i)
                    {
                        const label newFacei =
                            faceProcAddressing[meshi][fvp.start()+i];
                        const label patchFacei = newFacei-fvp0.start();
                        if
                        (
                            patchFacei >= 0
                         && patchFacei < fvp0.size()
                        )
                        {
                            addedToNew[i] = patchFacei;
                        }
                    }

                    bfld0[newPatchi].rmap(bfld[oldPatchi], addedToNew);
                }
                else
                {
                    WarningInFunction << "Ignoring old patch "
                        << bfld[oldPatchi].patch().name() << " on field "
                        << flds[meshi].name() << endl;  //exit(FatalError);
                }
            }
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapVolFields
(
    const UPtrList<fvMesh>& meshes,
    const labelList& oldPatchStarts0,
    const labelList& oldPatchSizes0,
    const labelListList& patchProcAddressing,
    const labelListList& cellProcAddressing,
    const labelListList& faceProcAddressing,
    const labelListList& pointProcAddressing,
    const bool fullyMapped
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fldType;

    if (meshes.size() == 0 || !meshes.set(0))
    {
        FatalErrorInFunction << "Not valid field at element 0"
            << " in mesh list of size " << meshes.size() << exit(FatalError);
    }
    const fvMesh& mesh0 = meshes[0];

    HashTable<const fldType*> fields
    (
        mesh0.objectRegistry::lookupClass<fldType>()
    );


    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    for (const auto& fld : fields)
    {
        DebugPout
            << "MapVolFields : Storing old time for " << fld->name()
            << endl;

        const_cast<fldType&>(*fld).storeOldTimes();
    }


    for (const auto& fld : fields)
    {
        const word& name0 = fld->name();

        DebugPout
            << "MapVolFields : mapping " << name0 << endl;

        UPtrList<fldType> meshToField(meshes.size());
        forAll(meshes, meshi)
        {
            if (meshes.set(meshi))
            {
                auto& meshFld = meshes[meshi].
                    objectRegistry::lookupObjectRef<fldType>(name0);
                meshToField.set(meshi, &meshFld);
            }
        }

        MapVolField
        (
            meshToField,
            oldPatchStarts0,
            oldPatchSizes0,
            patchProcAddressing,
            cellProcAddressing,
            faceProcAddressing,
            fullyMapped
        );
    }
}


template<class Type>
void Foam::fvMeshAdder::MapDimFields
(
    const UPtrList<fvMesh>& meshes,
    const labelListList& cellProcAddressing,
    const bool fullyMapped
)
{
    typedef DimensionedField<Type, volMesh> fldType;
    typedef GeometricField<Type, fvPatchField, volMesh> excludeType;

    if (meshes.size() == 0 || !meshes.set(0))
    {
        FatalErrorInFunction << "Not valid field at element 0"
            << " in mesh list of size " << meshes.size() << exit(FatalError);
    }
    const fvMesh& mesh0 = meshes[0];

    HashTable<const fldType*> fields
    (
        mesh0.objectRegistry::lookupClass<fldType>()
    );


    for (const auto& fld : fields)
    {
        if (!isA<excludeType>(*fld))
        {
            const word& name0 = fld->name();

            DebugPout
                << "MapDimFields : mapping " << name0 << endl;

            UPtrList<fldType> meshToField(meshes.size());
            forAll(meshes, meshi)
            {
                if (meshes.set(meshi))
                {
                    auto& meshFld = meshes[meshi].
                        objectRegistry::lookupObjectRef<fldType>(name0);
                    meshToField.set(meshi, &meshFld);
                }
            }

            MapDimField(meshToField, cellProcAddressing, fullyMapped);
        }
        else
        {
            DebugPout
                << "MapDimFields : ignoring " << fld->name() << endl;
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapSurfaceFields
(
    const UPtrList<fvMesh>& meshes,
    const labelList& oldFaceOwner0,
    const labelList& oldPatchStarts0,
    const labelList& oldPatchSizes0,
    const labelListList& patchProcAddressing,
    const labelListList& cellProcAddressing,
    const labelListList& faceProcAddressing,
    const labelListList& pointProcAddressing,
    const bool fullyMapped
)
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fldType;

    if (meshes.size() == 0 || !meshes.set(0))
    {
        FatalErrorInFunction << "Not valid field at element 0"
            << " in mesh list of size " << meshes.size() << exit(FatalError);
    }
    const auto& mesh0 = meshes[0];

    HashTable<const fldType*> fields
    (
        mesh0.objectRegistry::lookupClass<fldType>()
    );


    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    for (const auto& fld : fields)
    {
        DebugPout
            << "MapSurfaceFields : Storing old time for " << fld->name()
            << endl;

        const_cast<fldType&>(*fld).storeOldTimes();
    }


    for (const auto& fld : fields)
    {
        const word& name0 = fld->name();

        DebugPout
            << "MapSurfaceFields : Mapping " << fld->name() << endl;

        UPtrList<fldType> meshToField(meshes.size());
        forAll(meshes, meshi)
        {
            if (meshes.set(meshi))
            {
                auto& meshFld = meshes[meshi].
                    objectRegistry::lookupObjectRef<fldType>(name0);
                meshToField.set(meshi, &meshFld);
            }
        }

        MapSurfaceField
        (
            meshToField,
            oldFaceOwner0,
            oldPatchStarts0,
            oldPatchSizes0,
            patchProcAddressing,
            cellProcAddressing,
            faceProcAddressing,
            fullyMapped
        );
    }
}


// ************************************************************************* //
