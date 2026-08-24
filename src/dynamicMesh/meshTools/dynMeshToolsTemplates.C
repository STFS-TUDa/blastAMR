/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dynMeshTools.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "globalMeshData.H"
#include "contiguous.H"
#include "transform.H"
#include "IOobjectList.H"
#include "fvMesh.H"
#include "pointMesh.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Read and add fields to the database
template<class FieldType>
void Foam::meshTools::readGeoFields
(
    const fvMesh& mesh,
    const IOobjectList& objects
)
{

    IOobjectList fields = objects.lookupClass(FieldType::typeName);
    forAllIter(IOobjectList, fields, fieldIter)
    {
        if (!mesh.foundObject<FieldType>(fieldIter()->name()))
        {
            IOobject fieldTargetIOobject
            (
                fieldIter()->name(),
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            );

            if (fieldTargetIOobject.typeHeaderOk<FieldType>(true))
            {
                FieldType* fPtr
                (
                    new FieldType
                    (
                        fieldTargetIOobject,
                        mesh
                    )
                );
                fPtr->store(fPtr);
            }
        }
    }
}


//- Read and add fields to the database
template<class FieldType>
void Foam::meshTools::readPointFields
(
    const fvMesh& mesh,
    const IOobjectList& objects
)
{
    IOobjectList fields(objects.lookupClass(FieldType::typeName));
    forAllIter(IOobjectList, fields, fieldIter)
    {
        if (!mesh.foundObject<FieldType>(fieldIter()->name()))
        {
            IOobject fieldTargetIOobject
            (
                fieldIter()->name(),
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            );

            if (fieldTargetIOobject.typeHeaderOk<FieldType>(true))
            {
                FieldType* fPtr
                (
                    new FieldType
                    (
                        fieldTargetIOobject,
                        pointMesh::New(mesh)
                    )
                );
                fPtr->store(fPtr);
            }
        }
    }
}

template<class T>
void Foam::meshTools::mapNewInternalFaces
(
    fvMesh& mesh,
    const labelList& faceMap,
    GeometricField<T, fvsPatchField, surfaceMesh>& sFld
)
{
    typedef GeometricField<T, fvsPatchField, surfaceMesh> GeoField;

    // Flat field over internal + boundary faces, for ease of looping
    Field<T> tsFld(mesh.nFaces(), Zero);
    SubField<T>(tsFld, mesh.nInternalFaces()) = sFld.primitiveField();

    const typename GeoField::Boundary& bFld = sFld.boundaryField();
    forAll(bFld, patchi)
    {
        label facei = mesh.boundaryMesh()[patchi].start();
        for (const T& val : bFld[patchi])
        {
            tsFld[facei++] = val;
        }
    }

    const labelUList& owner = mesh.faceOwner();
    const labelUList& neighbour = mesh.faceNeighbour();
    const cellList& cells = mesh.cells();

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        if (faceMap[facei] != -1)
        {
            continue;
        }

        // Created out of nothing: average the faces of the owner and
        // neighbour cells that did come from somewhere
        T sum(pTraits<T>::zero);
        label counter = 0;

        for (const label ownFacei : cells[owner[facei]])
        {
            if (faceMap[ownFacei] != -1)
            {
                sum += tsFld[ownFacei];
                counter++;
            }
        }

        for (const label neiFacei : cells[neighbour[facei]])
        {
            if (faceMap[neiFacei] != -1)
            {
                sum += tsFld[neiFacei];
                counter++;
            }
        }

        if (counter > 0)
        {
            sFld[facei] = sum/counter;
        }
    }
}


template<class T>
void Foam::meshTools::mapNewInternalFaces
(
    fvMesh& mesh,
    const labelList& faceMap
)
{
    typedef GeometricField<T, fvsPatchField, surfaceMesh> GeoField;

    HashTable<GeoField*> flds(mesh.objectRegistry::lookupClass<GeoField>());

    forAllIters(flds, iter)
    {
        GeoField& sFld = *iter.val();

        if (sFld.is_oriented())
        {
            continue;
        }

        mapNewInternalFaces<T>(mesh, faceMap, sFld);
    }
}


// ************************************************************************* //
