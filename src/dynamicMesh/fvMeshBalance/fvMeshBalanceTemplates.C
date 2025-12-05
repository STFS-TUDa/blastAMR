/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
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

#include "volMesh.H"
#include "fvPatchField.H"
#include "surfaceFields.H"
#include "processorPolyPatch.H"


template<class GeoField>
void Foam::fvMeshBalance::correctBoundaries()
{
    HashTable<GeoField*> flds(mesh_.lookupClass<GeoField>());

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();

        // Skip if boundary field size doesn't match mesh boundary
        // This can happen during redistribution when patches are added/removed
        if (fld.boundaryField().size() != mesh_.boundaryMesh().size())
        {
            if (debug)
            {
                Pout<< "correctBoundaries: skipping " << fld.name()
                    << " - boundary size mismatch (field: "
                    << fld.boundaryField().size() << ", mesh: "
                    << mesh_.boundaryMesh().size() << ")" << endl;
            }
            continue;
        }

        //mimic "evaluate" but only for coupled patches (processor or cyclic)
        // and only for blocking or nonBlocking comms (no scheduled comms)
        if
        (
            Pstream::defaultCommsType == Pstream::commsTypes::blocking
         || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            label nReq = Pstream::nRequests();

            forAll(fld.boundaryField(), patchi)
            {
                if (patchi < mesh_.boundaryMesh().size()
                    && isA<processorPolyPatch>(mesh_.boundaryMesh()[patchi]))
                {
                    fld.boundaryFieldRef()[patchi].initEvaluate
                    (
                        Pstream::defaultCommsType
                    );
                }
            }

            // Block for any outstanding requests
            if
            (
                Pstream::parRun()
             && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
            )
            {
                Pstream::waitRequests(nReq);
            }

            forAll(fld.boundaryField(), patchi)
            {
                if (patchi < mesh_.boundaryMesh().size()
                    && isA<processorPolyPatch>(mesh_.boundaryMesh()[patchi]))
                {
                    fld.boundaryFieldRef()[patchi].evaluate
                    (
                        Pstream::defaultCommsType
                    );
                }
            }
        }
        else
        {
            //Scheduled patch updates not supported
            FatalErrorInFunction
                << "Unsuported communications type "
                << Pstream::commsTypeNames[Pstream::defaultCommsType]
                << exit(FatalError);
        }
    }
}

template<class GeoField>
void Foam::fvMeshBalance::syncOldTimeFields()
{
    HashTable<GeoField*> flds(mesh_.lookupClass<GeoField>());
    wordList localNames;
    forAllConstIter(typename HashTable<GeoField*>, flds, iter)
    {
        localNames.append(iter.key());
    }
    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = localNames;
    Pstream::allGatherList(allNames);
    wordHashSet allFieldNames;
    forAll(allNames, proci)
    {
        forAll(allNames[proci], i)
        {
            allFieldNames.insert(allNames[proci][i]);
        }
    }

    // For each field in the set, ensure it exists on this processor
    // by triggering oldTime() creation on any field that this processor has
    for (const word& fieldName : allFieldNames)
    {
        GeoField* fldPtr = mesh_.getObjectPtr<GeoField>(fieldName);

        if (!fldPtr)
        {
            // This processor doesn't have this field at all
            // This is pretty serious - the field should exist everywhere
            // For now, only handle oldTime fields (ending with _0)
            if (debug)
            {
                Pout<< "syncOldTimeFields: processor " << Pstream::myProcNo()
                    << " missing field " << fieldName << endl;
            }
            continue;
        }

        // assuming "oldTime" fields end with _0
        if (!fieldName.ends_with("_0"))
        {
            const word oldTimeName = fieldName + "_0";
            if (allFieldNames.found(oldTimeName))
            {
                if (debug)
                {
                    Pout<< "syncOldTimeFields: processor " << Pstream::myProcNo()
                        << " ensuring oldTime for " << fieldName << endl;
                }
                fldPtr->oldTime();
            }
        }
    }
}


template<class Type>
void Foam::fvMeshBalance::pushUntransformedData
(
    const polyMesh& mesh,
    Field<Type>& pointData
)
{
    const globalMeshData& gmd = mesh.globalData();
//     Field<scalar> nSharedPoints(pointData.size(), 1);
//     gmd.syncPointData
//     (
//         pointData,
//         plusEqOp<Type>(),
//         mapDistribute::transform()
//     );
//     gmd.syncPointData
//     (
//         nSharedPoints,
//         plusEqOp<scalar>(),
//         mapDistribute::transform()
//     );
//     pointData /= nSharedPoints;

    // Transfer onto coupled patch

    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, i)
    {
        elems[i] = pointData[meshPoints[i]];
    }

    // Combine master data with slave data
    forAll(slaves, i)
    {
        const labelList& slavePoints = slaves[i];

        // Copy master data to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elems[i];
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    forAll(meshPoints, i)
    {
        pointData[meshPoints[i]] = elems[i];
    }
}

// ************************************************************************* //
