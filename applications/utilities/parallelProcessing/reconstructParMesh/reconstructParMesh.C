/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    reconstructParMesh

Description
    Reconstructs a mesh using geometric information only. Writes
    point/face/cell procAddressing so afterwards reconstructPar can be used to
    reconstruct fields.

    Note:
    - uses geometric matching tolerance (set with -mergeTol option)

    If the parallel case does not have correct procBoundaries use the
    -fullMatch option which will check all boundary faces (bit slower).

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOobjectList.H"
#include "labelIOList.H"
#include "processorPolyPatch.H"
#include "mapAddedPolyMesh.H"
#include "polyMeshAdder.H"
#include "faceCoupleInfo.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Tolerance (as fraction of the bounding box). Needs to be fairly lax since
// usually meshes get written with limited precision (6 digits)
static const scalar defaultMergeTol = 1E-7;


static void renumber
(
    const labelList& map,
    labelList& elems
)
{
    forAll(elems, i)
    {
        if (elems[i] >= 0)
        {
            elems[i] = map[elems[i]];
        }
    }
}


// Determine which faces are coupled. Uses geometric merge distance.
// Looks either at all boundaryFaces (fullMatch) or only at the
// procBoundaries for procI. Assumes that masterMesh contains already merged
// all the processors < procI.
autoPtr<faceCoupleInfo> determineCoupledFaces
(
    const bool fullMatch,
    const label procI,
    const polyMesh& masterMesh,
    const polyMesh& meshToAdd,
    const scalar mergeDist
)
{
    if (fullMatch || masterMesh.nCells() == 0)
    {
        return autoPtr<faceCoupleInfo>
        (
            new faceCoupleInfo
            (
                masterMesh,
                meshToAdd,
                mergeDist,      // absolute merging distance
                true            // matching faces identical
            )
        );
    }
    else
    {
        // Pick up all patches on masterMesh ending in "toDDD" where DDD is
        // the processor number procI.

        const polyBoundaryMesh& masterPatches = masterMesh.boundaryMesh();

        const string toProcString("to" + name(procI));

        DynamicList<label> masterFaces
        (
            masterMesh.nFaces()
          - masterMesh.nInternalFaces()
        );

        forAll(masterPatches, patchI)
        {
            const polyPatch& pp = masterPatches[patchI];

            if
            (
                isA<processorPolyPatch>(pp)
             && (
                    pp.name().rfind(toProcString)
                 == (pp.name().size()-toProcString.size())
                )
            )
            {
                label meshFaceI = pp.start();
                forAll(pp, i)
                {
                    masterFaces.append(meshFaceI++);
                }
            }
        }
        masterFaces.shrink();


        // Pick up all patches on meshToAdd ending in "procBoundaryDDDtoYYY"
        // where DDD is the processor number procI and YYY is < procI.

        const polyBoundaryMesh& addPatches = meshToAdd.boundaryMesh();

        DynamicList<label> addFaces
        (
            meshToAdd.nFaces()
          - meshToAdd.nInternalFaces()
        );

        forAll(addPatches, patchI)
        {
            const polyPatch& pp = addPatches[patchI];

            if (isA<processorPolyPatch>(pp))
            {
                bool isConnected = false;

                for (label mergedProcI = 0; mergedProcI < procI; mergedProcI++)
                {
                    const string fromProcString
                    (
                        "procBoundary"
                      + name(procI)
                      + "to"
                      + name(mergedProcI)
                    );

                    if (pp.name() == fromProcString)
                    {
                        isConnected = true;
                        break;
                    }
                }

                if (isConnected)
                {
                    label meshFaceI = pp.start();
                    forAll(pp, i)
                    {
                        addFaces.append(meshFaceI++);
                    }
                }
            }
        }
        addFaces.shrink();

        return autoPtr<faceCoupleInfo>
        (
            new faceCoupleInfo
            (
                masterMesh,
                masterFaces,
                meshToAdd,
                addFaces,
                mergeDist,      // absolute merging distance
                true,           // matching faces identical?
                false,          // if perfectmatch are faces already ordered
                                // (e.g. processor patches)
                false           // are faces each on separate patch?
            )
        );
    }
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("mergeTol", "relative merge distance");
    argList::validOptions.insert("fullMatch", "");

#   include "addTimeOptions.H"
#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"

    Pout<< "This is an experimental tool which tries to merge"
        << " individual processor" << nl
        << "meshes back into one master mesh. Use it if the original"
        << " master mesh has" << nl
        << "been deleted or if the processor meshes have been modified"
        << " (topology change)." << nl
        << "This tool will write the resulting mesh to a new time step"
        << " and construct" << nl
        << "xxxxProcAddressing files in the processor meshes so"
        << " reconstructPar can be" << nl
        << "used to regenerate the fields on the master mesh." << nl
        << nl
        << "Not well tested & use at your own risk!" << nl
        << endl;


    word regionName = polyMesh::defaultRegion;
    fileName regionPrefix = "";
    if (args.options().found("region"))
    {
        regionName = args.options()["region"];
        regionPrefix = regionName;
        Info<< "Operating on region " << regionName << nl << endl;
    }

    scalar mergeTol = defaultMergeTol;
    if (args.options().found("mergeTol"))
    {
        mergeTol = readScalar(IStringStream(args.options()["mergeTol"])());
    }
    scalar writeTol = Foam::pow(10.0, -scalar(IOstream::defaultPrecision()));

    Pout<< "Merge tolerance : " << mergeTol << nl
        << "Write tolerance : " << writeTol << endl;

    if (runTime.writeFormat() == IOstream::ASCII && mergeTol < writeTol)
    {
        FatalErrorIn(args.executable())
            << "Your current settings specify ASCII writing with "
            << IOstream::defaultPrecision() << " digits precision." << endl
            << "Your merging tolerance (" << mergeTol << ") is finer than this."
            << endl
            << "Please change your writeFormat to binary"
            << " or increase the writePrecision" << endl
            << "or adjust the merge tolerance (-mergeTol)."
            << exit(FatalError);
    }


    const bool fullMatch = args.options().found("fullMatch");

    if (fullMatch)
    {
        Pout<< "Doing geometric matching on all boundary faces." << nl << endl;
    }
    else
    {
        Pout<< "Doing geometric matching on correct procBoundaries only."
            << nl << "This assumes a correct decomposition." << endl;
    }



    int nProcs = 0;

    while
    (
        exists
        (
            args.rootPath()
          / args.caseName()
          / fileName(word("processor") + name(nProcs))
        )
    )
    {
        nProcs++;
    }

    Pout<< "Found " << nProcs << " processor directories" << nl << endl;


    // Read all databases.
    PtrList<Time> databases(nProcs);

    forAll (databases, procI)
    {
        Pout<< "Reading database "
            << args.caseName()/fileName(word("processor") + name(procI))
            << endl;

        databases.set
        (
            procI,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/fileName(word("processor") + name(procI))
            )
        );

        Time& procTime = databases[procI];

        instantList Times = procTime.times();

        // set startTime and endTime depending on -time and -latestTime options
#       include "checkTimeOptions.H"

        procTime.setTime(Times[startTime], startTime);

        if (procI > 0 && databases[procI-1].value() != procTime.value())
        {
            FatalErrorIn(args.executable())
                << "Time not equal on processors." << nl
                << "Processor:" << procI-1
                << " time:" << databases[procI-1].value() << nl
                << "Processor:" << procI
                << " time:" << procTime.value()
                << exit(FatalError);
        }
    }

    // Set master time
    Pout<< "Setting master time to " << databases[0].timeName() << nl << endl;
    runTime.setTime(databases[0]);


    // Read point on individual processors to determine merge tolerance
    // (otherwise single cell domains might give problems)

    boundBox bb
    (
        point(GREAT, GREAT, GREAT),
        point(-GREAT, -GREAT, -GREAT)
    );

    for (label procI = 0; procI < nProcs; procI++)
    {
        fileName pointsInstance
        (
            databases[procI].findInstance
            (
                regionPrefix/polyMesh::meshSubDir,
                "points"
            )
        );

        if (pointsInstance != databases[procI].timeName())
        {
            FatalErrorIn(args.executable())
                << "Your time was specified as " << databases[procI].timeName()
                << " but there is no polyMesh/points in that time." << endl
                << "(there is a points file in " << pointsInstance
                << ")" << endl
                << "Please rerun with the correct time specified"
                << " (through the -constant, -time or -latestTime option)."
                << endl << exit(FatalError);
        }

        Pout<< "Reading points from "
            << databases[procI].caseName()
            << " for time = " << databases[procI].timeName()
            << nl << endl;

        pointIOField points
        (
            IOobject
            (
                "points",
                databases[procI].findInstance
                (
                    regionPrefix/polyMesh::meshSubDir,
                    "points"
                ),
                regionPrefix/polyMesh::meshSubDir,
                databases[procI],
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        boundBox domainBb(points, false);

        bb.min() = min(bb.min(), domainBb.min());
        bb.max() = max(bb.max(), domainBb.max());
    }
    const scalar mergeDist = mergeTol*mag(bb.max() - bb.min());

    Pout<< "Overall mesh bounding box  : " << bb << nl
        << "Relative tolerance         : " << mergeTol << nl
        << "Absolute matching distance : " << mergeDist << nl
        << endl;


    // Addressing from processor to reconstructed case
    labelListList cellProcAddressing(nProcs);
    labelListList faceProcAddressing(nProcs);
    labelListList pointProcAddressing(nProcs);
    labelListList boundaryProcAddressing(nProcs);

    // Internal faces on the final reconstructed mesh
    label masterInternalFaces;
    // Owner addressing on the final reconstructed mesh
    labelList masterOwner;

    {
        // Construct empty mesh.
        Pout<< "Constructing empty mesh to add to." << nl << endl;
        polyMesh masterMesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::NO_READ
            ),
            pointField(0),
            faceList(0),
            cellList(0)
        );

        for (label procI = 0; procI < nProcs; procI++)
        {
            Pout<< "Reading mesh to add from "
                << databases[procI].caseName()
                << " for time = " << databases[procI].timeName()
                << nl << endl;

            polyMesh meshToAdd
            (
                IOobject
                (
                    regionName,
                    databases[procI].timeName(),
                    databases[procI]
                )
            );

            // Initialize its addressing
            cellProcAddressing[procI] = identity(meshToAdd.nCells());
            faceProcAddressing[procI] = identity(meshToAdd.nFaces());
            pointProcAddressing[procI] = identity(meshToAdd.nPoints());
            boundaryProcAddressing[procI] =
                identity(meshToAdd.boundaryMesh().size());


            // Find geometrically shared points/faces.
            autoPtr<faceCoupleInfo> couples =
                determineCoupledFaces
                (
                    fullMatch,
                    procI,
                    masterMesh,
                    meshToAdd,
                    mergeDist
                );


            // Add elements to mesh
            Pout<< "Adding to master mesh" << nl << endl;

            autoPtr<mapAddedPolyMesh> map = polyMeshAdder::add
            (
                masterMesh,
                meshToAdd,
                couples
            );

            // Update all addressing so xxProcAddressing points to correct item
            // in masterMesh.

            // Processors that were already in masterMesh
            for (label mergedI = 0; mergedI < procI; mergedI++)
            {
                renumber(map().oldCellMap(), cellProcAddressing[mergedI]);
                renumber(map().oldFaceMap(), faceProcAddressing[mergedI]);
                renumber(map().oldPointMap(), pointProcAddressing[mergedI]);
                // Note: boundary is special since can contain -1.
                renumber(map().oldPatchMap(), boundaryProcAddressing[mergedI]);
            }

            // Added processor
            renumber(map().addedCellMap(), cellProcAddressing[procI]);
            renumber(map().addedFaceMap(), faceProcAddressing[procI]);
            renumber(map().addedPointMap(), pointProcAddressing[procI]);
            renumber(map().addedPatchMap(), boundaryProcAddressing[procI]);

            Pout<< endl;
        }


        // Save some properties on the reconstructed mesh
        masterInternalFaces = masterMesh.nInternalFaces();
        masterOwner = masterMesh.faceOwner();


        Pout<< "\nWriting merged mesh to "
            << runTime.path()/runTime.timeName()
            << nl << endl;

        if (!masterMesh.write())
        {
            FatalErrorIn(args.executable())
                << "Failed writing polyMesh."
                << exit(FatalError);
        }
    }


    // Write the addressing

    Pout<< "Reconstructing the addressing from the processor meshes"
        << " to the newly reconstructed mesh" << nl << endl;

    forAll(databases, procI)
    {
        Pout<< "Reading processor " << procI << " mesh from "
            << databases[procI].caseName() << endl;

        polyMesh procMesh
        (
            IOobject
            (
                regionName,
                databases[procI].timeName(),
                databases[procI]
            )
        );


        // From processor point to reconstructed mesh point

        Pout<< "Writing pointProcAddressing to "
            << databases[procI].caseName()
              /procMesh.facesInstance()
              /polyMesh::meshSubDir
            << endl;

        labelIOList
        (
            IOobject
            (
                "pointProcAddressing",
                procMesh.facesInstance(),
                polyMesh::meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false                       // do not register
            ),
            pointProcAddressing[procI]
        ).write();


        // From processor face to reconstructed mesh face

        Pout<< "Writing faceProcAddressing to "
            << databases[procI].caseName()
              /procMesh.facesInstance()
              /polyMesh::meshSubDir
            << endl;

        labelIOList faceProcAddr
        (
            IOobject
            (
                "faceProcAddressing",
                procMesh.facesInstance(),
                polyMesh::meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false                       // do not register
            ),
            faceProcAddressing[procI]
        );

        // Now add turning index to faceProcAddressing.
        // See reconstrurPar for meaning of turning index.
        forAll(faceProcAddr, procFaceI)
        {
            label masterFaceI = faceProcAddr[procFaceI];

            if
            (
               !procMesh.isInternalFace(procFaceI)
             && masterFaceI < masterInternalFaces
            )
            {
                // proc face is now external but used to be internal face.
                // Check if we have owner or neighbour.

                label procOwn = procMesh.faceOwner()[procFaceI];
                label masterOwn = masterOwner[masterFaceI];

                if (cellProcAddressing[procI][procOwn] == masterOwn)
                {
                    // No turning. Offset by 1.
                    faceProcAddr[procFaceI]++;
                }
                else
                {
                    // Turned face.
                    faceProcAddr[procFaceI] =
                        -1 - faceProcAddr[procFaceI];
                }
            }
            else
            {
                // No turning. Offset by 1.
                faceProcAddr[procFaceI]++;
            }
        }

        faceProcAddr.write();


        // From processor cell to reconstructed mesh cell

        Pout<< "Writing cellProcAddressing to "
            << databases[procI].caseName()
              /procMesh.facesInstance()
              /polyMesh::meshSubDir
            << endl;

        labelIOList
        (
            IOobject
            (
                "cellProcAddressing",
                procMesh.facesInstance(),
                polyMesh::meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false                       // do not register
            ),
            cellProcAddressing[procI]
        ).write();



        // From processor patch to reconstructed mesh patch

        Pout<< "Writing boundaryProcAddressing to "
            << databases[procI].caseName()
              /procMesh.facesInstance()
              /polyMesh::meshSubDir
            << endl;

        labelIOList
        (
            IOobject
            (
                "boundaryProcAddressing",
                procMesh.facesInstance(),
                polyMesh::meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false                       // do not register
            ),
            boundaryProcAddressing[procI]
        ).write();

        Pout<< endl;
    }

    Pout<< "End.\n" << endl;

    return 0;
}


// ************************************************************************* //
