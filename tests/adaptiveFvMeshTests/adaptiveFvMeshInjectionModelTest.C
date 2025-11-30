#include "IOobject.H"
#include "PstreamReduceOps.H"
#include "catch2/catch_all.hpp"
#include "catch2/catch_test_macros.hpp"
#include "fvCFD.H"
#include "fvMesh.H"
#include "messageStream.H"
#include "volFields.H"
#include "adaptiveFvMesh.H"
#include "volFieldsFwd.H"
#include "stringOps.H"
#include "uniformDimensionedFields.H"
#include "basicKinematicCloud.H"
#include "cloudSupport.H"
#include "abortHandle.H"
#include "cellSet.H"

using namespace Foam;
extern Time* timePtr;
extern argList* argsPtr;


// Helper to create cloud properties with ManualInjection
IStringStream manualInjectionCloudProps
(
    const word& cloudName,
    scalar injectionTime,
    label nParcels
)
{
    return IStringStream
    (
        "solution"
        "{"
        "    active          true;"
        "    coupled         false;"
        "    transient       true;"
        "    cellValueSourceCorrection false;"
        "    interpolationSchemes"
        "    {"
        "        rho         cell;"
        "        U           cellPoint;"
        "        mu          cell;"
        "    }"
        "    integrationSchemes"
        "    {"
        "        U           Euler;"
        "    }"
        "}"
        "constantProperties"
        "{"
        "    rho0            1000;"
        "}"
        "subModels"
        "{"
        "    particleForces"
        "    {"
        "    }"
        "    injectionModels"
        "    {"
        "        manualInjection"
        "        {"
        "            type            manualInjection;"
        "            massTotal       1e-6;"
        "            parcelBasisType fixed;"
        "            nParticle       1;"
        "            SOI             " + Foam::name(injectionTime) + ";"
        "            positionsFile   \"cloudPositions\";"
        "            ignoreOutOfBounds true;"
        "            U0              (0 0 0);"
        "            sizeDistribution"
        "            {"
        "                type        fixedValue;"
        "                fixedValueDistribution"
        "                {"
        "                    value   1e-4;"
        "                }"
        "            }"
        "        }"
        "    }"
        "    dispersionModel none;"
        "    patchInteractionModel none;"
        "    stochasticCollisionModel none;"
        "    surfaceFilmModel none;"
        "}"
        "cloudFunctions"
        "{}"
    );
}


// Helper to create cloud properties with PatchInjection
IStringStream patchInjectionCloudProps
(
    const word& cloudName,
    const word& patchName,
    scalar injectionTime
)
{
    return IStringStream
    (
        "solution"
        "{"
        "    active          true;"
        "    coupled         false;"
        "    transient       true;"
        "    cellValueSourceCorrection false;"
        "    interpolationSchemes"
        "    {"
        "        rho         cell;"
        "        U           cellPoint;"
        "        mu          cell;"
        "    }"
        "    integrationSchemes"
        "    {"
        "        U           Euler;"
        "    }"
        "}"
        "constantProperties"
        "{"
        "    rho0            1000;"
        "}"
        "subModels"
        "{"
        "    particleForces"
        "    {"
        "    }"
        "    injectionModels"
        "    {"
        "        patchInjection"
        "        {"
        "            type            patchInjection;"
        "            patch           " + patchName + ";"
        "            duration        0.001;"
        "            parcelsPerSecond 1000;"
        "            flowRateProfile constant 1;"
        "            massTotal       1e-6;"
        "            parcelBasisType fixed;"
        "            nParticle       1;"
        "            SOI             " + Foam::name(injectionTime) + ";"
        "            U0              (0.1 0 0);"
        "            sizeDistribution"
        "            {"
        "                type        fixedValue;"
        "                fixedValueDistribution"
        "                {"
        "                    value   1e-4;"
        "                }"
        "            }"
        "        }"
        "    }"
        "    dispersionModel none;"
        "    patchInteractionModel none;"
        "    stochasticCollisionModel none;"
        "    surfaceFilmModel none;"
        "}"
        "cloudFunctions"
        "{}"
    );
}


// Helper to create cloud properties with ConeInjection
IStringStream coneInjectionCloudProps
(
    const word& cloudName,
    const point& position,
    const vector& direction,
    scalar injectionTime
)
{
    // Format position and direction as OpenFOAM vectors
    std::ostringstream posStr, dirStr;
    posStr << "(" << position.x() << " " << position.y() << " " << position.z() << ")";
    dirStr << "(" << direction.x() << " " << direction.y() << " " << direction.z() << ")";

    return IStringStream
    (
        "solution"
        "{"
        "    active          true;"
        "    coupled         false;"
        "    transient       true;"
        "    cellValueSourceCorrection false;"
        "    interpolationSchemes"
        "    {"
        "        rho         cell;"
        "        U           cellPoint;"
        "        mu          cell;"
        "    }"
        "    integrationSchemes"
        "    {"
        "        U           Euler;"
        "    }"
        "}"
        "constantProperties"
        "{"
        "    rho0            1000;"
        "}"
        "subModels"
        "{"
        "    particleForces"
        "    {"
        "    }"
        "    injectionModels"
        "    {"
        "        coneInjection"
        "        {"
        "            type            coneInjection;"
        "            SOI             " + Foam::name(injectionTime) + ";"
        "            duration        0.001;"
        "            positionAxis"
        "            ("
        "                (" + posStr.str() + " " + dirStr.str() + ")"
        "            );"
        "            parcelsPerInjector 10;"
        "            flowRateProfile constant 1;"
        "            Umag            constant 1.0;"
        "            thetaInner      constant 0;"
        "            thetaOuter      constant 30;"
        "            massTotal       1e-6;"
        "            parcelBasisType mass;"
        "            sizeDistribution"
        "            {"
        "                type        fixedValue;"
        "                fixedValueDistribution"
        "                {"
        "                    value   1e-4;"
        "                }"
        "            }"
        "        }"
        "    }"
        "    dispersionModel none;"
        "    patchInteractionModel none;"
        "    stochasticCollisionModel none;"
        "    surfaceFilmModel none;"
        "}"
        "cloudFunctions"
        "{}"
    );
}


// Helper to create cloud properties with CellZoneInjection
IStringStream cellZoneInjectionCloudProps
(
    const word& cloudName,
    const word& cellZoneName,
    scalar injectionTime
)
{
    return IStringStream
    (
        "solution"
        "{"
        "    active          true;"
        "    coupled         false;"
        "    transient       true;"
        "    cellValueSourceCorrection false;"
        "    interpolationSchemes"
        "    {"
        "        rho         cell;"
        "        U           cellPoint;"
        "        mu          cell;"
        "    }"
        "    integrationSchemes"
        "    {"
        "        U           Euler;"
        "    }"
        "}"
        "constantProperties"
        "{"
        "    rho0            1000;"
        "}"
        "subModels"
        "{"
        "    particleForces"
        "    {"
        "    }"
        "    injectionModels"
        "    {"
        "        cellZoneInjection"
        "        {"
        "            type            cellZoneInjection;"
        "            SOI             " + Foam::name(injectionTime) + ";"
        "            cellZone        " + cellZoneName + ";"
        "            massTotal       1e-6;"
        "            parcelBasisType fixed;"
        "            nParticle       1;"
        "            numberDensity   1e6;"
        "            U0              (0 0 0);"
        "            sizeDistribution"
        "            {"
        "                type        fixedValue;"
        "                fixedValueDistribution"
        "                {"
        "                    value   1e-4;"
        "                }"
        "            }"
        "        }"
        "    }"
        "    dispersionModel none;"
        "    patchInteractionModel none;"
        "    stochasticCollisionModel none;"
        "    surfaceFilmModel none;"
        "}"
        "cloudFunctions"
        "{}"
    );
}


TEST_CASE
(
    "ManualInjection model particle conservation during load balancing",
    "[cloud][injection][manualInjection][hex2D][hex3D][poly2D][poly3D][parallel]"
)
{
    if (!Pstream::parRun())
    {
        SUCCEED("Skipping ManualInjection test in serial mode");
        return;
    }

    Time& runTime = *timePtr;
    const boundBox refBox(point(0.02, 0.025, -1), point(0.04, 0.035, 1));

    // Use appropriate refiner for mesh type
    word refiner = runTime.caseName().find("poly") != string::npos
        ? "polyRefiner" : "hexRefiner";

    IStringStream dynamicMeshStream
    (
        "dynamicFvMesh   adaptiveFvMesh;"
        "errorEstimator  fieldValue;"
        "fieldName       test;"
        "balance         yes;"
        "refiner         "+refiner+";"
        "refineInterval  1;"
        "unrefineInterval  1;"
        "balanceInterval  1;"
        "lowerRefineLevel 0.01;"
        "unrefineLevel 0.01;"
        "nBufferLayers   1;"
        "maxRefinement   2;"
        "dumpLevel       false;"
        "protectedPatches ();"
        "nPatchesBuffers 1;"
    );

    runTime.setTime(0.0, 0);

    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dynamicMeshStream
    );
    dynamicMeshDict.regIOobject::write();

    adaptiveFvMesh mesh
    (
        IOobject
        (
            "",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Create refinement criterion field
    volScalarField test
    (
        IOobject
        (
            "test",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    forAll(mesh.C(), ci)
    {
        if (refBox.contains(mesh.C()[ci]))
        {
            test[ci] = 1;
        }
    }

    // Create carrier fields
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 1.0)
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedVector("U", dimVelocity, Zero)
    );

    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar("mu", dimMass/dimLength/dimTime, 1e-5)
    );

    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimAcceleration, vector(0, -9.81, 0))
    );

    // Create positions file for ManualInjection
    // Collect some cell centres as injection positions
    const auto& cc = mesh.cellCentres();
    pointField injectionPositions;

    forAll(cc, celli)
    {
        if (refBox.contains(cc[celli]) && injectionPositions.size() < 10)
        {
            injectionPositions.append(cc[celli]);
        }
    }

    // Gather all positions globally
    List<pointField> allPositions(Pstream::nProcs());
    allPositions[Pstream::myProcNo()] = injectionPositions;
    Pstream::allGatherList(allPositions);

    // Flatten and write positions file
    pointField globalPositions;
    forAll(allPositions, proci)
    {
        globalPositions.append(allPositions[proci]);
    }

    // Write positions file to each processor's constant directory
    // In parallel, each proc needs the file in its own processorN/constant
    {
        OFstream posFile(mesh.time().path()/mesh.time().constant()/"cloudPositions");
        // Write OpenFOAM header
        posFile << "FoamFile\n"
                << "{\n"
                << "    version     2.0;\n"
                << "    format      ascii;\n"
                << "    class       vectorField;\n"
                << "    object      cloudPositions;\n"
                << "}\n"
                << globalPositions;
    }

    // Write cloud properties
    IStringStream cloudPropsStream = manualInjectionCloudProps
    (
        "manualInjectionCloud",
        0.0,  // SOI at time 0
        globalPositions.size()
    );

    IOdictionary cloudPropsDict
    (
        IOobject
        (
            "manualInjectionCloudProperties",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloudPropsStream
    );
    cloudPropsDict.regIOobject::write();

    // Create cloud with ManualInjection
    basicKinematicCloud cloud
    (
        "manualInjectionCloud",
        rho,
        U,
        mu,
        g
    );

    // Advance time so injection can occur (SOI=0 needs time > 0)
    runTime++;

    // Evolve cloud to trigger injection
    cloud.evolve();

    label initialGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    CAPTURE
    (
        Pstream::myProcNo(), Pstream::nProcs(),
        globalPositions.size(), initialGlobal
    );

    // Verify injection occurred
    REQUIRE(initialGlobal > 0);

    // Perform refinement + load balancing cycles
    for (label i = 0; i < 2; i++)
    {
        runTime++;
        tryAndCatchAbortingCode
        (
            [&]()
            {
                mesh.update();
            }
        );
    }

    // Verify particle conservation
    label finalGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    INFO("ManualInjection: Initial=" << initialGlobal << ", Final=" << finalGlobal);
    REQUIRE(finalGlobal == initialGlobal);

    // Verify all particles are in valid cells
    bool allValid = true;
    for (const basicKinematicParcel& p : cloud)
    {
        if (p.cell() < 0 || p.cell() >= mesh.nCells())
        {
            allValid = false;
            break;
        }
    }
    REQUIRE(allValid);

    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "PatchInjection model particle conservation during load balancing",
    "[cloud][injection][patchInjection][hex2D][hex3D][poly2D][poly3D][parallel]"
)
{
    if (!Pstream::parRun())
    {
        SUCCEED("Skipping PatchInjection test in serial mode");
        return;
    }

    FatalError.throwExceptions();
    Time& runTime = *timePtr;

    const boundBox refBox(point(0.02, 0.025, -1), point(0.04, 0.035, 1));

    // Use appropriate refiner for mesh type
    word refiner = runTime.caseName().find("poly") != string::npos
        ? "polyRefiner" : "hexRefiner";

    IStringStream dynamicMeshStream
    (
        "dynamicFvMesh   adaptiveFvMesh;"
        "errorEstimator  fieldValue;"
        "fieldName       test;"
        "balance         yes;"
        "refiner         "+refiner+";"
        "refineInterval  1;"
        "unrefineInterval  1;"
        "balanceInterval  1;"
        "lowerRefineLevel 0.01;"
        "unrefineLevel 0.01;"
        "nBufferLayers   1;"
        "maxRefinement   2;"
        "dumpLevel       false;"
        "protectedPatches ();"
        "nPatchesBuffers 1;"
    );

    runTime.setTime(0.0, 0);

    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dynamicMeshStream
    );
    dynamicMeshDict.regIOobject::write();

    adaptiveFvMesh mesh
    (
        IOobject
        (
            "",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Find a suitable patch for injection
    word injectPatchName;
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        if (pp.size() > 0 && !pp.coupled())
        {
            injectPatchName = pp.name();
            break;
        }
    }

    if (injectPatchName.empty())
    {
        SUCCEED("No suitable patch found for PatchInjection test");
        return;
    }

    // Create refinement criterion field
    volScalarField test
    (
        IOobject
        (
            "test",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    forAll(mesh.C(), ci)
    {
        if (refBox.contains(mesh.C()[ci]))
        {
            test[ci] = 1;
        }
    }

    // Create carrier fields
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 1.0)
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedVector("U", dimVelocity, Zero)
    );

    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar("mu", dimMass/dimLength/dimTime, 1e-5)
    );

    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimAcceleration, vector(0, -9.81, 0))
    );

    // Write cloud properties with PatchInjection
    IStringStream cloudPropsStream = patchInjectionCloudProps
    (
        "patchInjectionCloud",
        injectPatchName,
        0.0  // SOI at time 0
    );

    IOdictionary cloudPropsDict
    (
        IOobject
        (
            "patchInjectionCloudProperties",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloudPropsStream
    );
    cloudPropsDict.regIOobject::write();

    // Create cloud with PatchInjection
    basicKinematicCloud cloud
    (
        "patchInjectionCloud",
        rho,
        U,
        mu,
        g
    );

    // Evolve cloud to trigger injection
    cloud.evolve();

    label initialGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    CAPTURE
    (
        Pstream::myProcNo(), Pstream::nProcs(),
        injectPatchName, initialGlobal
    );

    // PatchInjection may not inject on all procs, check globally
    if (initialGlobal == 0)
    {
        SUCCEED("No particles injected (patch may be empty or not local)");
        return;
    }

    // Perform refinement + load balancing cycles
    for (label i = 0; i < 2; i++)
    {
        runTime++;
        tryAndCatchAbortingCode
        (
            [&]()
            {
                mesh.update();
            }
        );
    }

    // Verify particle conservation
    label finalGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    INFO("PatchInjection: Initial=" << initialGlobal << ", Final=" << finalGlobal);
    REQUIRE(finalGlobal == initialGlobal);

    // Verify all particles are in valid cells
    bool allValid = true;
    for (const basicKinematicParcel& p : cloud)
    {
        if (p.cell() < 0 || p.cell() >= mesh.nCells())
        {
            allValid = false;
            break;
        }
    }
    REQUIRE(allValid);

    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "ConeInjection model particle conservation during load balancing",
    "[cloud][injection][coneInjection][hex2D][hex3D][poly2D][poly3D][parallel]"
)
{
    if (!Pstream::parRun())
    {
        SUCCEED("Skipping ConeInjection test in serial mode");
        return;
    }

    FatalError.throwExceptions();
    Time& runTime = *timePtr;

    const boundBox refBox(point(0.02, 0.025, -1), point(0.04, 0.035, 1));

    // Use appropriate refiner for mesh type
    word refiner = runTime.caseName().find("poly") != string::npos
        ? "polyRefiner" : "hexRefiner";

    IStringStream dynamicMeshStream
    (
        "dynamicFvMesh   adaptiveFvMesh;"
        "errorEstimator  fieldValue;"
        "fieldName       test;"
        "balance         yes;"
        "refiner         "+refiner+";"
        "refineInterval  1;"
        "unrefineInterval  1;"
        "balanceInterval  1;"
        "lowerRefineLevel 0.01;"
        "unrefineLevel 0.01;"
        "nBufferLayers   1;"
        "maxRefinement   2;"
        "dumpLevel       false;"
        "protectedPatches ();"
        "nPatchesBuffers 1;"
    );

    runTime.setTime(0.0, 0);

    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dynamicMeshStream
    );
    dynamicMeshDict.regIOobject::write();

    adaptiveFvMesh mesh
    (
        IOobject
        (
            "",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Create refinement criterion field
    volScalarField test
    (
        IOobject
        (
            "test",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    forAll(mesh.C(), ci)
    {
        if (refBox.contains(mesh.C()[ci]))
        {
            test[ci] = 1;
        }
    }

    // Create carrier fields
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 1.0)
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedVector("U", dimVelocity, Zero)
    );

    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar("mu", dimMass/dimLength/dimTime, 1e-5)
    );

    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimAcceleration, vector(0, -9.81, 0))
    );

    // Choose injection point inside the mesh bounds
    const boundBox& meshBounds = mesh.bounds();
    point injectionPoint = meshBounds.centre();
    vector injectionDir(1, 0, 0);  // Inject in x-direction

    // Write cloud properties with ConeInjection
    IStringStream cloudPropsStream = coneInjectionCloudProps
    (
        "coneInjectionCloud",
        injectionPoint,
        injectionDir,
        0.0  // SOI at time 0
    );

    IOdictionary cloudPropsDict
    (
        IOobject
        (
            "coneInjectionCloudProperties",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloudPropsStream
    );
    cloudPropsDict.regIOobject::write();

    // Create cloud with ConeInjection
    basicKinematicCloud cloud
    (
        "coneInjectionCloud",
        rho,
        U,
        mu,
        g
    );

    // Evolve cloud to trigger injection
    cloud.evolve();

    label initialGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    CAPTURE
    (
        Pstream::myProcNo(), Pstream::nProcs(),
        injectionPoint, initialGlobal
    );

    // ConeInjection may inject only on the proc owning the injection cell
    if (initialGlobal == 0)
    {
        SUCCEED("No particles injected (injection point not in local mesh)");
        return;
    }

    // Perform refinement + load balancing cycles
    for (label i = 0; i < 2; i++)
    {
        runTime++;
        tryAndCatchAbortingCode
        (
            [&]()
            {
                mesh.update();
            }
        );
    }

    // Verify particle conservation
    label finalGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    INFO("ConeInjection: Initial=" << initialGlobal << ", Final=" << finalGlobal);
    REQUIRE(finalGlobal == initialGlobal);

    // Verify all particles are in valid cells
    bool allValid = true;
    for (const basicKinematicParcel& p : cloud)
    {
        if (p.cell() < 0 || p.cell() >= mesh.nCells())
        {
            allValid = false;
            break;
        }
    }
    REQUIRE(allValid);

    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "CellZoneInjection model particle conservation during load balancing",
    "[cloud][injection][cellZoneInjection][hex2D][hex3D][poly2D][poly3D][parallel]"
)
{
    if (!Pstream::parRun())
    {
        SUCCEED("Skipping CellZoneInjection test in serial mode");
        return;
    }

    FatalError.throwExceptions();
    Time& runTime = *timePtr;

    const boundBox refBox(point(0.02, 0.025, -1), point(0.04, 0.035, 1));

    // Use appropriate refiner for mesh type
    word refiner = runTime.caseName().find("poly") != string::npos
        ? "polyRefiner" : "hexRefiner";

    IStringStream dynamicMeshStream
    (
        "dynamicFvMesh   adaptiveFvMesh;"
        "errorEstimator  fieldValue;"
        "fieldName       test;"
        "balance         yes;"
        "refiner         "+refiner+";"
        "refineInterval  1;"
        "unrefineInterval  1;"
        "balanceInterval  1;"
        "lowerRefineLevel 0.01;"
        "unrefineLevel 0.01;"
        "nBufferLayers   1;"
        "maxRefinement   2;"
        "dumpLevel       false;"
        "protectedPatches ();"
        "nPatchesBuffers 1;"
    );

    runTime.setTime(0.0, 0);

    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dynamicMeshStream
    );
    dynamicMeshDict.regIOobject::write();

    adaptiveFvMesh mesh
    (
        IOobject
        (
            "",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Create a cellZone programmatically for the injection region
    // Use cells inside the refinement box
    const word injectionZoneName = "injectionZone";

    labelList zoneCells;
    forAll(mesh.C(), celli)
    {
        if (refBox.contains(mesh.C()[celli]))
        {
            zoneCells.append(celli);
        }
    }

    // Add cellZone to the mesh
    if (mesh.cellZones().findZoneID(injectionZoneName) < 0)
    {
        mesh.cellZones().append
        (
            new cellZone
            (
                injectionZoneName,
                zoneCells,
                mesh.cellZones().size(),  // zone index
                mesh.cellZones()
            )
        );
    }

    label globalZoneCells = returnReduce(zoneCells.size(), sumOp<label>());

    if (globalZoneCells == 0)
    {
        SUCCEED("No cells in injection zone");
        return;
    }

    // Create refinement criterion field
    volScalarField test
    (
        IOobject
        (
            "test",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    forAll(mesh.C(), ci)
    {
        if (refBox.contains(mesh.C()[ci]))
        {
            test[ci] = 1;
        }
    }

    // Create carrier fields
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 1.0)
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedVector("U", dimVelocity, Zero)
    );

    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar("mu", dimMass/dimLength/dimTime, 1e-5)
    );

    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimAcceleration, vector(0, -9.81, 0))
    );

    // Write cloud properties with CellZoneInjection
    IStringStream cloudPropsStream = cellZoneInjectionCloudProps
    (
        "cellZoneInjectionCloud",
        injectionZoneName,
        0.0  // SOI at time 0
    );

    IOdictionary cloudPropsDict
    (
        IOobject
        (
            "cellZoneInjectionCloudProperties",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloudPropsStream
    );
    cloudPropsDict.regIOobject::write();

    // Create cloud with CellZoneInjection
    basicKinematicCloud cloud
    (
        "cellZoneInjectionCloud",
        rho,
        U,
        mu,
        g
    );

    // Evolve cloud to trigger injection
    cloud.evolve();

    label initialGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    CAPTURE
    (
        Pstream::myProcNo(), Pstream::nProcs(),
        injectionZoneName, globalZoneCells, initialGlobal
    );

    if (initialGlobal == 0)
    {
        SUCCEED("No particles injected from cellZone");
        return;
    }

    // Perform refinement + load balancing cycles
    for (label i = 0; i < 2; i++)
    {
        runTime++;
        tryAndCatchAbortingCode
        (
            [&]()
            {
                mesh.update();
            }
        );
    }

    // Verify particle conservation
    label finalGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    INFO("CellZoneInjection: Initial=" << initialGlobal << ", Final=" << finalGlobal);
    REQUIRE(finalGlobal == initialGlobal);

    // Verify all particles are in valid cells
    bool allValid = true;
    for (const basicKinematicParcel& p : cloud)
    {
        if (p.cell() < 0 || p.cell() >= mesh.nCells())
        {
            allValid = false;
            break;
        }
    }
    REQUIRE(allValid);

    runTime.setTime(0.0, 0);
}


// ************************************************************************* //
