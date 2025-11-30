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
#include "boxToCell.H"
#include "cellSet.H"
#include "stringOps.H"
#include "uniformDimensionedFields.H"
#include "passiveParticleCloud.H"
#include "cloudSupport.H"
#include "basicKinematicCloud.H"
#include "basicKinematicCollidingCloud.H"
#include "basicKinematicMPPICCloud.H"
#include "abortHandle.H"

using namespace Foam;
extern Time* timePtr;
extern argList* argsPtr;


// Helper to create kinematicCloudProperties IStringStream for basic kinematic cloud
IStringStream kinematicCloudPropsStream(const word& cloudName)
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
        "        none"
        "        {"
        "            type    none;"
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

// Helper to create kinematicCloudProperties IStringStream for MPPIC cloud
IStringStream mppicCloudPropsStream(const word& cloudName)
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
        "        none"
        "        {"
        "            type    none;"
        "        }"
        "    }"
        "    dispersionModel none;"
        "    patchInteractionModel none;"
        "    stochasticCollisionModel none;"
        "    surfaceFilmModel none;"
        "    packingModel none;"
        "    dampingModel none;"
        "    isotropyModel none;"
        "}"
        "cloudFunctions"
        "{}"
    );
}

// Helper to create kinematicCloudProperties IStringStream for colliding cloud
IStringStream collidingCloudPropsStream(const word& cloudName)
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
        "        none"
        "        {"
        "            type    none;"
        "        }"
        "    }"
        "    dispersionModel none;"
        "    patchInteractionModel none;"
        "    stochasticCollisionModel none;"
        "    surfaceFilmModel none;"
        "    collisionModel none;"
        "}"
        "cloudFunctions"
        "{}"
    );
}


TEST_CASE
(
    "basicKinematicCloud particle conservation during load balancing",
    "[cloud][kinematicCloud][loadBalance][hex2D][hex3D][poly2D][poly3D][parallel]"
)
{
    if (!Pstream::parRun())
    {
        SUCCEED("Skipping load balancing test in serial mode");
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

    // Create carrier fields required by basicKinematicCloud
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

    // Create gravity
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

    // Write cloud properties using IStringStream pattern
    IStringStream cloudPropsStream = kinematicCloudPropsStream("kinematicCloud");
    IOdictionary cloudPropsDict
    (
        IOobject
        (
            "kinematicCloudProperties",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloudPropsStream
    );
    cloudPropsDict.regIOobject::write();

    // Create basicKinematicCloud
    basicKinematicCloud cloud
    (
        "kinematicCloud",
        rho,
        U,
        mu,
        g
    );

    // Add particles manually across the mesh
    const auto& cc = mesh.cellCentres();
    label targetParticles = 30;
    label added = 0;

    forAll(cc, celli)
    {
        if (added < targetParticles)
        {
            basicKinematicParcel* p = new basicKinematicParcel
            (
                mesh,
                cc[celli],
                celli
            );
            p->d() = 1e-4;
            p->nParticle() = 1.0;
            p->rho() = 1000.0;
            p->U() = Zero;

            cloud.addParticle(p);
            added++;
        }
    }

    label initialGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    CAPTURE
    (
        Pstream::myProcNo(), Pstream::nProcs(),
        runTime.caseName(),
        initialGlobal
    );

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

    // Verify global particle conservation
    label finalGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    INFO("basicKinematicCloud: Initial=" << initialGlobal << ", Final=" << finalGlobal);
    REQUIRE(finalGlobal == initialGlobal);

    // Verify particles are in valid cells
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
    "basicKinematicCollidingCloud particle conservation during load balancing",
    "[cloud][collidingCloud][loadBalance][hex2D][hex3D][poly2D][poly3D][parallel]"
)
{
    if (!Pstream::parRun())
    {
        SUCCEED("Skipping load balancing test in serial mode");
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

    // Write cloud properties with colliding model support
    IStringStream cloudPropsStream = collidingCloudPropsStream("collidingCloud");
    IOdictionary cloudPropsDict
    (
        IOobject
        (
            "collidingCloudProperties",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloudPropsStream
    );
    cloudPropsDict.regIOobject::write();

    // Create basicKinematicCollidingCloud
    basicKinematicCollidingCloud cloud
    (
        "collidingCloud",
        rho,
        U,
        mu,
        g
    );

    // Add particles
    const auto& cc = mesh.cellCentres();
    label targetParticles = 25;
    label added = 0;

    forAll(cc, celli)
    {
        if (added < targetParticles)
        {
            basicKinematicCollidingParcel* p = new basicKinematicCollidingParcel
            (
                mesh,
                cc[celli],
                celli
            );
            p->d() = 1e-4;
            p->nParticle() = 1.0;
            p->rho() = 1000.0;
            p->U() = Zero;

            cloud.addParticle(p);
            added++;
        }
    }

    label initialGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    CAPTURE(Pstream::myProcNo(), Pstream::nProcs(), initialGlobal);
    REQUIRE(initialGlobal > 0);

    // Perform refinement + load balancing
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

    label finalGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    INFO("basicKinematicCollidingCloud: Initial=" << initialGlobal << ", Final=" << finalGlobal);
    REQUIRE(finalGlobal == initialGlobal);

    // Verify valid cells
    bool allValid = true;
    for (const basicKinematicCollidingParcel& p : cloud)
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
    "basicKinematicMPPICCloud particle conservation during load balancing",
    "[cloud][mppicCloud][loadBalance][hex2D][hex3D][poly2D][poly3D][parallel]"
)
{
    if (!Pstream::parRun())
    {
        SUCCEED("Skipping load balancing test in serial mode");
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

    // Write cloud properties with MPPIC model support
    IStringStream cloudPropsStream = mppicCloudPropsStream("mppicCloud");
    IOdictionary cloudPropsDict
    (
        IOobject
        (
            "mppicCloudProperties",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloudPropsStream
    );
    cloudPropsDict.regIOobject::write();

    // Create basicKinematicMPPICCloud
    basicKinematicMPPICCloud cloud
    (
        "mppicCloud",
        rho,
        U,
        mu,
        g
    );

    // Add particles
    const auto& cc = mesh.cellCentres();
    label targetParticles = 35;
    label added = 0;

    forAll(cc, celli)
    {
        if (added < targetParticles)
        {
            basicKinematicMPPICParcel* p = new basicKinematicMPPICParcel
            (
                mesh,
                cc[celli],
                celli
            );
            p->d() = 1e-4;
            p->nParticle() = 1.0;
            p->rho() = 1000.0;
            p->U() = Zero;

            cloud.addParticle(p);
            added++;
        }
    }

    label initialGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    CAPTURE(Pstream::myProcNo(), Pstream::nProcs(), initialGlobal);
    REQUIRE(initialGlobal > 0);

    // Perform refinement + load balancing
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

    label finalGlobal = returnReduce(cloud.nParcels(), sumOp<label>());

    INFO("basicKinematicMPPICCloud: Initial=" << initialGlobal << ", Final=" << finalGlobal);
    REQUIRE(finalGlobal == initialGlobal);

    // Verify valid cells
    bool allValid = true;
    for (const basicKinematicMPPICParcel& p : cloud)
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
    "Multiple cloud types during load balancing",
    "[cloud][multiType][loadBalance][hex2D][hex3D][poly2D][poly3D][parallel]"
)
{
    if (!Pstream::parRun())
    {
        SUCCEED("Skipping multi-type test in serial mode");
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

    // Create shared carrier fields
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

    // Also create a passive cloud for comparison
    passiveParticleCloud passiveCloud(mesh, Foam::zero{}, "passiveCloud");

    // Write cloud properties for kinematic cloud
    IStringStream cloudPropsStream = kinematicCloudPropsStream("kinematicCloud2");
    IOdictionary cloudPropsDict
    (
        IOobject
        (
            "kinematicCloud2Properties",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloudPropsStream
    );
    cloudPropsDict.regIOobject::write();

    basicKinematicCloud kinematicCloud
    (
        "kinematicCloud2",
        rho,
        U,
        mu,
        g
    );

    // Add particles to both clouds
    const auto& cc = mesh.cellCentres();

    // Passive cloud: 20 particles
    label added = 0;
    forAll(cc, celli)
    {
        if (added < 20)
        {
            passiveCloud.addParticle(new passiveParticle(mesh, cc[celli], celli));
            added++;
        }
    }

    // Kinematic cloud: 15 particles in different cells
    added = 0;
    forAll(cc, celli)
    {
        if (celli >= 5 && added < 15)
        {
            basicKinematicParcel* p = new basicKinematicParcel
            (
                mesh,
                cc[celli],
                celli
            );
            p->d() = 1e-4;
            p->nParticle() = 1.0;
            p->rho() = 1000.0;
            p->U() = Zero;

            kinematicCloud.addParticle(p);
            added++;
        }
    }

    label initialPassive = returnReduce(passiveCloud.size(), sumOp<label>());
    label initialKinematic = returnReduce(kinematicCloud.nParcels(), sumOp<label>());

    CAPTURE
    (
        Pstream::myProcNo(), Pstream::nProcs(),
        initialPassive, initialKinematic
    );

    REQUIRE(initialPassive > 0);
    REQUIRE(initialKinematic > 0);

    // Perform refinement + load balancing
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

    // Verify both clouds independently
    label finalPassive = returnReduce(passiveCloud.size(), sumOp<label>());
    label finalKinematic = returnReduce(kinematicCloud.nParcels(), sumOp<label>());

    INFO("passiveCloud: " << initialPassive << " -> " << finalPassive);
    INFO("kinematicCloud: " << initialKinematic << " -> " << finalKinematic);

    REQUIRE(finalPassive == initialPassive);
    REQUIRE(finalKinematic == initialKinematic);

    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "Particle cell indices updated after load balancing",
    "[cloud][cellIndex][loadBalance][hex2D][hex3D][poly2D][poly3D][parallel]"
)
{
    if (!Pstream::parRun())
    {
        SUCCEED("Skipping cell index test in serial mode");
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

    IStringStream cloudPropsStream = kinematicCloudPropsStream("cellIndexCloud");
    IOdictionary cloudPropsDict
    (
        IOobject
        (
            "cellIndexCloudProperties",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloudPropsStream
    );
    cloudPropsDict.regIOobject::write();

    basicKinematicCloud cloud
    (
        "cellIndexCloud",
        rho,
        U,
        mu,
        g
    );

    // Add particles
    const auto& cc = mesh.cellCentres();
    label added = 0;

    forAll(cc, celli)
    {
        if (added < 40)
        {
            basicKinematicParcel* p = new basicKinematicParcel
            (
                mesh,
                cc[celli],
                celli
            );
            p->d() = 1e-4;
            p->nParticle() = 1.0;
            p->rho() = 1000.0;
            p->U() = Zero;

            cloud.addParticle(p);
            added++;
        }
    }

    label initialGlobal = returnReduce(cloud.nParcels(), sumOp<label>());
    REQUIRE(initialGlobal > 0);

    // Perform refinement + load balancing
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

    // After load balancing, verify each particle's cell index is valid
    // AND the particle position is actually inside that cell
    bool allCellsValid = true;
    bool allPositionsMatch = true;

    for (const basicKinematicParcel& p : cloud)
    {
        // Check cell index is in range
        if (p.cell() < 0 || p.cell() >= mesh.nCells())
        {
            allCellsValid = false;
            break;
        }

        // Check that particle position is inside its claimed cell
        // Use a tolerance-based point-in-cell check
        const point& pos = p.position();
        const label cellI = p.cell();

        // The cell should contain the particle position
        if (!mesh.pointInCell(pos, cellI))
        {
            allPositionsMatch = false;
            Pout << "Particle at " << pos << " claims cell " << cellI
                 << " but pointInCell returned false" << endl;
        }
    }

    REQUIRE(allCellsValid);
    REQUIRE(allPositionsMatch);

    // Verify total count unchanged
    label finalGlobal = returnReduce(cloud.nParcels(), sumOp<label>());
    REQUIRE(finalGlobal == initialGlobal);

    runTime.setTime(0.0, 0);
}


// ************************************************************************* //
