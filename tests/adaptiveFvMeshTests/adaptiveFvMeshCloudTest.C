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

// Lagrangian includes
#include "passiveParticleCloud.H"
#include "cloudSupport.H"

#include <csetjmp>
#include <csignal>
#include <cstdlib>
#include <functional>

using namespace Foam;
extern Time* timePtr;
extern argList* argsPtr;

// This is to make sure OpenFOAM doesn't bail on us at any point
// Typically, some functions with abort even if you set throwExceptions; which is not cool
// also, code might run into deadlocks or infinite loops, timeout will kill try to kill the process
// This section protects against these scenatios
// ********************************************************************************************** //
jmp_buf cloud_amr_env;
void onSigabrtCloud(int signum)
{
  signal (signum, SIG_DFL);
  longjmp (cloud_amr_env, 1);
}
void tryAndCatchAbortingCodeCloud(std::function<void(void)> func)
{
    FatalError.dontThrowExceptions();
    if (setjmp (cloud_amr_env) == 0) {
        signal(SIGUSR1, &onSigabrtCloud);
        signal(SIGUSR2, &onSigabrtCloud);
        signal(SIGABRT, &onSigabrtCloud);
        signal(SIGTERM, &onSigabrtCloud);
        signal(SIGQUIT, &onSigabrtCloud);
        func();
        signal(SIGQUIT, SIG_DFL);
        signal(SIGUSR1, SIG_DFL);
        signal(SIGUSR2, SIG_DFL);
        signal(SIGABRT, SIG_DFL);
        signal(SIGTERM, SIG_DFL);
    }
    else {
        Pout<< "Either this code tried to abort or there was"
            " an attempt to terminate it (e.g. with a timeout) on " <<
            Pstream::myProcNo() << "..." << endl;
        // Need to fail the test case now
        bool abortedOrTerminated = true;
        REQUIRE(abortedOrTerminated == false);
    }
}
// ********************************************************************************************** //

// Helper function to add particles to a cloud at specified positions
void addParticlesInBox
(
    passiveParticleCloud& cloud,
    const fvMesh& mesh,
    const boundBox& box,
    label nParticles
)
{
    const auto& cc = mesh.cellCentres();
    label added = 0;

    forAll(cc, celli)
    {
        if (box.contains(cc[celli]) && added < nParticles)
        {
            // Add particle at cell centre
            cloud.addParticle(new passiveParticle(mesh, cc[celli], celli));
            added++;
        }
    }
}


TEST_CASE
(
    "Check cloud handling during mesh refinement",
    "[cloud][hex2D][hex3D][poly2D][poly3D][serial][parallel]"
)
{
    FatalError.throwExceptions();
    Time& runTime = *timePtr;
    argList& args = *argsPtr;

    // Refinement box - same as in adaptiveFvMeshTest.C
    const word boxString = "(0.02 0.025 -1) (0.04 0.035 1)";
    const boundBox refBox(point(0.02, 0.025, -1), point(0.04, 0.035, 1));

    // Tested variables' matrix
    word refiner = GENERATE("polyRefiner", "hexRefiner");
    if (refiner == "hexRefiner" && runTime.caseName().find("poly") != string::npos) {
        REQUIRE(true);
        return;
    }
    word balance = Pstream::parRun() ? GENERATE("no") : GENERATE("no");
    label nBufferLayers = GENERATE(1);
    label maxRefL = GENERATE(2);

    // Supported constant/dynamicMeshDict entries
    IStringStream is
    (
        "dynamicFvMesh   adaptiveFvMesh;"
        "errorEstimator  fieldValue;"
        "fieldName       test;"
        "balance         "+balance+";"
        "refiner         "+refiner+";"
        "refineInterval  1;"
        "unrefineInterval  1;"
        "balanceInterval  1;"
        "lowerRefineLevel 0.01;"
        "unrefineLevel 0.01;"
        "nBufferLayers   "+Foam::name(nBufferLayers)+";"
        "maxRefinement   "+Foam::name(maxRefL)+";"
        "dumpLevel       false;"
        "protectedPatches ();"
        "nPatchesBuffers 1;"
    );

    // Reset time
    runTime.setTime(0.0, 0);

    // Create mesh object
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
        is
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

    // Create AMR criterion field
    volScalarField test
    (
        IOobject
        (
            "test",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Set test field to 1 inside refinement box
    forAll(mesh.C(), ci)
    {
        const auto& center = mesh.C()[ci];
        if (refBox.contains(center))
        {
            test[ci] = 1;
        }
    }

    // Create a passive particle cloud
    passiveParticleCloud cloud(mesh, Foam::zero{}, "testCloud");

    // Add particles inside the refinement box
    const label targetParticles = 10;
    addParticlesInBox(cloud, mesh, refBox, targetParticles);

    // Get initial particle count (global)
    label initialLocalCount = cloud.size();
    label initialGlobalCount = returnReduce(initialLocalCount, sumOp<label>());

    // Capture test variables for sane reports
    CAPTURE
    (
        Pstream::myProcNo(), Pstream::nProcs(),
        runTime.caseName(),
        mesh.nCells(), balance,
        refiner, maxRefL, nBufferLayers,
        initialGlobalCount
    );

    REQUIRE(initialGlobalCount > 0);

    // Perform refinement cycles
    for (label i = 0; i < maxRefL; i++)
    {
        runTime++;
        tryAndCatchAbortingCodeCloud
        (
            [&]()
            {
                mesh.update();
            }
        );
    }

    // Check particle count after refinement
    label finalLocalCount = cloud.size();
    label finalGlobalCount = returnReduce(finalLocalCount, sumOp<label>());

    // Particles should be conserved during refinement
    // Note: Some particles may be lost if they end up outside the mesh
    // but for internal particles this should not happen
    INFO("Initial particles: " << initialGlobalCount << ", Final particles: " << finalGlobalCount);
    REQUIRE(finalGlobalCount == initialGlobalCount);

    // Verify particles are still in valid cells
    bool allValid = true;
    for (const passiveParticle& p : cloud)
    {
        if (p.cell() < 0 || p.cell() >= mesh.nCells())
        {
            allValid = false;
            break;
        }
    }
    REQUIRE(allValid);

    // Reset time for good measure
    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "Check cloud handling during load balancing",
    "[cloud][cloudLB][hex2D][hex3D][poly2D][poly3D][parallel]"
)
{
    // Skip if not running in parallel
    if (!Pstream::parRun())
    {
        SUCCEED("Skipping load balancing test in serial mode");
        return;
    }

    FatalError.throwExceptions();
    Time& runTime = *timePtr;
    argList& args = *argsPtr;

    // Refinement box
    const word boxString = "(0.02 0.025 -1) (0.04 0.035 1)";
    const boundBox refBox(point(0.02, 0.025, -1), point(0.04, 0.035, 1));

    // Tested variables' matrix - with load balancing enabled
    word refiner = GENERATE("polyRefiner", "hexRefiner");
    if (refiner == "hexRefiner" && runTime.caseName().find("poly") != string::npos) {
        REQUIRE(true);
        return;
    }
    word balance = "yes";
    label nBufferLayers = 1;
    label maxRefL = 2;

    // Supported constant/dynamicMeshDict entries
    IStringStream is
    (
        "dynamicFvMesh   adaptiveFvMesh;"
        "errorEstimator  fieldValue;"
        "fieldName       test;"
        "balance         "+balance+";"
        "refiner         "+refiner+";"
        "refineInterval  1;"
        "unrefineInterval  1;"
        "balanceInterval  1;"
        "lowerRefineLevel 0.01;"
        "unrefineLevel 0.01;"
        "nBufferLayers   "+Foam::name(nBufferLayers)+";"
        "maxRefinement   "+Foam::name(maxRefL)+";"
        "dumpLevel       false;"
        "protectedPatches ();"
        "nPatchesBuffers 1;"
    );

    // Reset time
    runTime.setTime(0.0, 0);

    // Create mesh object
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
        is
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

    // Create AMR criterion field - set to 1 everywhere to trigger significant refinement
    volScalarField test
    (
        IOobject
        (
            "test",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Set test field to 1 inside refinement box
    forAll(mesh.C(), ci)
    {
        const auto& center = mesh.C()[ci];
        if (refBox.contains(center))
        {
            test[ci] = 1;
        }
    }

    // Create a passive particle cloud
    passiveParticleCloud cloud(mesh, Foam::zero{}, "testCloud");

    // Add particles distributed across the mesh (not just in refinement box)
    // to better test load balancing
    const boundBox meshBox(mesh.bounds());
    const label targetParticles = 20;
    addParticlesInBox(cloud, mesh, meshBox, targetParticles);

    // Get initial particle count (global)
    label initialLocalCount = cloud.size();
    label initialGlobalCount = returnReduce(initialLocalCount, sumOp<label>());

    // Capture test variables
    CAPTURE
    (
        Pstream::myProcNo(), Pstream::nProcs(),
        runTime.caseName(),
        mesh.nCells(), balance,
        refiner, maxRefL,
        initialGlobalCount, initialLocalCount
    );

    REQUIRE(initialGlobalCount > 0);

    // Perform refinement cycles - this should trigger load balancing
    for (label i = 0; i < maxRefL; i++)
    {
        runTime++;
        tryAndCatchAbortingCodeCloud
        (
            [&]()
            {
                mesh.update();
            }
        );
    }

    // Check particle count after refinement + load balancing
    label finalLocalCount = cloud.size();
    label finalGlobalCount = returnReduce(finalLocalCount, sumOp<label>());

    // Global particle count must be conserved
    INFO("Initial global: " << initialGlobalCount << ", Final global: " << finalGlobalCount);
    REQUIRE(finalGlobalCount == initialGlobalCount);

    // Local counts may differ due to redistribution - just verify they sum correctly
    labelList localCounts(Pstream::nProcs(), 0);
    localCounts[Pstream::myProcNo()] = finalLocalCount;
    reduce(localCounts, sumOp<labelList>());

    label sumLocal = 0;
    forAll(localCounts, proci)
    {
        sumLocal += localCounts[proci];
    }
    REQUIRE(sumLocal == finalGlobalCount);

    // Verify all particles are in valid cells on their respective processors
    bool allValid = true;
    for (const passiveParticle& p : cloud)
    {
        if (p.cell() < 0 || p.cell() >= mesh.nCells())
        {
            allValid = false;
            break;
        }
    }
    REQUIRE(allValid);

    // Reset time
    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "Check cloud utility functions",
    "[cloud][hex2D][serial][parallel]"
)
{
    FatalError.throwExceptions();
    Time& runTime = *timePtr;

    // Only run on hex2D case for simplicity
    if (runTime.caseName().find("hex2D") == string::npos)
    {
        SUCCEED("Skipping utility test for non-hex2D case");
        return;
    }

    // Reset time
    runTime.setTime(0.0, 0);

    // Create a simple fvMesh (not adaptiveFvMesh)
    fvMesh mesh
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

    // Create a passive particle cloud
    passiveParticleCloud cloud(mesh, Foam::zero{}, "utilTestCloud");

    // Add some particles
    const boundBox meshBox(mesh.bounds());
    addParticlesInBox(cloud, mesh, meshBox, 15);

    label cloudSize = cloud.size();
    label globalSize = returnReduce(cloudSize, sumOp<label>());

    CAPTURE(Pstream::myProcNo(), cloudSize, globalSize);

    // Test cloudSupport::countParticles
    label countedParticles = cloudSupport::countParticles(mesh);
    label globalCounted = returnReduce(countedParticles, sumOp<label>());
    REQUIRE(globalCounted == globalSize);

    // Test cloudSupport::cloudNames
    wordList names = cloudSupport::cloudNames(mesh);
    bool foundCloud = false;
    forAll(names, i)
    {
        if (names[i] == "utilTestCloud")
        {
            foundCloud = true;
            break;
        }
    }
    REQUIRE(foundCloud);

    // Test cloudSupport::particlesPerCell
    tmp<labelField> ppc = cloudSupport::particlesPerCell(mesh);
    label sumPPC = 0;
    forAll(ppc(), i)
    {
        sumPPC += ppc()[i];
    }
    REQUIRE(sumPPC == cloudSize);

    // Reset time
    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "Check load balancing with oldTime field mismatch",
    "[cloud][cloudLB][oldTime][poly2D][poly3D][parallel]"
)
{
    // Skip if not running in parallel
    if (!Pstream::parRun())
    {
        SUCCEED("Skipping oldTime mismatch test in serial mode");
        return;
    }

    FatalError.throwExceptions();
    Time& runTime = *timePtr;

    // Refinement box
    const boundBox refBox(point(0.02, 0.025, -1), point(0.04, 0.035, 1));

    // Only use polyRefiner (hexRefiner has separate bugs)
    word refiner = "polyRefiner";
    word balance = "yes";
    label nBufferLayers = 1;
    label maxRefL = 2;

    // Supported constant/dynamicMeshDict entries
    IStringStream is
    (
        "dynamicFvMesh   adaptiveFvMesh;"
        "errorEstimator  fieldValue;"
        "fieldName       test;"
        "balance         "+balance+";"
        "refiner         "+refiner+";"
        "refineInterval  1;"
        "unrefineInterval  1;"
        "balanceInterval  1;"
        "lowerRefineLevel 0.01;"
        "unrefineLevel 0.01;"
        "nBufferLayers   "+Foam::name(nBufferLayers)+";"
        "maxRefinement   "+Foam::name(maxRefL)+";"
        "dumpLevel       false;"
        "protectedPatches ();"
        "nPatchesBuffers 1;"
    );

    // Reset time
    runTime.setTime(0.0, 0);

    // Create mesh object
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
        is
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

    // Create AMR criterion field
    volScalarField test
    (
        IOobject
        (
            "test",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Set test field to 1 inside refinement box
    forAll(mesh.C(), ci)
    {
        const auto& center = mesh.C()[ci];
        if (refBox.contains(center))
        {
            test[ci] = 1;
        }
    }

    // Create a velocity field U - this mimics what a real solver would have
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, Zero)
    );

    // KEY TEST: Create oldTime only on some processors (simulating kinematicCloud behavior)
    // In kinematicCloud, U.oldTime() is called when particles hit walls with fixesValue BCs.
    // This happens inconsistently across processors based on particle locations.
    if (Pstream::myProcNo() == 0 || Pstream::myProcNo() == 1)
    {
        // Trigger oldTime creation on processors 0 and 1 only
        U.oldTime();
    }

    // Create a passive particle cloud
    passiveParticleCloud cloud(mesh, Foam::zero{}, "testCloud");

    // Add particles distributed across the mesh
    const boundBox meshBox(mesh.bounds());
    const label targetParticles = 20;
    addParticlesInBox(cloud, mesh, meshBox, targetParticles);

    // Get initial particle count (global)
    label initialGlobalCount = returnReduce(cloud.size(), sumOp<label>());

    CAPTURE
    (
        Pstream::myProcNo(), Pstream::nProcs(),
        runTime.caseName(),
        mesh.nCells(), balance,
        initialGlobalCount
    );

    REQUIRE(initialGlobalCount > 0);

    // Perform refinement cycles - this should trigger load balancing
    // The syncOldTimeFields function should handle the U_0 mismatch
    for (label i = 0; i < maxRefL; i++)
    {
        runTime++;
        tryAndCatchAbortingCodeCloud
        (
            [&]()
            {
                mesh.update();
            }
        );
    }

    // Check particle count after refinement + load balancing
    label finalGlobalCount = returnReduce(cloud.size(), sumOp<label>());

    // Global particle count must be conserved
    INFO("Initial global: " << initialGlobalCount << ", Final global: " << finalGlobalCount);
    REQUIRE(finalGlobalCount == initialGlobalCount);

    // Verify all particles are in valid cells
    bool allValid = true;
    for (const passiveParticle& p : cloud)
    {
        if (p.cell() < 0 || p.cell() >= mesh.nCells())
        {
            allValid = false;
            break;
        }
    }
    REQUIRE(allValid);

    // Reset time
    runTime.setTime(0.0, 0);
}
