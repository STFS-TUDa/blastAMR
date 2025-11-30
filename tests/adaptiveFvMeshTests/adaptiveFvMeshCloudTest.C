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
#include "passiveParticleCloud.H"
#include "cloudSupport.H"
#include "abortHandle.H"

using namespace Foam;
extern Time* timePtr;
extern argList* argsPtr;

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
            cloud.addParticle(new passiveParticle(mesh, cc[celli], celli));
            added++;
        }
    }
}


TEST_CASE
(
    "Check cloud handling during mesh refinement/lb",
    "[cloud][hex2D][hex3D][poly2D][poly3D][serial][parallel]"
)
{
    FatalError.throwExceptions();
    Time& runTime = *timePtr;
    argList& args = *argsPtr;

    // Refinement box
    const word boxString = "(0.02 0.025 -1) (0.04 0.035 1)";
    const boundBox refBox(point(0.02, 0.025, -1), point(0.04, 0.035, 1));

    // Tested variables' matrix
    word refiner = GENERATE("polyRefiner", "hexRefiner");
    // Skip hexRefiner for polyhedral meshes (hexRefiner only works on hex meshes)
    if (refiner == "hexRefiner" &&
        runTime.caseName().find("poly") != string::npos) {
        SUCCEED("Skipping hexRefiner for polyhedral mesh");
        return;
    }
    word balance = Pstream::parRun() ? GENERATE("yes", "no") : GENERATE("no");
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
        tryAndCatchAbortingCode
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
    "Check load balancing with oldTime field mismatch",
    "[cloud][cloudLB][oldTime][hex2D][hex3D][poly2D][poly3D][parallel]"
)
{
    FatalError.throwExceptions();
    Time& runTime = *timePtr;

    // Refinement box
    const boundBox refBox(point(0.02, 0.025, -1), point(0.04, 0.035, 1));

    // Use polyRefiner for poly cases, hexRefiner for hex cases
    word refiner = runTime.caseName().find("poly") != string::npos
        ? "polyRefiner" : "hexRefiner";
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
        tryAndCatchAbortingCode
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
