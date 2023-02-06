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
jmp_buf amr_env;
void onSigabrt(int signum)
{
  signal (signum, SIG_DFL);
  longjmp (amr_env, 1);
}
void tryAndCatchAbortingCode(std::function<void(void)> func)
{
    FatalError.throwExceptions();
    if (setjmp (amr_env) == 0) {
        signal(SIGUSR1, &onSigabrt);
        signal(SIGUSR2, &onSigabrt);
        signal(SIGABRT, &onSigabrt);
        signal(SIGTERM, &onSigabrt);
        signal(SIGQUIT, &onSigabrt);
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

TEST_CASE("Check polyhedral AMR/LB functionality for adaptiveFvMesh", "[hex2D][hex3D][poly2D][poly3D][serial][parallel][!throws]")
{
    FatalError.throwExceptions();
    Time& runTime = *timePtr;
    argList& args = *argsPtr;

    // Refinement box
    const word boxString = "(0.02 0.025 -1) (0.04 0.035 1)";
    const std::vector<scalar> boxBounds = {0.02, 0.04, 0.025, 0.035};

    // Tested variables' matrix 
    word refiner = GENERATE("polyRefiner");
    word balance = Pstream::parRun() ? GENERATE("no", "yes") : GENERATE("no");
    label maxRefL = GENERATE(1, 3);

    // Supported constant/dynamicMeshDict entries
    IStringStream is
    (
        "dynamicFvMesh   adaptiveFvMesh;"
        "errorEstimator  delta;"
        "deltaField      test;"
        "balance         "+balance+";"
        "refiner         "+refiner+";"
        "refineInterval  1;"
        "unrefineInterval  1;"
        "balanceInterval  1;"
        "lowerRefineLevel 0.01;"
        "unrefineLevel 0.01;"
        "nBufferLayers   3;"
        "maxRefinement   "+Foam::name(maxRefL)+";"
        "dumpLevel       false;"
    );

    SECTION("Box is refined incrementally " + Foam::name(maxRefL) + " times")
    {
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

        // Capture test variables for sane reports
        CAPTURE
        (
            Pstream::myProcNo(), Pstream::nProcs(),
            mesh.nCells(), balance,
            refiner, maxRefL
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
            dimensionedScalar("zero",dimless,0.0)
        );

        // Set test to 1 inside the refined box
        forAll(mesh.C(), ci)
        {
            const auto& center = mesh.C()[ci];
            if 
            (
                center.x() >= boxBounds[0] && center.x() <= boxBounds[1]
                && center.y() >= boxBounds[2] && center.y() <= boxBounds[3]
            )
            {
                test[ci] = 1;
            } else {
                test[ci] = 0;
            }
        }

        // Cell set for the refined box
        IStringStream refBoxIs(boxString);
        boxToCell refBox(mesh, refBoxIs);
        cellSet refBoxCells(mesh, "refBoxCells", IOobject::NO_READ);
        refBox.applyToSet(topoSetSource::setAction::NEW, refBoxCells);

        // Record initial number of cells
        label origRefBoxNCells = returnReduce(refBoxCells.size(), sumOp<label>());
        label newRefBoxNCells = -1;

        // run a refinement cycle, minor unrefinement should happen here, but this not checked.
        for(label i = 0; i < maxRefL; i++)
        {
            runTime++;
            // Update mesh
            tryAndCatchAbortingCode
            (
                [&]()
                {
                  mesh.update();
                }
            );
        }

        // Check the refined box again with a cellSet
        IStringStream newRefBoxIs(boxString);
        boxToCell newRefBox(mesh, newRefBoxIs);
        cellSet newRefBoxCells(mesh, "newRefBoxCells", IOobject::NO_READ);
        newRefBox.applyToSet(topoSetSource::setAction::NEW, newRefBoxCells);

        newRefBoxNCells = returnReduce(newRefBoxCells.size(), sumOp<label>());

        // Require that the new box gets refined to at least 4xmaxRefLx the original cell count
        // This is an inaccurate check; but there is nothing we can do.
        CHECK(newRefBoxNCells >= 4 * origRefBoxNCells);
    }

    SECTION("Box is unrefined, by moving the refinement region to somewhere else.")
    {
        // Reset time
        runTime.setTime(0.0, 0);

        // Create the mesh
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

        // Capture test variables for sane reports
        CAPTURE
        (
            Pstream::myProcNo(), Pstream::nProcs(),
            mesh.nCells(), balance,
            refiner, maxRefL
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
            dimensionedScalar("zero",dimless,0.0)
        );

        // Set test to 1 inside the refined box
        forAll(mesh.C(), ci)
        {
            const auto& center = mesh.C()[ci];
            if 
            (
                center.x() >= boxBounds[0] && center.x() <= boxBounds[1]
                && center.y() >= boxBounds[2] && center.y() <= boxBounds[3]
            )
            {
                test[ci] = 1;
            } else {
                test[ci] = 0;
            }
        }

        for(label i = 0; i < maxRefL; i++)
        {
            runTime++;
            // Update mesh
            tryAndCatchAbortingCode
            (
                [&]()
                {
                    mesh.update();
                }
            );
        }

        IStringStream boxIs(boxString);
        boxToCell box(mesh, boxIs);
        cellSet boxCells(mesh, "boxCells", IOobject::NO_READ);
        box.applyToSet(topoSetSource::setAction::NEW, boxCells);

        label origBoxNCells = returnReduce(boxCells.size(), sumOp<label>());
        label newBoxNCells = -1;

        // nuke test field, this must trigger unrefinement on next mesh update
        test = 0.0;
        runTime++;
        tryAndCatchAbortingCode
        (
            [&]()
            {
                mesh.update();
            }
        );

        // Check the refinement box for cell count again
        IStringStream newBoxIs(boxString);
        boxToCell newBox(mesh, newBoxIs);
        cellSet newBoxCells(mesh, "newBoxCells", IOobject::NO_READ);
        newBox.applyToSet(topoSetSource::setAction::NEW, newBoxCells);

        newBoxNCells = returnReduce(newBoxCells.size(), sumOp<label>());

        // Require that the old box gets unrefined
        CHECK(newBoxNCells < origBoxNCells);
    }

    // Reset time for good mesure
    runTime.setTime(0.0, 0);
}
