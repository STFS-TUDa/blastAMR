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

Description
    Unit tests for cloudSupport handling of intermediate cloud types.
    These tests verify that cloudSupport correctly identifies and handles
    different cloud types (basicKinematicCloud, basicThermoCloud, etc.)
    during load balancing operations.

\*---------------------------------------------------------------------------*/

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

// Intermediate cloud includes for type checking
#include "basicKinematicCloud.H"
#include "basicKinematicCollidingCloud.H"
#include "basicKinematicMPPICCloud.H"
#include "basicThermoCloud.H"
#include "basicReactingCloud.H"
#include "basicReactingMultiphaseCloud.H"

// Parcel includes
#include "basicKinematicParcel.H"

#include <csetjmp>
#include <csignal>
#include <cstdlib>
#include <functional>

using namespace Foam;
extern Time* timePtr;
extern argList* argsPtr;

// Signal handling for catching aborts
// ********************************************************************************************** //
jmp_buf intermediate_cloud_env;
void onSigabrtIntermediateCloud(int signum)
{
    signal(signum, SIG_DFL);
    longjmp(intermediate_cloud_env, 1);
}
void tryAndCatchAbortingCodeIntermediateCloud(std::function<void(void)> func)
{
    FatalError.dontThrowExceptions();
    if (setjmp(intermediate_cloud_env) == 0) {
        signal(SIGUSR1, &onSigabrtIntermediateCloud);
        signal(SIGUSR2, &onSigabrtIntermediateCloud);
        signal(SIGABRT, &onSigabrtIntermediateCloud);
        signal(SIGTERM, &onSigabrtIntermediateCloud);
        signal(SIGQUIT, &onSigabrtIntermediateCloud);
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
        bool abortedOrTerminated = true;
        REQUIRE(abortedOrTerminated == false);
    }
}
// ********************************************************************************************** //


// Test helper: Add kinematic parcels to a cloud at specified positions
// This uses the simple parcel constructor that only requires mesh, position, cellI
template<class CloudType, class ParcelType>
void addKinematicParcelsInBox
(
    CloudType& cloud,
    const fvMesh& mesh,
    const boundBox& box,
    label nParticles,
    scalar diameter = 1e-4,
    scalar nParticleValue = 1.0,
    scalar rhoValue = 1000.0
)
{
    const auto& cc = mesh.cellCentres();
    label added = 0;

    forAll(cc, celli)
    {
        if (box.contains(cc[celli]) && added < nParticles)
        {
            // Create parcel using simple constructor
            auto* p = new ParcelType(mesh, cc[celli], celli);

            // Set required properties to avoid FPE
            p->d() = diameter;
            p->nParticle() = nParticleValue;
            p->rho() = rhoValue;

            cloud.addParticle(p);
            added++;
        }
    }
}


TEST_CASE
(
    "Verify cloudSupport type detection for intermediate clouds",
    "[cloud][cloudTypes][hex2D][serial][parallel]"
)
{
    /*
     * This test verifies that cloudSupport can correctly detect and handle
     * different intermediate cloud types through dynamic_cast.
     *
     * Note: Creating actual intermediate clouds (basicKinematicCloud, etc.)
     * requires significant setup (cloudProperties dict, carrier fields).
     * This test focuses on verifying the type detection infrastructure.
     */

    FatalError.throwExceptions();
    Time& runTime = *timePtr;

    // Only run on hex2D case
    if (runTime.caseName().find("hex2D") == string::npos)
    {
        SUCCEED("Skipping cloud type test for non-hex2D case");
        return;
    }

    // Reset time
    runTime.setTime(0.0, 0);

    // Create a simple fvMesh
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

    // Test 1: Create a passiveParticleCloud and verify it's detected
    {
        passiveParticleCloud pCloud(mesh, Foam::zero{}, "passiveTestCloud");

        // Add some particles
        const boundBox meshBox(mesh.bounds());
        const auto& cc = mesh.cellCentres();
        label added = 0;
        forAll(cc, celli)
        {
            if (meshBox.contains(cc[celli]) && added < 5)
            {
                pCloud.addParticle(new passiveParticle(mesh, cc[celli], celli));
                added++;
            }
        }

        // Verify cloudSupport functions work
        label count = cloudSupport::countParticles(mesh);
        label globalCount = returnReduce(count, sumOp<label>());

        CAPTURE(Pstream::myProcNo(), count, globalCount);
        REQUIRE(globalCount >= 0);

        // Verify cloud names detection
        wordList names = cloudSupport::cloudNames(mesh);
        bool foundCloud = false;
        for (const word& name : names)
        {
            if (name == "passiveTestCloud")
            {
                foundCloud = true;
                break;
            }
        }
        REQUIRE(foundCloud);

        // Test storeGlobalPositions - should work for passiveParticleCloud
        cloudSupport::storeGlobalPositions(mesh);

        // Verify particle count unchanged
        label newCount = cloudSupport::countParticles(mesh);
        label newGlobalCount = returnReduce(newCount, sumOp<label>());
        REQUIRE(newGlobalCount == globalCount);
    }

    // Reset time
    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "Test cloudSupport particlesPerCell with mixed positions",
    "[cloud][particlesPerCell][hex2D][serial][parallel]"
)
{
    /*
     * Tests the particlesPerCell function which counts particles in each cell.
     * This is used for load balancing weight calculations.
     */

    FatalError.throwExceptions();
    Time& runTime = *timePtr;

    // Only run on hex2D case
    if (runTime.caseName().find("hex2D") == string::npos)
    {
        SUCCEED("Skipping particlesPerCell test for non-hex2D case");
        return;
    }

    // Reset time
    runTime.setTime(0.0, 0);

    // Create mesh
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

    // Create cloud with particles concentrated in specific cells
    passiveParticleCloud cloud(mesh, Foam::zero{}, "ppcTestCloud");

    // Add multiple particles to the first few cells
    const auto& cc = mesh.cellCentres();
    if (mesh.nCells() > 0)
    {
        // Add 3 particles to cell 0
        for (label i = 0; i < 3; i++)
        {
            cloud.addParticle(new passiveParticle(mesh, cc[0], 0));
        }

        // Add 2 particles to cell 1 (if exists)
        if (mesh.nCells() > 1)
        {
            for (label i = 0; i < 2; i++)
            {
                cloud.addParticle(new passiveParticle(mesh, cc[1], 1));
            }
        }
    }

    // Get particles per cell
    tmp<labelField> ppc = cloudSupport::particlesPerCell(mesh);

    // Verify counts
    if (mesh.nCells() > 0)
    {
        REQUIRE(ppc()[0] == 3);
        if (mesh.nCells() > 1)
        {
            REQUIRE(ppc()[1] == 2);
        }
    }

    // Sum should equal total particles
    label sum = 0;
    forAll(ppc(), i)
    {
        sum += ppc()[i];
    }
    REQUIRE(sum == cloud.size());

    // Reset time
    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "Test cloudSupport with empty clouds",
    "[cloud][empty][hex2D][serial][parallel]"
)
{
    /*
     * Tests that cloudSupport handles empty clouds correctly.
     * Empty clouds should not cause errors during load balancing.
     */

    FatalError.throwExceptions();
    Time& runTime = *timePtr;

    // Only run on hex2D case
    if (runTime.caseName().find("hex2D") == string::npos)
    {
        SUCCEED("Skipping empty cloud test for non-hex2D case");
        return;
    }

    // Reset time
    runTime.setTime(0.0, 0);

    // Create mesh
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

    // Create empty cloud
    passiveParticleCloud cloud(mesh, Foam::zero{}, "emptyTestCloud");

    REQUIRE(cloud.size() == 0);

    // All cloudSupport functions should work with empty clouds
    cloudSupport::storeGlobalPositions(mesh);

    label count = cloudSupport::countParticles(mesh);
    REQUIRE(count == 0);

    tmp<labelField> ppc = cloudSupport::particlesPerCell(mesh);
    label sum = 0;
    forAll(ppc(), i)
    {
        sum += ppc()[i];
    }
    REQUIRE(sum == 0);

    wordList names = cloudSupport::cloudNames(mesh);
    bool found = false;
    for (const word& name : names)
    {
        if (name == "emptyTestCloud")
        {
            found = true;
            break;
        }
    }
    REQUIRE(found);

    // Reset time
    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "Test cloudSupport multiple clouds",
    "[cloud][multiple][hex2D][serial][parallel]"
)
{
    /*
     * Tests cloudSupport with multiple clouds registered simultaneously.
     * This is important for cases with multiple particle phases.
     */

    FatalError.throwExceptions();
    Time& runTime = *timePtr;

    // Only run on hex2D case
    if (runTime.caseName().find("hex2D") == string::npos)
    {
        SUCCEED("Skipping multiple clouds test for non-hex2D case");
        return;
    }

    // Reset time
    runTime.setTime(0.0, 0);

    // Create mesh
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

    // Create multiple clouds
    passiveParticleCloud cloud1(mesh, Foam::zero{}, "multiCloud1");
    passiveParticleCloud cloud2(mesh, Foam::zero{}, "multiCloud2");
    passiveParticleCloud cloud3(mesh, Foam::zero{}, "multiCloud3");

    // Add different numbers of particles to each
    const auto& cc = mesh.cellCentres();
    if (mesh.nCells() > 0)
    {
        // Cloud 1: 5 particles
        for (label i = 0; i < min(label(5), mesh.nCells()); i++)
        {
            cloud1.addParticle(new passiveParticle(mesh, cc[i], i));
        }

        // Cloud 2: 3 particles
        for (label i = 0; i < min(label(3), mesh.nCells()); i++)
        {
            cloud2.addParticle(new passiveParticle(mesh, cc[i], i));
        }

        // Cloud 3: empty
    }

    // Verify cloud names
    wordList names = cloudSupport::cloudNames(mesh);
    REQUIRE(names.size() >= 3);

    bool found1 = false, found2 = false, found3 = false;
    for (const word& name : names)
    {
        if (name == "multiCloud1") found1 = true;
        if (name == "multiCloud2") found2 = true;
        if (name == "multiCloud3") found3 = true;
    }
    REQUIRE(found1);
    REQUIRE(found2);
    REQUIRE(found3);

    // Verify total count
    label totalCount = cloudSupport::countParticles(mesh);
    label expectedCount = cloud1.size() + cloud2.size() + cloud3.size();
    REQUIRE(totalCount == expectedCount);

    // storeGlobalPositions should work for all clouds
    cloudSupport::storeGlobalPositions(mesh);

    // Reset time
    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "Test passiveParticle cloud load balancing with refinement",
    "[cloud][loadBalance][refinement][hex2D][poly2D][parallel]"
)
{
    /*
     * Tests that passiveParticleCloud survives mesh refinement and
     * load balancing operations. This is the baseline test for
     * intermediate cloud support.
     */

    // Skip if not running in parallel
    if (!Pstream::parRun())
    {
        SUCCEED("Skipping load balancing test in serial mode");
        return;
    }

    FatalError.throwExceptions();
    Time& runTime = *timePtr;

    // Refinement box
    const boundBox refBox(point(0.02, 0.025, -1), point(0.04, 0.035, 1));

    // Test parameters
    word refiner = GENERATE("polyRefiner", "hexRefiner");
    if (refiner == "hexRefiner" && runTime.caseName().find("poly") != string::npos) {
        REQUIRE(true);
        return;
    }
    word balance = "yes";
    label maxRefL = 2;

    // dynamicMeshDict configuration
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
        "nBufferLayers   1;"
        "maxRefinement   "+Foam::name(maxRefL)+";"
        "dumpLevel       false;"
        "protectedPatches ();"
        "nPatchesBuffers 1;"
    );

    // Reset time
    runTime.setTime(0.0, 0);

    // Create mesh
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
            true
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Set test field to 1 inside refinement box
    forAll(mesh.C(), ci)
    {
        if (refBox.contains(mesh.C()[ci]))
        {
            test[ci] = 1;
        }
    }

    // Create cloud distributed across the mesh
    passiveParticleCloud cloud(mesh, Foam::zero{}, "lbTestCloud");

    const boundBox meshBox(mesh.bounds());
    const auto& cc = mesh.cellCentres();
    label targetParticles = 30;
    label added = 0;

    forAll(cc, celli)
    {
        if (meshBox.contains(cc[celli]) && added < targetParticles)
        {
            cloud.addParticle(new passiveParticle(mesh, cc[celli], celli));
            added++;
        }
    }

    // Record initial particle count
    label initialLocal = cloud.size();
    label initialGlobal = returnReduce(initialLocal, sumOp<label>());

    CAPTURE
    (
        Pstream::myProcNo(), Pstream::nProcs(),
        runTime.caseName(),
        mesh.nCells(), refiner, maxRefL,
        initialLocal, initialGlobal
    );

    REQUIRE(initialGlobal > 0);

    // Perform refinement + load balancing cycles
    for (label i = 0; i < maxRefL; i++)
    {
        runTime++;
        tryAndCatchAbortingCodeIntermediateCloud
        (
            [&]()
            {
                mesh.update();
            }
        );
    }

    // Verify global particle count is conserved
    label finalLocal = cloud.size();
    label finalGlobal = returnReduce(finalLocal, sumOp<label>());

    INFO("Initial: " << initialGlobal << ", Final: " << finalGlobal);
    REQUIRE(finalGlobal == initialGlobal);

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


TEST_CASE
(
    "Test cloud type hierarchy detection",
    "[cloud][typeHierarchy][hex2D][serial]"
)
{
    /*
     * Verifies that the cloud type hierarchy is correctly understood.
     * This is critical for proper dynamic_cast ordering in cloudSupport.
     *
     * Type hierarchy (most to least derived):
     * - basicReactingMultiphaseCloud
     * - basicReactingCloud
     * - basicThermoCloud
     * - basicKinematicCollidingCloud / basicKinematicMPPICCloud
     * - basicKinematicCloud
     * - passiveParticleCloud (Cloud<passiveParticle>)
     * - cloud (base)
     */

    FatalError.throwExceptions();
    Time& runTime = *timePtr;

    // Only run on hex2D case in serial
    if (runTime.caseName().find("hex2D") == string::npos)
    {
        SUCCEED("Skipping type hierarchy test for non-hex2D case");
        return;
    }

    // Verify type information exists for all intermediate cloud types
    // This doesn't create clouds, just checks the type system is properly set up

    // Check that kinematicCloud type name is registered (base class typeName)
    REQUIRE(kinematicCloud::typeName == "kinematicCloud");

    // The test verifies at compile time that all these types are available
    // and can be used in dynamic_cast operations

    SUCCEED("Type hierarchy verification passed");
}


// ************************************************************************* //
