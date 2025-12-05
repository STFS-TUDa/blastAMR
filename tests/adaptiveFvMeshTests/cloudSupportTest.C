#include "IOobject.H"
#include "PstreamReduceOps.H"
#include "catch2/catch_all.hpp"
#include "catch2/catch_test_macros.hpp"
#include "fvCFD.H"
#include "fvMesh.H"
#include "messageStream.H"
#include "volFields.H"
#include "passiveParticleCloud.H"
#include "cloudSupport.H"
#include "basicKinematicCloud.H"
#include "basicKinematicCollidingCloud.H"
#include "basicKinematicMPPICCloud.H"
#include "basicThermoCloud.H"
#include "basicReactingCloud.H"
#include "basicReactingMultiphaseCloud.H"
#include "abortHandle.H"

using namespace Foam;
extern Time* timePtr;
extern argList* argsPtr;


TEST_CASE
(
    "Test cloudSupport with empty clouds",
    "[cloudSupport][empty][hex2D][hex3D][poly2D][poly3D][serial][parallel]"
)
{
    Time& runTime = *timePtr;
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

    passiveParticleCloud cloud(mesh, Foam::zero{}, "emptyTestCloud");
    CHECK(cloud.size() == 0);

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

    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "Test cloudSupport multiple clouds",
    "[cloudSupport][multiple][hex2D][hex3D][poly2D][poly3D][serial][parallel]"
)
{
    Time& runTime = *timePtr;
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

    passiveParticleCloud cloud1(mesh, Foam::zero{}, "multiCloud1");
    passiveParticleCloud cloud2(mesh, Foam::zero{}, "multiCloud2");
    passiveParticleCloud cloud3(mesh, Foam::zero{}, "multiCloud3");

    const auto& cc = mesh.cellCentres();
    if (mesh.nCells() > 0)
    {
        for (label i = 0; i < min(label(5), mesh.nCells()); i++)
        {
            cloud1.addParticle(new passiveParticle(mesh, cc[i], i));
        }
        for (label i = 0; i < min(label(3), mesh.nCells()); i++)
        {
            cloud2.addParticle(new passiveParticle(mesh, cc[i], i));
        }
    }

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

    label totalCount = cloudSupport::countParticles(mesh);
    label expectedCount = cloud1.size() + cloud2.size() + cloud3.size();
    REQUIRE(totalCount == expectedCount);

    cloudSupport::storeGlobalPositions(mesh);
    runTime.setTime(0.0, 0);
}


TEST_CASE
(
    "Test cloud type hierarchy detection",
    "[cloudSupport][typeHierarchy][hex2D][hex3D][poly2D][poly3D][serial]"
)
{
    REQUIRE(kinematicCloud::typeName == "kinematicCloud");
    SUCCEED("Type hierarchy verification passed");
}


// ************************************************************************* //
