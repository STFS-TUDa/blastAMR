Info<< nl << "Reading thermophysicalProperties" << endl;
autoPtr<hCombustionThermo> thermo
(
    hCombustionThermo::New(mesh)
);

combustionMixture& composition = thermo->composition();
PtrList<volScalarField>& Y = composition.Y();

word inertSpecie(thermo->lookup("inertSpecie"));

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo->rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField& p = thermo->p();
const volScalarField& psi = thermo->psi();
const volScalarField& T = thermo->T();
volScalarField& h = thermo->h();


#include "compressibleCreatePhi.H"

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo()
    )
);

Info<< "Creating field DpDt\n" << endl;
volScalarField DpDt =
    fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

Info << "Constructing chemical mechanism" << endl;
chemistryModel chemistry
(
    thermo(),
    rho
);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

for(label i=0; i<Y.size(); i++)
{
    fields.add(Y[i]);
}
fields.add(h);

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
);
