#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.get<word>("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
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

volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"
#include "createFvOptions.H"

/* Additional fields specific to the solver */
// Mean species molecular weight
volScalarField WMean
(
    IOobject
    (
        "WMean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimMoles, Zero)
);

// Species molecular weights
PtrList<dimensionedScalar> W_k(Y.size());
forAll(Y, i){
    W_k.set
    (
        i,
        new dimensionedScalar
        (
            "W_"+Y[i].name(),
            dimMass/dimMoles,
            composition.W(i)
        )
    );
}

// Species specific sensible enthalpies
PtrList<volScalarField> hs(Y.size());
forAll(hs, i)
{
    hs.set
    (
        i,
        new volScalarField 
    	(
            IOobject
            (
                "hs_"+Y[i].name(),
        	    runTime.timeName(),
	            mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimEnergy/dimMass,0.0),
            "zeroGradient"
        )
    );
}

volScalarField hSpeciesDiffusion
(
    IOobject
    (
        "hSpeciesDiffusion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("hSpeciesDiffusion", dimEnergy/dimTime/dimVolume, 0.0)
);

volScalarField hSpeciesCorrection
(
    IOobject
    (
        "hSpeciesCorrection",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("hSpeciesCorrection", dimEnergy/dimTime/dimVolume, 0.0)
);

// Correction velocity momentum flux
surfaceScalarField phi_corr
(
    IOobject
    (
        "phi_corr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("phi_corr", dimMass/dimTime, 0.0)
);

// Diffusive fluxes
PtrList<surfaceScalarField> J_(Y.size());
forAll(Y, i)
{
    J_.set (
	i,
	new surfaceScalarField
	(
	    IOobject
	    (
	 	"J_" + Y[i].name(),
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("J", dimensionSet(1, -0, -1, 0, 0), 0.0)
	)
    );
}

// Dictionary which stores Lewis and Prandtl number information
IOdictionary constLewisProperties
(
  IOobject
   (
    "constLewisProperties",
    runTime.constant(),
    runTime,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
   )
);

// Prandtl number of mixture
const scalar Pr_mix(readScalar(constLewisProperties.lookup("Prandtl")));

// Species Lewis numbers
scalarList Le_k(Y.size());
forAll(Y, i){
    Le_k[i] = readScalar(constLewisProperties.subDict(Y[i].name()).lookup("Lewis"));
}

