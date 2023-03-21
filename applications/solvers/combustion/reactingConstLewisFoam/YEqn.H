tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);


{
    /* Additional calculations specific to this solver */
    // Reset to zero
    WMean *= 0.0;
    hSpeciesCorrection *= 0.0;
    hSpeciesDiffusion *= 0.0;
    phi_corr *= 0.0;

    forAll(Y, i){
        J_[i] *= 0.0;
    }

    // Calculate species specific enthalpies
    forAll(Y, i){
        WMean += Y[i] * W_k[i];
    } 

    // Calculate species specific sensible enthalpies
    forAll(Y, i){
        forAll(mesh.cells(), cell){
            hs[i].primitiveFieldRef()[cell] = composition.Hs(i, p[cell], T[cell]);
        }
        hs[i].correctBoundaryConditions();
    }

    // Calculate gradient of mean molecular weight of mixture
    surfaceScalarField gradWMean(fvc::snGrad(dimensionedScalar("",dimless,1.0)*WMean)*mesh.magSf());

    // Calculate laminar diffusive fluxes for correction velocity
    forAll(Y, i){
        volScalarField& Yi = Y[i];
        J_[i] = -fvc::interpolate(thermo.mu()/(Le_k[i] * Pr_mix)) * (fvc::snGrad(Yi) * mesh.magSf()) - fvc::interpolate((thermo.mu()/(Le_k[i] * Pr_mix)) * Yi / WMean) * gradWMean;
    }

    // Calculate correction velocity
    forAll(Y,i){
        phi_corr -= J_[i];
    }
    forAll(Y,i){
        J_[i] += phi_corr*fvc::interpolate(Y[i]);
    }

    // Calculate enthalpy term due to species diffusion
    forAll(Y,i)
    {
        volScalarField& Yi = Y[i];
        volScalarField& hi = hs[i];
        hSpeciesCorrection += fvc::laplacian(thermo.mu()/Pr_mix*hi, Yi);
        hSpeciesDiffusion += fvc::div(J_[i], hi, "div(J_k,hs_k)");
    }

    hSpeciesCorrection.correctBoundaryConditions();
    hSpeciesDiffusion.correctBoundaryConditions();
    /* End of additional calculations*/

    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi + phi_corr, Yi)
              - fvm::laplacian(thermo.mu()/(Le_k[i] * Pr_mix), Yi)
	          - fvm::div(fvc::interpolate((thermo.mu()/(Le_k[i] * Pr_mix))/WMean)*gradWMean, Yi)
             ==
                reaction->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
