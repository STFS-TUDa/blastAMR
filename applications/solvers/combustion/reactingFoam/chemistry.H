{
    Info << "Solving chemistry" << endl;

    chemistry.solve
    (
        runTime.value() - runTime.deltaTValue(),
        runTime.deltaTValue()
    );

    // turbulent time scale
    if (turbulentReaction)
    {
        volScalarField tk =
                Cmix*sqrt(turbulence->muEff()/rho/turbulence->epsilon());
        volScalarField tc = chemistry.tc();

        // Chalmers PaSR model
        kappa = (runTime.deltaT() + tc)/(runTime.deltaT() + tc + tk);
    }
    else
    {
        kappa = 1.0;
    }
}
