{
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
//      - fvm::laplacian(turbulence->muEff(), hs)  // unit lewis no.
     ==
        dpdt
      - (fvc::ddt(rho, K) + fvc::div(phi, K))
      + reaction->Sh()
    );

    hsEqn.relax();
    hsEqn.solve();

    thermo.correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
