Info << nl << "Solving surface heat transfer ..." << endl;

volScalarField kappaSur = mixture.kappa();
//surfaceScalarField kappaSur = mixture.kappaf();
Info << "Average overall kappa of the surface domain: " << average(kappaSur) << endl;

if (heatPhases == "all")
{
    for (label i = 0; i <= nCorrTSur; i++)
    {
        fvScalarMatrix TSurEqn
        (
            fvm::ddt(rhoCp, TSur)
          + fvm::div(rhoCpPhi, TSur, "div(rhoCpPhi,T)")
          - fvm::laplacian(kappaSur, TSur, "laplacian(kappa,T)")
          ==
            fvOptions(rhoCp, TSur)
        );

        TSurEqn.relax();
        fvOptions.constrain(TSurEqn);
        TSurEqn.solve();
    }
}
else if (heatPhases == "phase1")
{
    kappaSur *= pos0(alpha1 - 0.99);

    for (label i = 0; i <= nCorrTSur; i++)
    {
        fvScalarMatrix TSurEqn
        (
            fvm::ddt(rho1*Cp1, TSur)
          + fvm::div(rho1*Cp1*alphaPhi10, TSur, "div(rhoCpPhi,T)")
          - fvm::laplacian(kappaSur, TSur, "laplacian(kappa,T)")
          ==
            fvOptions(rho1*Cp1*alpha1, TSur)
        );

        TSurEqn.relax();
        fvOptions.constrain(TSurEqn);
        TSurEqn.solve();
    }
}
else
{
    FatalError
        << "Incorrect phase name. Heat transfer in the surface domain" << nl
        << "can be solved for either both phases or phase1" << exit(FatalError);
}
