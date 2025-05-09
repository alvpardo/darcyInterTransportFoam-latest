Info << nl << "Solving subsurface heat transfer ..." << endl;

for (label i = 0; i <= nCorrTSub; i++)
{
    fvScalarMatrix TSubEqn
    (
        fvm::ddt(theta*rho_f*Cp_f + (1-theta)*rho_s*Cp_s, TSub)
      + fvm::div(rho_f*Cp_f*phi_USub, TSub, "div(rhoCpPhi,T)")
      - fvm::laplacian(theta*kappa_f + (1-theta)*kappa_s, TSub, "laplacian(kappa,T)")
      ==
        fvOptions(theta*rho_f*Cp_f + (1-theta)*rho_s*Cp_s, TSub)
    );

    TSubEqn.relax();
    fvOptions.constrain(TSubEqn);
    TSubEqn.solve();
}
