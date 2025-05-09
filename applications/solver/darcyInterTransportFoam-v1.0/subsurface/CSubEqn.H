Info << nl << "Solving subsurface scalar transport ..." << endl;

// Hydrodynamic dispersion coefficient =
if (transportProperties.readIfPresent("constDSub",constDSub))
{
    // = constant hydrodynamic dispersion coefficient
    DSub == constDSub*tensor::one;
}
else
{
    // = effective molecular diffusion + mechanical dispersion
    DSub == (molDSub/(1-2*log(theta)))*tensor::one+alpha_T*mag(USub)*tensor(I)+(alpha_L-alpha_T)*((USub*USub)/(mag(USub)+USmall));
}

// Optional: kill the diffusive flux across the SWI
if (ignoreDiffusiveFluxCoupling)
{
    DSub.boundaryFieldRef()[subSurfaceTopPatchID] *= 0.0;
}

for (label i = 0; i <= nCorrCSub; i++)
{
    fvScalarMatrix CSubEqn
    (
        fvm::ddt(CSub)
      + fvm::div(phi_USub, CSub, "div(phi,C)")
      - fvm::laplacian(DSub, CSub, "laplacian(D,C)")
      ==
        fvOptions(CSub)
    );

    CSubEqn.relax();
    fvOptions.constrain(CSubEqn);
    CSubEqn.solve();

    // Clip the negative value of concentration
    CSub.max(0.0);
}

// Solute mass flux
phiCSub = phi_q*fvc::interpolate(CSub);
