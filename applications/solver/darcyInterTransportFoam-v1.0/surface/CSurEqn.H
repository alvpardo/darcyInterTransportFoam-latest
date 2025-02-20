Info << nl << "Solving surface scalar transport ..." << endl;

// Surface diffusion coefficient =
if (mixture.readIfPresent("constDSur",constDSur))
{
    // = constant surface diffusion coefficient
    DSur == constDSur;
}
else
{
    // = physical diffusivity
    DSur == molDSur;

    if (turbModelType != "laminar")
    {
        // = physical diffusivity + turbulent diffusivity
        DSur += turbulence->nut()/Sct;
    }
}

// Optional: kill the diffusive flux across the SWI
if (ignoreDiffusiveFluxCoupling)
{
    DSur.boundaryFieldRef()[surfaceBottomPatchID] *= 0.0;
}

if (solutePhases == "all")
{
    // If Henry's law solubility constant Hcc=1.0 (default), this term is zero.
    surfaceScalarField phiCiHaroun = ((fvc::interpolate(DSur) * (1.0-Hcc) / 
        (fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*Hcc + SMALL)) * fvc::snGrad(alpha1)) * mesh.magSf();

    for (label i = 0; i <= nCorrCSur; i++)
    {
        fvScalarMatrix CSurEqn
        (
            fvm::ddt(CSur)
          + fvm::div(phi, CSur, "div(phi,C)")
          - fvm::laplacian(DSur, CSur, "laplacian(D,C)")
          + fvm::div(phiCiHaroun, CSur, "div(phi,C)")
          ==
            fvOptions(CSur)
        );

        CSurEqn.relax();
        fvOptions.constrain(CSurEqn);
        CSurEqn.solve();

        // Clip the negative value of concentration
        CSur.max(0.0);
    }

    // Solute mass flux
    phiCSur = phi*fvc::interpolate(CSur);
}
else if (solutePhases == "phase1")
{
    DSur *= pos0(alpha1 - 0.99);

    // Reset D dimensions consistent with alphaPhiUn
    DSur.dimensions().reset(alphaPhiUn.dimensions()/dimLength);

    for (label i = 0; i <= nCorrCSur; i++)
    {
        fvScalarMatrix CSurEqn
        (
            fvm::ddt(CSur)
          + fvm::div(alphaPhiUn, CSur, "div(phi,C)")
          - fvm::laplacian(DSur, CSur, "laplacian(D,C)")
          ==
            fvOptions(alpha1, CSur)
        );

        CSurEqn.relax();
        fvOptions.constrain(CSurEqn);
        CSurEqn.solve();
    }

    // Solute mass flux
    phiCSur = alphaPhiUn*fvc::interpolate(CSur);
}
else
{
    FatalError
        << "Incorrect phase name. Solute transport in the surface domain" << nl
        << "can be solved for either both phases or phase1" << exit(FatalError);
}
