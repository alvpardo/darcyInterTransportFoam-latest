Info << nl << "Solving surface flow ..." << endl;

#include "alphaControls.H"
#include "alphaEqnSubCycle.H"

mixture.correct();

#include "UEqn.H"

// Pressure corrector loop
while (pimple.correct())
{
    #include "pEqn.H"
}

if (pimple.turbCorr())
{
    turbulence->correct();
}

// Water depth (L)
dSur == p/(rho*mag(g));

// Piezometric level (L)
piezoSur == dSur + mesh.C().component(vector::Z);

// Hydraulic head (L)
hSur == (0.5*magSqr(U)/mag(g)) + piezoSur;

#include "patchHydrostaticHead.H"

// Dry faces at the bottom patch of the surface domain (-)
forAll(dryBottomSur,i)
{
    dryBottomSur[i] = dryBottomSur_orig[i];

    if (dryBottomSur[i] != -1)
    {
        dryBottomSur[i] = neg(alpha1[i] - alphaWaterThres);
    }
}

dryBottomSur.correctBoundaryConditions();

if (alphaAirClip)
{
    dSur *= pos0(alpha1 - alphaWaterThres);
    piezoSur *= pos0(alpha1 - alphaWaterThres);
    hSur *= pos0(alpha1 - alphaWaterThres);
}
