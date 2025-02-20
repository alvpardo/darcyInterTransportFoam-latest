if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        mesh.newIOobject("rhoPhiSum"),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    surfaceScalarField rhoCpPhiSum = 0.0*rhoCpPhi;

    tmp<volScalarField> trSubDeltaT;

    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
    }

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        rhoCpPhiSum += (runTime.deltaT()/totalDeltaT)*rhoCpPhi;
    }

    rhoPhi = rhoPhiSum;
    rhoCpPhi = rhoCpPhiSum;
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;
rhoCp == alpha1*rho1*Cp1 + alpha2*rho2*Cp2;
