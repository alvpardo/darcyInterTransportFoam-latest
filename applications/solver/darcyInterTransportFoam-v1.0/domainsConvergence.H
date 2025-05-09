// Flow convergence check
if (flow or (flowSur and flowSub))
{
    scalarField phiBottomSur = phi.boundaryField()[surfaceBottomPatchID];
    scalarField phiTopSub = phi_q.boundaryField()[subSurfaceTopPatchID];

    scalarField couplePhiDiff = 0.*phiBottomSur;

    forAll (couplePhiDiff,i)
    {
        couplePhiDiff[i] = phiBottomSur[i] - phiTopSub[i];
    }

    couplePhiMaxMismatch = gMax(mag(couplePhiDiff));
}
else
{
    couplePhiMaxMismatch = 0.0;
}

if (couplePhiMaxMismatch > couplePhiTol)
{
    Info << "Maximum coupled flow phi mismatch: " << couplePhiMaxMismatch << " > " << couplePhiTol << endl;
}
else
{
    Info << "Maximum coupled flow phi mismatch: " << couplePhiMaxMismatch << " < " << couplePhiTol << " OK" << endl;
}


// Advective scalar convergence check
if (soluteTransport or (soluteTransportSur and soluteTransportSub))
{
    scalarField solutePhiBottomSur = phiCSur.boundaryField()[surfaceBottomPatchID];
    scalarField solutePhiTopSub = phiCSub.boundaryField()[subSurfaceTopPatchID];

    scalarField coupleSolutePhiDiff = 0.*solutePhiBottomSur;

    forAll (coupleSolutePhiDiff,i)
    {
        coupleSolutePhiDiff[i] = solutePhiBottomSur[i] - solutePhiTopSub[i];
    }

    coupleSolutePhiMaxMismatch = gMax(mag(coupleSolutePhiDiff));
}
else
{
    coupleSolutePhiMaxMismatch = 0.0;
}

if (coupleSolutePhiMaxMismatch > couplePhiSoluteTol)
{
    Info << "Maximum coupled solute phi mismatch: " << coupleSolutePhiMaxMismatch << " > " << couplePhiSoluteTol << "\n" << endl;
}
else
{
    Info << "Maximum coupled solute phi mismatch: " << coupleSolutePhiMaxMismatch << " < " << couplePhiSoluteTol << " OK\n" << endl;
}
