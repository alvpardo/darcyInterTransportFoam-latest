Info << nl << "Solving subsurface flow ..." << endl;

// Read gradP calculated from surface domain
uniformDimensionedVectorField gradP
(
    IOobject
    (
        "deltaTgradP",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedVector(dimensionSet(1,-2,-2,0,0,0,0),vector::zero)
);

Info << "gradP = " << gradP.value() << endl;

// Predicted velocity from gradP, not the pressure distribution yet
UgradP == gradPMulti*(K&gradP/(rho1*mag(g)));

// Flow flux due to the predicted velocity
phiUgradP = fvc::flux(UgradP);

// Non-orthogonal corrector loop
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix hSubEqn
    (
        fvm::laplacian(K,hSub) == fvc::div(phiUgradP)
    );

    hSubEqn.solve();

    if (simple.finalNonOrthogonalIter())
    {
        // Darcy flux (returned from hSubEqn.flux())
        phi_q = phiUgradP - hSubEqn.flux();

        hSub.relax();

        // Specific discharge (L/T)
        q = UgradP - (K&fvc::grad(hSub));
        q.correctBoundaryConditions();
    }
}

// Porewater flux (needed in solute transport equation)
phi_USub = phi_q/linearInterpolate(theta);

Info<< "Continuity error of phi_USub= "
    << mag(fvc::div(phi_USub))().weightedAverage(meshSub.V()).value() << endl;

// Porewater/average linear velocity (=specific discharge/porosity) (L/T)
USub == q/theta;

// Static pressure (pressure head) (M/L*T^-2)
pSub == rho1*mag(g)*(hSub-meshSub.C().component(vector::Z));

// Water depth (L)
dSub == pSub/(rho1*mag(g));
