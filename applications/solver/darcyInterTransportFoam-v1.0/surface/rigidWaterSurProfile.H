// Shape of the rigid water surface
word waterShape
(
    mixture.getOrDefault<word>("waterShape","flat")
);

// Mean water depth (L)
dimensionedScalar dMean
(
    dimLength,mixture.getOrDefault<scalar>("dMean",0.0)
);

IOdictionary waterSurfaceProfileCoeff
(
    IOobject
    (
        "waterSurfaceProfileCoeff",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

surfaceScalarField waterSurfaceProfile
(
    IOobject
    (
        "waterSurfaceProfile",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength,Zero)
);

    // Cell centre Z coordinates of the surface domain
    volScalarField zSur
    (
        IOobject
        (
            "zSur",
            runTime.timeName(),
            meshSub,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh.C().component(vector::Z)
    );

    // Z coordinates of the surface bottom patch
    fvPatchScalarField zSurBottom = zSur.boundaryField()[surfaceBottomPatchID];

if (rigidWaterSurface)
{
    // Coefficient a0
    dimensionedScalar a0
    (
        dimLength,waterSurfaceProfileCoeff.get<scalar>("a0")
    );

    // Coefficient a1
    dimensionedScalar a1
    (
        dimLength,waterSurfaceProfileCoeff.get<scalar>("a1")
    );

    // Coefficient b1
    dimensionedScalar b1
    (
        dimLength,waterSurfaceProfileCoeff.get<scalar>("b1")
    );

    // Coefficient a2
    dimensionedScalar a2
    (
        dimLength,waterSurfaceProfileCoeff.get<scalar>("a2")
    );

    // Coefficient b2
    dimensionedScalar b2
    (
        dimLength,waterSurfaceProfileCoeff.get<scalar>("b2")
    );

    // Coefficient a3
    dimensionedScalar a3
    (
        dimLength,waterSurfaceProfileCoeff.get<scalar>("a3")
    );

    // Coefficient b3
    dimensionedScalar b3
    (
        dimLength,waterSurfaceProfileCoeff.get<scalar>("b3")
    );

    // Coefficient a4
    dimensionedScalar a4
    (
        dimLength,waterSurfaceProfileCoeff.get<scalar>("a4")
    );

    // Coefficient b4
    dimensionedScalar b4
    (
        dimLength,waterSurfaceProfileCoeff.get<scalar>("b4")
    );

    // Coefficient w
    dimensionedScalar w
    (
        dimLength,waterSurfaceProfileCoeff.get<scalar>("w")
    );

    // Cell centre X coordinates of the surface domain
    volScalarField xSur
    (
        IOobject
        (
            "xSur",
            runTime.timeName(),
            meshSub,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh.C().component(vector::X)
    );

    // X coordinates of the surface bottom patch
    fvPatchScalarField xSurBottom = xSur.boundaryField()[surfaceBottomPatchID];

    // Constant, curved free surface (rigid-lid)
    forAll (xSurBottom,i)
    {
        waterSurfaceProfile = a0 + a1*cos(xSurBottom[i]*w) + b1*sin(xSurBottom[i]*w) + a2*cos(2*xSurBottom[i]*w)
            + b2*sin(2*xSurBottom[i]*w) + a3*cos(3*xSurBottom[i]*w) + b3*sin(3*xSurBottom[i]*w) + a4*cos(4*xSurBottom[i]*w) + b4*sin(4*xSurBottom[i]*w);
    }
}
