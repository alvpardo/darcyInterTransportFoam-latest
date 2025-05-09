#include "createRDeltaT.H"

Info << "Surface domain setup\n" << endl;

// Modified pressure = static pressure - hydrostatic pressure (M/L*T^-2)
// p_rgh = p - rho*(g & z)
Info << "Reading the surface flow field p_rgh [M/L*T^-2]" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Flow velocity (L/T)
Info << "Reading the surface flow field U [L/T]\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createSurPhi.H"

// Static pressure (M/L*T^-2)
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(1,-1,-2,0,0,0,0),Zero)
);

// Water depth (L)
volScalarField dSur
(
    IOobject
    (
        "d",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength,Zero)
);

// Piezometric level (L)
volScalarField piezoSur
(
    IOobject
    (
        "piezo",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength,Zero)
);

// Hydraulic head (L)
volScalarField hSur
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength,Zero)
);


// Reading the surface transportProperties dictionary
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

// Whether to solve flow in the surface domain
bool flowSur
(
    mixture.getOrDefault<Switch>("flow",false)
);

Info << "\nSingle simulation of flow at the surface domain: " << flowSur << "\n" << endl;

// Whether to clip the air fraction
bool alphaAirClip
(
    mixture.getOrDefault<Switch>("alphaAirClip",false)
);

// Water fraction threshold
scalar alphaWaterThres
(
    mixture.getOrDefault<scalar>("alphaWaterThres",0.5)
);

if (rigidWaterSurface)
{
    Info << "Rigid water surface " << rigidWaterSurface << "\n" << endl;
}
else
{
    Info << "Clipping of the air fraction: " << alphaAirClip << endl;
    Info << "Water fraction threshold: " << alphaWaterThres << "\n" << endl;
}

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

// Density of phase1 (M/L^3)
const dimensionedScalar& rho1 = mixture.rho1();

// Density of phase2 (M/L^3)
const dimensionedScalar& rho2 = mixture.rho2();

// Need to store rho for ddt(rho, U) (M/L^3)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Mass flux (M/T)
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi // By default. It provides the same results without initializing it
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());


// Hydrostatic head at the inlet and outlet patches of the surface domain (L)
IOdictionary inOutStatic_hSur
(
    IOobject
    (
        "inOutStatic_h",
        runTime.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
);

// Face centres of the bottom patch of the surface domain
fvPatchVectorField faceCtr_surBottom = mesh.C().boundaryField()[surfaceBottomPatchID];

// Set to zero the third dimension of faceCtr_surBottom
forAll (faceCtr_surBottom,i)
{
    faceCtr_surBottom[i].component(vector::Z) *= 0.0;
}

// Dry cells at the bottom patch of the surface domain (-)
volScalarField dryBottomSur
(
    IOobject
    (
        "dryBottom",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless,Zero)
);

// Keeping the original values of dryBottomSur for the latter
// definition of the dry regions at the surface bottom patch
volScalarField dryBottomSur_orig(dryBottomSur);

// Whether to solve heat transfer in the surface domain
bool heatTransferSur
(
    mixture.getOrDefault<Switch>("heatTransfer",false)
);

Info << nl << "Single simulation of heat transfer at the surface domain: " << heatTransferSur << "\n" << endl;

// Surface phases where to solve heat transfer. Either "all" or "phase1"
word heatPhases
(
    mixture.getOrDefault<word>("heatPhases","phase1")
);

Info << "Surface phases where to solve heat transfer: " << heatPhases << endl;

// Number of TSur correcting iterations
scalar nCorrTSur
(
    mixture.getOrDefault<scalar>("nCorrT",0)
);

Info << "Number of surface temperature correcting iterations: " << nCorrTSur << "\n" << endl;

// Surface water temperature (THETA)
Info << "Reading the surface heat field T [THETA]" << endl;
volScalarField TSur
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimTemperature,Zero)
);

// Surface phases
word phase1Name(mixture.get<wordList>("phases")[0]);
word phase2Name(mixture.get<wordList>("phases")[1]);

// Specific heat capacity of phase1 (L^2/T^2*C)
const dimensionedScalar& Cp1 = mixture.Cp1();
Info << "Reading the specific heat capacity of " << phase1Name << " = " << Cp1.value() << " [L^2/T^2*C]" << endl;

// Specific heat capacity of phase2 (L^2/T^2*C)
const dimensionedScalar& Cp2 = mixture.Cp2();
Info << "Reading the specific heat capacity of " << phase2Name << " = " << Cp2.value() << " [L^2/T^2*C]" << endl;

// Prandtl number of phase1 (-)
const dimensionedScalar& Pr1 = mixture.Pr1();
Info << "Reading the Prandtl number of " << phase1Name << " = " << Pr1.value() << " [-]" << endl;

// Prandtl number of phase2 (-)
const dimensionedScalar& Pr2 = mixture.Pr2();
Info << "Reading the Prandtl number of " << phase2Name << " = " << Pr2.value() << " [-]" << endl;

dictionary phase1=mixture.subDict(phase1Name);

// Kinematic viscosity of phase1 (L^2/T)
dimensionedScalar nu1
(
    dimViscosity,phase1.get<scalar>("nu")
);

Info << "Reading the kinematic viscosity of " << phase1Name << " = " << nu1.value() << " [L^2/T]\n" << endl;

// Need to store rhoCp for ddt(rhoCp, T) (M/L*T^2*K)
volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1*Cp1 + alpha2*rho2*Cp2,
    alpha1.boundaryField().types()
);
rhoCp.oldTime();

// Heat flux (M*L^2/T^3*K)
surfaceScalarField rhoCpPhi
(
    IOobject
    (
        "rhoCpPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rhoCp)*phi
);


// Whether to solve solute transport in the surface domain
bool soluteTransportSur
(
    mixture.getOrDefault<Switch>("soluteTransport",false)
);

Info << nl << "Single simulation of solute transport at the surface domain: " << soluteTransportSur << "\n" << endl;

// Surface phases where to solve solute transport. Either "all" or "phase1"
word solutePhases
(
    mixture.getOrDefault<word>("solutePhases","phase1")
);

Info << "Surface phases where to solve solute transport: " << solutePhases << endl;

// Surface transport field name
word soluteNameSur
(
    mixture.getOrDefault<word>("soluteName","C")
);

// Number of CSur correcting iterations
scalar nCorrCSur
(
    mixture.getOrDefault<scalar>("nCorrC",0)
);

Info << "Number of surface solute correcting iterations: " << nCorrCSur << "\n" << endl;

// Surface solute concentration (M/L^3)
Info << "Reading the surface solute field " << soluteNameSur << " [M/L^3]" << endl;
volScalarField CSur
(
    IOobject
    (
        soluteNameSur,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimDensity,Zero)
);

mesh.setFluxRequired(CSur.name());

// Surface diffusion coefficient (L^2/T)
Info << "Reading the surface diffusion coefficient field D_" << soluteNameSur << " [L^2/T]" << endl;
volScalarField DSur
(
    IOobject
    (
        "D_" + soluteNameSur,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions()/dimLength,Zero)
);

// Constant surface diffusion coefficient (L^2/T)
dimensionedScalar constDSur
(
    phi.dimensions()/dimLength,Zero
);

// Surface molecular diffusion coefficient - physical diffusivity (L^2/T)
dimensionedScalar molDSur
(
    phi.dimensions()/dimLength,mixture.getOrDefault<scalar>("molDSur",Zero)
);

Info << "Reading the surface molecular diffusion coefficient = " << molDSur.value() << " [L^2/T]" << endl;

// Reading the surface turbulenceProperties dictionary
IOdictionary turbulenceProperties
(
    IOobject
    (
        "turbulenceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Turbulence model type
word turbModelType
(
    turbulenceProperties.get<word>("simulationType")
);

// Turbulent Schmidt number (-)
dimensionedScalar Sct
(
    dimless,mixture.getOrDefault<scalar>("Sct",1.0)
);

Info << "Reading the turbulent Schmidt number = " << Sct.value() << " [-]" << endl;

// Henry's law solubility constant (-)
dimensionedScalar Hcc
(
    dimless,mixture.getOrDefault<scalar>("Hcc",1.0)
);

Info << "Reading the Henry's law solubility constant = " << Hcc.value() << " [-]\n" << endl;

#include "createMRF.H"
#include "createFvOptions.H"

Info << "\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" << endl;
