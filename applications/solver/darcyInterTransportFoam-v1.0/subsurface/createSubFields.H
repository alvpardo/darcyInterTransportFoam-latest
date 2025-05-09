Info << "Subsurface domain setup\n" << endl;

// Hydraulic head (L)
Info << "Reading the subsurface flow field h [L]" << endl;
volScalarField hSub
(
    IOobject
    (
        "h",
        runTime.timeName(),
        meshSub,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSub
);

meshSub.setFluxRequired(hSub.name());

// Hydraulic conductivity (L/T)
Info << "Reading the subsurface flow field K [L/T]" << endl;
volTensorField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        meshSub,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSub
);

// Keeping the original values of K for latter assignment
// of impervious regions within the aquifer
volTensorField K_orig(K);

// Porosity (-)
Info << "Reading the subsurface flow field theta [-]\n" << endl;
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    meshSub,
    dimensionedScalar(dimless,1)
);

// Keeping the original values of theta for latter assignment
// of impervious regions within the aquifer
volScalarField theta_orig(theta);

// Specific discharge (L/T)
volVectorField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    meshSub,
    dimensionedVector(dimVelocity,vector::zero)
);

// Porewater/average linear velocity (=specific discharge/porosity) (L/T)
volVectorField USub
(
    IOobject
    (
        "U",
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    meshSub,
    dimensionedVector(dimVelocity,vector::zero)
);

#include "createSubPhi.H"

// Static pressure (pressure head) (M/L*T^-2)
volScalarField pSub
(
    IOobject
    (
        "p",
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    meshSub,
    dimensionedScalar(dimensionSet(1,-1,-2,0,0,0,0),Zero)
);

// Water depth (L)
volScalarField dSub
(
    IOobject
    (
        "d",
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    meshSub,
    dimensionedScalar(dimLength,Zero)
);


// Reading the subsurface transportProperties dictionary
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        meshSub,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Whether to solve flow in the subsurface domain
bool flowSub
(
    transportProperties.getOrDefault<Switch>("flow",false)
);

Info << "Single simulation of flow at the subsurface domain: " << flowSub << "\n" << endl;

// Whether to apply the surface pressure gradient in the subsurface domain
bool pGradSur
(
    transportProperties.getOrDefault<Switch>("pGradSur",false)
);

Info << "Application of the surface pressure gradient in the subsurface domain: " << pGradSur << "\n" << endl;

// Multiplier for gradP
scalar gradPMulti = 0.0;
if (pGradSur) gradPMulti = 1.0;

// Predicted velocity from gradP, not the pressure distribution yet (L/T)
volVectorField UgradP
(
    IOobject
    (
        "UgradP",
        runTime.timeName(),
        meshSub,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshSub,
    dimensionedVector(dimVelocity,vector::zero)
);

// Flow flux due to the predicted velocity (L^3/T)
surfaceScalarField phiUgradP
(
    IOobject
    (
        "phiUgradP",
        runTime.timeName(),
        meshSub,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshSub,
    dimensionedScalar(dimensionSet(0,3,-1,0,0,0,0),Zero)
);


// Reading the imperRegionProperties dictionary
IOdictionary imperRegionProperties
(
    IOobject
    (
        "imperRegionProperties",
        runTime.constant(),
        meshSub,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Whether to generate impervious regions within the subsurface domain
bool imperRegion
(
    imperRegionProperties.getOrDefault<Switch>("imperRegion",false)
);

// Simulation end time
scalar endTime
(
    readScalar(runTime.controlDict().lookup("endTime"))
);

// Simulation end time to update the extension of the subsurface impervious regions
scalar endTimeUpdateRegion_
(
    imperRegionProperties.getOrDefault<scalar>("endTimeUpdateRegion",endTime)
);

word endTimeUpdateRegion = name(endTimeUpdateRegion_);

// Property/ies used to define the subsurface impervious regions. Either "all", "K" or "theta"
word imperProperty
(
    imperRegionProperties.getOrDefault<word>("imperProperty","all")
);

// Depth of the impervious regions within the subsurface domain. Must be negative (L)
dimensionedScalar d_imper
(
    dimLength,imperRegionProperties.get<scalar>("d_imper")
);

// Hydraulic conductivity of the impervious regions within the subsurface domain (L/T)
dimensionedScalar K_imper
(
    dimVelocity,imperRegionProperties.getOrDefault<scalar>("K_imper",1e-20)
);

// Porosity of the impervious regions within the subsurface domain (-)
dimensionedScalar theta_imper
(
    dimless,imperRegionProperties.getOrDefault<scalar>("theta_imper",1e-3)
);

// Solute concentration in the impervious regions of the subsurface domain (M/L^3)
dimensionedScalar CSub_imper
(
    dimDensity,imperRegionProperties.getOrDefault<scalar>("CSub_imper",0.0)
);

// Longitudinal dispersivity of the impervious regions within the subsurface domain (L)
dimensionedScalar alpha_L_imper
(
    dimLength,imperRegionProperties.getOrDefault<scalar>("alpha_L_imper",1e-7)
);

// Transverse dispersivity of the impervious regions within the subsurface domain (L)
dimensionedScalar alpha_T_imper
(
    dimLength,imperRegionProperties.getOrDefault<scalar>("alpha_T_imper",1e-7)
);

// Water temperature in the impervious regions of the subsurface domain (THETA)
dimensionedScalar TSub_imper
(
    dimTemperature,imperRegionProperties.getOrDefault<scalar>("TSub_imper",0.0)
);

if (imperRegion)
{
    Info << "Subsurface impervious region properties:" << endl;
    Info << "    Simulation end time to update the extension of the subsurface impervious regions = " << endTimeUpdateRegion << " seconds" << endl;
    Info << "    Property/ies used to define the subsurface impervious regions = " << imperProperty << endl;
    Info << "    Depth of the impervious regions within the subsurface domain = " << d_imper.value() << " [L]" << endl;
    Info << "    Hydraulic conductivity tensor of the impervious regions within the subsurface domain = " << K_imper.value()*tensor(I) << " [L/T]" << endl;
    Info << "    Porosity of the impervious regions within the subsurface domain = " << theta_imper.value() << " [-]" << endl;
    Info << "    Solute concentration in the impervious regions of the subsurface domain = " << CSub_imper.value() << " [M/L^3]" << endl;
    Info << "    Longitudinal dispersivity of the impervious regions within the subsurface domain = " << alpha_L_imper.value() << " [L]" << endl;
    Info << "    Transverse dispersivity of the impervious regions within the subsurface domain = " << alpha_T_imper.value() << " [L]" << endl;
    Info << "    Water temperature in the impervious regions of the subsurface domain = " << TSub_imper.value() << " [THETA]\n" << endl;
}

// Cell centres of the subsurface domain
pointField cellCtr_meshSub = meshSub.cellCentres();

// Set to zero the third dimension of cellCtr_meshSub
forAll (cellCtr_meshSub,i)
{
    cellCtr_meshSub[i].component(vector::Z) *= 0.0;
}

// Face centres of the top patch of the subsurface domain
fvPatchVectorField faceCtr_subTop = meshSub.C().boundaryField()[subSurfaceTopPatchID];

// Set to zero the third dimension of faceCtr_subTop
forAll (faceCtr_subTop,i)
{
    faceCtr_subTop[i].component(vector::Z) *= 0.0;
}

// Face centres of the bottom patch of the subsurface domain
fvPatchVectorField faceCtr_subBottom = meshSub.C().boundaryField()[subSurfaceBottomPatchID];

// Set to zero the third dimension of faceCtr_subBottom
forAll (faceCtr_subBottom,i)
{
    faceCtr_subBottom[i].component(vector::Z) *= 0.0;
}


// Whether to solve heat transfer in the subsurface domain
bool heatTransferSub
(
    transportProperties.getOrDefault<Switch>("heatTransfer",false)
);

Info << nl << "Single simulation of heat transfer at the subsurface domain: " << heatTransferSub << "\n" << endl;

// Number of TSub correcting iterations
scalar nCorrTSub
(
    transportProperties.getOrDefault<scalar>("nCorrT",0)
);

Info << "Number of subsurface temperature correcting iterations: " << nCorrTSub << "\n" << endl;

// Subsurface water temperature (THETA)
Info << "Reading the subsurface heat field T [THETA]" << endl;
volScalarField TSub
(
    IOobject
    (
        "T",
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    meshSub,
    dimensionedScalar(dimTemperature,Zero)
);

// Keeping the original values of TSub for latter assignment
// of impervious regions within the aquifer
volScalarField TSub_orig(TSub);

// Subsurface media
word medium1Name(transportProperties.get<wordList>("media")[0]);
word medium2Name(transportProperties.get<wordList>("media")[1]);

dictionary medium1=transportProperties.subDict(medium1Name);
dictionary medium2=transportProperties.subDict(medium2Name);

// Kinematic viscosity of the subsurface fluid fraction (L^2/T)
dimensionedScalar nu_f
(
    dimViscosity,medium1.getOrDefault<scalar>("nu",nu1.value())
);

Info << "Reading the kinematic viscosity of the subsurface fluid fraction = " << nu_f.value() << " [L^2/T]" << endl;

// Density of the subsurface fluid fraction (M/L^3)
dimensionedScalar rho_f
(
    dimDensity,medium1.getOrDefault<scalar>("rho",rho1.value())
);

Info << "Reading the density of the subsurface fluid fraction = " << rho_f.value() << " [M/L^3]" << endl;

// Density of the subsurface solid fraction (M/L^3)
dimensionedScalar rho_s
(
    dimDensity,medium2.get<scalar>("rho")
);

Info << "Reading the density of the subsurface solid fraction = " << rho_s.value() << " [M/L^3]" << endl;

// Specific heat capacity of the subsurface fluid fraction (L^2/T^2*C)
dimensionedScalar Cp_f
(
    dimensionSet(0,2,-2,-1,0,0,0),medium1.getOrDefault<scalar>("Cp",Cp1.value())
);

Info << "Reading the specific heat capacity of the subsurface fluid fraction = " << Cp_f.value() << " [L^2/T^2*C]" << endl;

// Specific heat capacity of the subsurface solid fraction (L^2/T^2*C)
dimensionedScalar Cp_s
(
    dimensionSet(0,2,-2,-1,0,0,0),medium2.get<scalar>("Cp")
);

Info << "Reading the specific heat capacity of the subsurface solid fraction = " << Cp_s.value() << " [L^2/T^2*C]" << endl;

// Prandtl number of the subsurface fluid fraction (-)
dimensionedScalar Pr_f
(
    dimless,medium1.getOrDefault<scalar>("Pr",Pr1.value())
);

Info << "Reading the Prandtl number of the subsurface fluid fraction = " << Pr_f.value() << " [-]" << endl;

// Thermal conductivity of the subsurface fluid fraction (M*L/T^3*C)
dimensionedScalar kappa_f
(
    dimensionSet(1,1,-3,-1,0,0,0),medium1.getOrDefault<scalar>("kappa",Cp_f.value()*nu_f.value()*rho_f.value()/Pr_f.value())
);

Info << "Reading the thermal conductivity of the subsurface fluid fraction = " << kappa_f.value() << " [M*L/T^3*C]" << endl;

// Thermal conductivity of the subsurface solid fraction (M*L/T^3*C)
dimensionedScalar kappa_s
(
    dimensionSet(1,1,-3,-1,0,0,0),medium2.get<scalar>("kappa")
);

Info << "Reading the thermal conductivity of the subsurface solid fraction = " << kappa_s.value() << " [M*L/T^3*C]\n" << endl;


// Whether to solve solute transport in the subsurface domain
bool soluteTransportSub
(
    transportProperties.getOrDefault<Switch>("soluteTransport",false)
);

Info << nl << "Single simulation of solute transport at the subsurface domain: " << soluteTransportSub << "\n" << endl;

// Subsurface transport field name
word soluteNameSub
(
    transportProperties.getOrDefault<word>("soluteName","C")
);

// Number of CSub correcting iterations
scalar nCorrCSub
(
    transportProperties.getOrDefault<scalar>("nCorrC",0)
);

Info << "Number of subsurface solute correcting iterations: " << nCorrCSub << "\n" << endl;

// Subsurface solute concentration (M/L^3)
Info << "Reading the subsurface solute field " << soluteNameSub << " [M/L^3]" << endl;
volScalarField CSub
(
    IOobject
    (
        soluteNameSub,
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    meshSub,
    dimensionedScalar(dimDensity,Zero)
);

// Keeping the original values of CSub for latter assignment
// of impervious regions within the aquifer
volScalarField CSub_orig(CSub);

// Hydrodynamic dispersion coefficient (L^2/T)
Info << "Reading the hydrodynamic dispersion coefficient field D_" << soluteNameSub << " [L^2/T]" << endl;
volTensorField DSub
(
    IOobject
    (
        "D_" + soluteNameSub,
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    meshSub,
    dimensionedTensor(phi.dimensions()/dimLength,tensor(I))
);

// Constant hydrodynamic dispersion coefficient (L^2/T)
dimensionedScalar constDSub
(
    phi.dimensions()/dimLength,Zero
);

// Small velocity for division
dimensionedScalar USmall
(
    dimVelocity,VSMALL
);

// Subsurface molecular diffusion coefficient (L^2/T)
dimensionedScalar molDSub
(
    phi.dimensions()/dimLength,transportProperties.getOrDefault<scalar>("molDSub",molDSur.value())
);

Info << "Reading the subsurface molecular diffusion coefficient = " << molDSub.value() << " [L^2/T]" << endl;

// Longitudinal dispersivity (L)
Info << "Reading the longitudinal dispersivity field [L]" << endl;
volScalarField alpha_L
(
    IOobject
    (
        "alpha_L",
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    meshSub,
    dimensionedScalar(dimLength,Zero)
);

// Keeping the original values of alpha_L for latter assignment
// of impervious regions within the aquifer
volScalarField alpha_L_orig(alpha_L);

// Transverse dispersivity (L)
Info << "Reading/calculating the transverse dispersivity field [L]\n" << endl;
volScalarField alpha_T
(
    IOobject
    (
        "alpha_T",
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.1*alpha_L
);

// Keeping the original values of alpha_T for latter assignment
// of impervious regions within the aquifer
volScalarField alpha_T_orig(alpha_T);

Info << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" << endl;
