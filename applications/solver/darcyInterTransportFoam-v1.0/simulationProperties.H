Info << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" << endl;
Info << "Reading simulationProperties:" << endl;

dictionary simulationProperties(IFstream("constant/simulationProperties")());

// Whether to solve flow in both domains. Also possible to be activated separately
// for each domain at the corresponding transportProperties dictionary
// true: flow will be simulated in both domains
// false: flow will not be simulated in any domain (an option for using pre-simulated flow fields)
bool flow
(
    simulationProperties.getOrDefault<Switch>("flow",true)
);

if (flow)
{
    Info << "    Flow will be simulated in both domains" << endl;
}

// Whether to solve heat transfer in both domains. Also possible to be activated
// separately for each domain at the corresponding transportProperties dictionary
bool heatTransfer
(
    simulationProperties.getOrDefault<Switch>("heatTransfer",true)
);

if (heatTransfer)
{
    Info << "    Heat transfer will be simulated in both domains" << endl;
}

// Whether to solve solute transport in both domains. Also possible to be activated
// separately for each domain at the corresponding transportProperties dictionary
bool soluteTransport
(
    simulationProperties.getOrDefault<Switch>("soluteTransport",true)
);

if (soluteTransport)
{
    Info << "    Solute transport will be simulated in both domains" << endl;
}

// Fluid flux tolerance for the coupling of both domains
/*scalar couplePhiTol
(
    simulationProperties.getOrDefault<scalar>("couplePhiTol",1e-3)
);

Info << "    Fluid flux tolerance for the coupling of both domains = " << couplePhiTol << endl;

// Advective scalar flux tolerance for the coupling of both domains
scalar couplePhiSoluteTol
(
    simulationProperties.getOrDefault<scalar>("couplePhiSoluteTol",1e-5)
);

Info << "    Advective scalar flux tolerance for the coupling of both domains = " << couplePhiSoluteTol << endl;*/

// Solute coupling: whether to kill the diffusive flux
bool ignoreDiffusiveFluxCoupling
(
    simulationProperties.getOrDefault<Switch>("ignoreDiffusiveFluxCoupling",false)
);

if (ignoreDiffusiveFluxCoupling)
{
    Info << "    Diffusive flux is ignored at the interface" << endl;
}
else
{
    Info << "    Diffusive flux is not ignored at the interface" << endl;
}

// Name of the interface patch
word interPatchName
(
    simulationProperties.get<word>("interPatchName")
);

Info << "    Interface patch name: " << interPatchName << endl;

// Bottom patch of the surface domain
label surfaceBottomPatchID = mesh.boundaryMesh().findPatchID(interPatchName);

// Top patch of the subsurface domain
label subSurfaceTopPatchID = meshSub.boundaryMesh().findPatchID(interPatchName);

if(subSurfaceTopPatchID == -1 || surfaceBottomPatchID == -1)
{
    FatalErrorIn(args.executable())
        << "Didn't find the patch named " << interPatchName << " in the subsurface domain"
        << " or the patch named " << interPatchName << " in the surface domain"
        << abort(FatalError);
}

// Name of the subsurface bottom patch
word subBottomPatchName
(
    simulationProperties.get<word>("subBottomPatchName")
);

Info << "    Subsurface bottom patch name: " << subBottomPatchName << "\n" << endl;

// Bottom patch of the subsurface domain
label subSurfaceBottomPatchID = meshSub.boundaryMesh().findPatchID(subBottomPatchName);

Info << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" << endl;
