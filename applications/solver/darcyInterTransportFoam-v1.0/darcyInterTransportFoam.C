/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

Application
    darcyInterTransportFoam-v1.0

Description
    A solver for hyporheic flow by combining interFoam and porousMediaFoam:
    1. interFoam: Solver for two incompressible, isothermal immiscible fluids
                  using VOF phase-fraction based interface capturing approach.
    2. porousMediaFoam: Solves unsteady saturated porous media flow.

    Fluid is two-way coupled, i.e., the fluid flux across the SWI is 
    continuous. This is important to find conservative fluid fluxes for
    scalar transport.

    Solute transport is added to the solver for both domains. The coupling 
    of solute transport is through the boundary conditions:
    1. surInterfaceSolute
    2. subInterfaceSolute

    Only one conservative scalar is solved. But this can be changed later 
    for cases like non-conservative (reactions) and multispecies.

    Heat transport is added to the solver for both domains. The coupling 
    of heat transport is through the boundary conditions:
    1. surInterfaceHeat
    2. subInterfaceHeat

    Limitations:
        The subsurface domain is assumed to be fully saturated.

Author
    Alvaro Pardo-Alvarez, 09.09.2025
    University of Neuchatel, Centre for Hydrogeology and Geothermics

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture_.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "IFstream.H"

#include "simpleControl.H"

#include "surInterfaceSoluteFvPatchScalarField.H"
#include "subInterfaceSoluteFvPatchScalarField.H"
#include "surInterfaceHeatFvPatchScalarField.H"
#include "subInterfaceHeatFvPatchScalarField.H"
#include "darcyGradHeadFvPatchScalarField.H"
#include "surfaceHydrostaticHeadFvPatchScalarField.H"
#include "zeroNetFluxVelocityFvPatchVectorField.H"
#include "flowRateInletParabolicVelocityFvPatchVectorField.H"

#include "primitivePatchInterpolation.H"
#include "patchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "A solver for hyporheic flow by combining interFoam and porousMediaFoam:\n"
        "1. interFoam: Solver for two incompressible, isothermal immiscible fluids\n"
        "              using VOF phase-fraction based interface capturing approach.\n"
        "2. porousMediaFoam: Solves unsteady saturated porous media flow."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"

    #include "initContinuityErrs.H"
    pimpleControl pimple(mesh);
    simpleControl simple(meshSub);
    #include "createTimeControls.H"

    #include "simulationProperties.H"

    #include "resetDomains.H"

    #include "createSurFields.H"
    #include "createSubFields.H"

    #include "singleDomainSimulation.H"

    #include "createAlphaFluxes.H"
    #include "correctPhi.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    label nDeltaT = 1;
    while (runTime.run())
    {
        if (nDeltaT != 1)
        {
            Info << "\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" << endl;
        }

        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info << "nDeltaT = " << nDeltaT << endl;
        Info << "Time = " << runTime.timeName() << nl << endl;

        // Surface/subsurface coupling iteration loop
        scalar couplePhiMaxMismatch = 1.0;
        scalar coupleSolutePhiMaxMismatch = 1.0;
        label coupleIterNum = 1;
        while (couplePhiMaxMismatch > couplePhiTol and coupleSolutePhiMaxMismatch > couplePhiSoluteTol)
        {
            Info << "Coupling iteration #" << coupleIterNum++ << " ..." << endl;

            Info << nl << "Surface domain ...\n" << endl;

            // Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                if (flow or flowSur)
                {
                    // Solve the surface flow equation
                    #include "interFoam.H"
                }

                if (heatTransfer or heatTransferSur)
                {
                    // Solve the surface heat transfer equation
                    #include "TSurEqn.H"
                }

                if (soluteTransport or soluteTransportSur)
                {
                    // Solve the surface solute transport equation
                    #include "CSurEqn.H"
                }

                Info << endl;
            }

            // Mapping the surface bottom variables to the subsurface top patch
            #include "mappingSurToSub.H"

            Info << "Subsurface domain ..." << endl;

            // Setting of the impervious region within the subsurface domain
            #include "impermeableRegion.H"

            if (flow or flowSub)
            {
                // Solve the subsurface flow equation
                #include "darcyFoam.H"
            }

            if (heatTransfer or heatTransferSub)
            {
                // Solve the subsurface heat transfer equation
                #include "TSubEqn.H"
            }

            if (soluteTransport or soluteTransportSub)
            {
                // Solve the subsurface solute transport equation
                #include "CSubEqn.H"
            }

            Info << endl;

            // Mapping the subsurface top variables to the surface bottom patch
            #include "mappingSubToSur.H"

            // Domains convergence check
            #include "domainsConvergence.H"
        }

        runTime.write();

        runTime.printExecutionTime(Info);

        nDeltaT++;
    }

    Info << "\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" << endl;

    #include "restoreInitialFields.H"

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
