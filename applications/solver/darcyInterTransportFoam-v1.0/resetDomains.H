dictionary resetOptions(IFstream("constant/resetOptions")());

// Whether to reset all the surface domain fields
bool resetSur
(
    resetOptions.getOrDefault<Switch>("resetSur",false)
);

// Whether to reset all the subsurface domain fields
bool resetSub
(
    resetOptions.getOrDefault<Switch>("resetSub",false)
);

// Whether to reset the heat transfer fields at the surface domain
bool resetHeatSur
(
    resetOptions.getOrDefault<Switch>("resetHeatSur",false)
);

// Whether to reset the heat transfer fields at the subsurface domain
bool resetHeatSub
(
    resetOptions.getOrDefault<Switch>("resetHeatSub",false)
);

// Whether to reset the solute transport fields at the surface domain
bool resetSoluteSur
(
    resetOptions.getOrDefault<Switch>("resetSoluteSur",false)
);

// Whether to reset the solute transport fields at the subsurface domain
bool resetSoluteSub
(
    resetOptions.getOrDefault<Switch>("resetSoluteSub",false)
);

bool reset = false;
word initialTime;
if (resetSur or resetSub or resetHeatSur or resetHeatSub or resetSoluteSur or resetSoluteSub)
{
    reset = true;

    Info << "Reading resetOptions:" << endl;

    // List of simulated times
    instantList times = runTime.times();

    if (resetOptions.found("initialTime"))
    {
        // Initial time of the new simulation
        scalar initialTime_
        (
            resetOptions.get<scalar>("initialTime")
        );

        initialTime = name(initialTime_);

        // Whether to remove the timesteps between the selected initial time
        // and the latest
        bool rmInterTimes
        (
            resetOptions.getOrDefault<Switch>("rmInterTimes",false)
        );

        // Removal of the timesteps between the selected initial time and the latest
        for (label t=times.size()-1; t>=Time::findClosestTimeIndex(times,initialTime_); t--)
        {
            runTime.setTime(times[t], 0);
            if (rmInterTimes and t != Time::findClosestTimeIndex(times,initialTime_))
            {
                rmDir(mesh.time().path()/runTime.timeName());
            }
        }
    }
    else
    {
        // Selection of the latest simulated time
        runTime.setTime(times[times.size() - 1], 0);
        initialTime = runTime.timeName();
    }

    Info << "    Initial time of the new simulation: " << initialTime << endl;

    // Simulation time from which the reset values are taken
    scalar resetTime_
    (
        resetOptions.getOrDefault<scalar>("resetTime",0)
    );

    word resetTime = name(resetTime_);

    Info << "    Simulation time from which the reset values are taken: " << resetTime << endl;

    // Temporary copy of the initial simulated time of the surface domain
    if (resetSur or resetHeatSur or resetSoluteSur)
    {
        cp(mesh.time().path()/initialTime/"surface", mesh.time().path()/"." + initialTime/"surface");
    }

    // Temporary copy of the initial simulated time of the subsurface domain
    if (resetSub or resetHeatSub or resetSoluteSub)
    {
        cp(meshSub.time().path()/initialTime/"subsurface", meshSub.time().path()/"." + initialTime/"subsurface");
    }


    // Reset of all the surface domain fields
    if (resetSur)
    {
        Info << "    Reset of all the surface domain fields" << endl;
        rmDir(mesh.time().path()/initialTime/"surface");
        cp(mesh.time().path()/resetTime/"surface", mesh.time().path()/initialTime/"surface");
    }

    // Reset of all the subsurface domain fields
    if (resetSub)
    {
        Info << "    Reset of all the subsurface domain fields" << endl;
        rmDir(meshSub.time().path()/initialTime/"subsurface");
        cp(meshSub.time().path()/resetTime/"subsurface", meshSub.time().path()/initialTime/"subsurface");
    }


    // Reset of the heat transfer fields at the surface domain
    if (resetHeatSur)
    {
        Info << "    Reset of the heat transfer fields at the surface domain" << endl;
        cp(mesh.time().path()/resetTime/"surface"/"T", mesh.time().path()/initialTime/"surface"/"T");
    }

    // Reset of the heat transfer fields at the subsurface domain
    if (resetHeatSub)
    {
        Info << "    Reset of the heat transfer fields at the subsurface domain" << endl;
        cp(meshSub.time().path()/resetTime/"subsurface"/"T", meshSub.time().path()/initialTime/"subsurface"/"T");
    }


    // Reset of the solute transport fields at the surface domain
    if (resetSoluteSur)
    {
        Info << "    Reset of the solute transport fields at the surface domain" << endl;
        dictionary soluteSurProp(IFstream("constant/surface/transportProperties")());
        word soluteNameSur(soluteSurProp.getOrDefault<word>("soluteName","C"));

        cp(mesh.time().path()/resetTime/"surface"/soluteNameSur, mesh.time().path()/initialTime/"surface"/soluteNameSur);
        cp(mesh.time().path()/resetTime/"surface"/"D_" + soluteNameSur, mesh.time().path()/initialTime/"surface"/"D_" + soluteNameSur);

        // phiC file is removed because it is not always present in the time directories
        rm(mesh.time().path()/initialTime/"surface"/"phiC");
        cp(mesh.time().path()/resetTime/"surface"/"phiC", mesh.time().path()/initialTime/"surface"/"phiC");
    }

    // Reset of the solute transport fields at the subsurface domain
    if (resetSoluteSub)
    {
        Info << "    Reset of the solute transport fields at the subsurface domain" << endl;
        dictionary soluteSubProp(IFstream("constant/subsurface/transportProperties")());
        word soluteNameSub(soluteSubProp.getOrDefault<word>("soluteName","C"));

        cp(meshSub.time().path()/resetTime/"subsurface"/soluteNameSub, meshSub.time().path()/initialTime/"subsurface"/soluteNameSub);
        cp(meshSub.time().path()/resetTime/"subsurface"/"D_" + soluteNameSub, meshSub.time().path()/initialTime/"subsurface"/"D_" + soluteNameSub);

        // phiC file is removed because it is not always present in the time directories
        rm(meshSub.time().path()/initialTime/"subsurface"/"phiC");
        cp(meshSub.time().path()/resetTime/"subsurface"/"phiC", meshSub.time().path()/initialTime/"subsurface"/"phiC");
    }

    Info << "\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" << endl;
}
