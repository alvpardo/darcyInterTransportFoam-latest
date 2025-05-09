bool single = false;

if (!flow or !(flowSur and flowSub))
{
    Info << "Single domain simulation\n" << endl;
}

// Velocity field at the surface bottom patch
if (flowSur and !flowSub and !flow)
{
    single = true;
    Info << "Flow at the surface domain" << endl;

    // Velocity field at the interface patch (surface bottom patch)
    word interUSur
    (
        mixture.getOrDefault<word>("interUSur","latestTime")
    );

    Info << "Velocity field at the interface patch (surface bottom patch): ";

    if (interUSur == "constU")
    {
        // Constant velocity at the surface bottom patch
        vector Ubottom
        (
            mixture.getOrDefault<vector>("Ubottom",vector::zero)
        );

        U.boundaryFieldRef()[surfaceBottomPatchID] == Ubottom;

        Info << "homogeneous\n" << endl;
    }
    else
    {
        Info << "the one resulting from the surface-subsurface simulation\n" << endl;
    }
}

// Temperature field at the surface bottom patch
if (heatTransferSur and !heatTransferSub and !heatTransfer)
{
    single = true;
    Info << "Heat transfer at the surface domain" << endl;

    // Temperature field at the interface patch (surface bottom patch)
    word interTSur
    (
        mixture.getOrDefault<word>("interTSur","zeroGradient")
    );

    Info << "Temperature field at the interface patch (surface bottom patch): ";

    if (interTSur == "latestTime" or interTSur == "constT")
    {
        scalarField Ttmp_sur(TSur.boundaryField()[surfaceBottomPatchID].size(),Zero);

        if (interTSur == "latestTime")
        {
            Ttmp_sur = TSur.boundaryField()[surfaceBottomPatchID];
            Info << "the one resulting from the surface-subsurface simulation\n" << endl;
        }
        else
        {
            // Constant temperature at the surface bottom patch
            scalar TSurBottom
            (
                mixture.getOrDefault<scalar>("TSurBottom",Zero)
            );

            Ttmp_sur = TSurBottom;
            Info << "homogeneous\n" << endl;
        }

        surInterfaceHeatFvPatchScalarField& tempTSurp =
            dynamic_cast<surInterfaceHeatFvPatchScalarField&>(TSur.boundaryFieldRef()[surfaceBottomPatchID]);

        tempTSurp.TFromSub() = Ttmp_sur;
    }
    else
    {
        Info << "zeroGradient\n" << endl;
    }
}

// Solute field at the surface bottom patch
if (soluteTransportSur and !soluteTransportSub and !soluteTransport)
{
    single = true;
    Info << "Solute transport at the surface domain" << endl;

    // Solute field at the interface patch (surface bottom patch)
    word interCSur
    (
        mixture.getOrDefault<word>("interCSur","zeroGradient")
    );

    Info << "Solute field at the interface patch (surface bottom patch): ";

    if (interCSur == "latestTime" or interCSur == "constC")
    {
        scalarField Ctmp_sur(CSur.boundaryField()[surfaceBottomPatchID].size(),Zero);

        if (interCSur == "latestTime")
        {
            Ctmp_sur = CSur.boundaryField()[surfaceBottomPatchID];
            Info << "the one resulting from the surface-subsurface simulation\n" << endl;
        }
        else
        {
            // Constant solute at the surface bottom patch
            scalar CSurBottom
            (
                mixture.getOrDefault<scalar>("CSurBottom",Zero)
            );

            Ctmp_sur = CSurBottom;
            Info << "homogeneous\n" << endl;
        }

        surInterfaceSoluteFvPatchScalarField& tempCSurp =
            dynamic_cast<surInterfaceSoluteFvPatchScalarField&>(CSur.boundaryFieldRef()[surfaceBottomPatchID]);

        tempCSurp.CFromSub() = Ctmp_sur;
    }
    else
    {
        Info << "zeroGradient\n" << endl;
    }
}


// Hydraulic head field at the subsurface top patch
if (flowSub and !flowSur and !flow)
{
    single = true;
    Info << "Flow at the subsurface domain" << endl;

    // Hydraulic head field at the interface patch (subsurface top patch)
    word interhSub
    (
        transportProperties.getOrDefault<word>("interhSub","latestTime")
    );

    Info << "Hydraulic head field at the interface patch (subsurface top patch): ";

    if (interhSub == "consth")
    {
        // Constant hydraulic head at the subsurface top patch
        scalar hSubTop
        (
            transportProperties.getOrDefault<scalar>("hSubTop",Zero)
        );

        hSub.boundaryFieldRef()[subSurfaceTopPatchID] == hSubTop;

        Info << "homogeneous\n" << endl;
    }
    else
    {
        Info << "the one resulting from the surface-subsurface simulation\n" << endl;
    }
}

// Temperature field at the subsurface top patch
if (heatTransferSub and !heatTransferSur and !heatTransfer)
{
    single = true;
    Info << "Heat transfer at the subsurface domain" << endl;

    // Temperature field at the interface patch (subsurface top patch)
    word interTSub
    (
        transportProperties.getOrDefault<word>("interTSub","zeroGradient")
    );

    Info << "Temperature field at the interface patch (subsurface top patch): ";

    if (interTSub == "latestTime" or interTSub == "constT")
    {
        scalarField Ttmp_sub(TSub.boundaryField()[subSurfaceTopPatchID].size(),Zero);

        if (interTSub == "latestTime")
        {
            Ttmp_sub = TSub.boundaryField()[subSurfaceTopPatchID];
            Info << "the one resulting from the surface-subsurface simulation\n" << endl;
        }
        else
        {
            // Constant temperature at the subsurface top patch
            scalar TSubTop
            (
                transportProperties.getOrDefault<scalar>("TSubTop",Zero)
            );

            Ttmp_sub = TSubTop;
            Info << "homogeneous\n" << endl;
        }

        subInterfaceHeatFvPatchScalarField& tempTSubp =
            dynamic_cast<subInterfaceHeatFvPatchScalarField&>(TSub.boundaryFieldRef()[subSurfaceTopPatchID]);

        tempTSubp.TFromSur() = Ttmp_sub;
    }
    else
    {
        Info << "zeroGradient\n" << endl;
    }
}

// Solute field at the subsurface top patch
if (soluteTransportSub and !soluteTransportSur and !soluteTransport)
{
    single = true;
    Info << "Solute transport at the subsurface domain" << endl;

    // Solute field at the interface patch (subsurface top patch)
    word interCSub
    (
        transportProperties.getOrDefault<word>("interCSub","zeroGradient")
    );

    Info << "Solute field at the interface patch (subsurface top patch): ";

    if (interCSub == "latestTime" or interCSub == "constC")
    {
        scalarField Ctmp_sub(CSub.boundaryField()[subSurfaceTopPatchID].size(),Zero);

        if (interCSub == "latestTime")
        {
            Ctmp_sub = CSub.boundaryField()[subSurfaceTopPatchID];
            Info << "the one resulting from the surface-subsurface simulation\n" << endl;
        }
        else
        {
            // Constant solute at the subsurface top patch
            scalar CSubTop
            (
                transportProperties.getOrDefault<scalar>("CSubTop",Zero)
            );

            Ctmp_sub = CSubTop;
            Info << "homogeneous\n" << endl;
        }

        subInterfaceSoluteFvPatchScalarField& tempCSubp =
            dynamic_cast<subInterfaceSoluteFvPatchScalarField&>(CSub.boundaryFieldRef()[subSurfaceTopPatchID]);

        tempCSubp.CFromSur() = Ctmp_sub;
    }
    else
    {
        Info << "zeroGradient\n" << endl;
    }
}

if (single)
{
    Info << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" << endl;
}
