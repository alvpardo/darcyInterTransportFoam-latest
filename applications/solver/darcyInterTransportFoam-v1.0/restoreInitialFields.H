if (reset)
{
    Info << "Reestablishment of the original..." << endl;
}

// Reestablishment of the original surface domain fields previous to their reset
if (resetSur)
{
    Info << "    ...surface domain fields..." << endl;

    rmDir(mesh.time().path()/initialTime/"surface");
    cp(mesh.time().path()/"." + initialTime/"surface", mesh.time().path()/initialTime/"surface");
}

// Reestablishment of the original subsurface domain fields previous to their reset
if (resetSub)
{
    Info << "    ...subsurface domain fields..." << endl;

    rmDir(meshSub.time().path()/initialTime/"subsurface");
    cp(meshSub.time().path()/"." + initialTime/"subsurface", meshSub.time().path()/initialTime/"subsurface");
}


// Reestablishment of the original surface heat transfer fields previous to their reset
if (resetHeatSur)
{
    Info << "    ...surface heat transfer fields..." << endl;

    cp(mesh.time().path()/"." + initialTime/"surface"/"T", mesh.time().path()/initialTime/"surface"/"T");
}

// Reestablishment of the original subsurface heat transfer fields previous to their reset
if (resetHeatSub)
{
    Info << "    ...subsurface heat transfer fields..." << endl;

    cp(meshSub.time().path()/"." + initialTime/"subsurface"/"T", meshSub.time().path()/initialTime/"subsurface"/"T");
}


// Reestablishment of the original surface solute transport fields previous to their reset
if (resetSoluteSur)
{
    Info << "    ...surface solute transport fields..." << endl;

    cp(mesh.time().path()/"." + initialTime/"surface"/soluteNameSur, mesh.time().path()/initialTime/"surface"/soluteNameSur);
    cp(mesh.time().path()/"." + initialTime/"surface"/"D_" + soluteNameSur, mesh.time().path()/initialTime/"surface"/"D_" + soluteNameSur);

    // phiC file is removed because it is not always present in the time directories
    rm(mesh.time().path()/initialTime/"surface"/"phiC");
    cp(mesh.time().path()/"." + initialTime/"surface"/"phiC", mesh.time().path()/initialTime/"surface"/"phiC");
}

// Reestablishment of the original subsurface solute transport fields previous to their reset
if (resetSoluteSub)
{
    Info << "    ...subsurface solute transport fields..." << endl;

    cp(meshSub.time().path()/"." + initialTime/"subsurface"/soluteNameSub, meshSub.time().path()/initialTime/"subsurface"/soluteNameSub);
    cp(meshSub.time().path()/"." + initialTime/"subsurface"/"D_" + soluteNameSub, meshSub.time().path()/initialTime/"subsurface"/"D_" + soluteNameSub);

    // phiC file is removed because it is not always present in the time directories
    rm(meshSub.time().path()/initialTime/"subsurface"/"phiC");
    cp(meshSub.time().path()/"." + initialTime/"subsurface"/"phiC", meshSub.time().path()/initialTime/"subsurface"/"phiC");
}

// Removal of the initial time temporary copy
if (reset)
{
    Info << "...previous to their reset" << endl;
    Info << "\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" << endl;
    rmDir(mesh.time().path()/"." + initialTime);
}
