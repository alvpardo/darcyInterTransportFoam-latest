/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    restoreZRefMesh

Description
    Optionally (1) transforms the mesh points in the polyMesh directory according
    to the original location previously modified by changeZRefMesh, (2) restores
    the hRef value to the original value provided and, (3) updates the hydrological
    outputs to the new mesh location. Applicable to one (surface) or two (surface
    and subsurface) domains.

Usage
    Example usage:
    \verbatim
    restoreZRefMesh -numDomains 2 -time '20' -transformZmesh -waterHeightOutlet -updateHydroFields
    \endverbatim

Author
    Alvaro Pardo-Alvarez, modified from transformPoints.C, 23.08.2022

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "numDomains",
        "scalar",
        "Number of mesh domains. Either 1 or 2 (default)"
    );

    argList::addOption
    (
        "time",
        "time",
        "Specify the time to search from and apply the transformation"
        " (default is latest)"
    );

    argList::addBoolOption
    (
        "transformZmesh",
        "bool",
        "Transforms the Z mesh coordinates according to zRef"
    );

    argList::addBoolOption
    (
        "waterHeightOutlet",
        "bool",
        "Restores the original water height at the outlet of the domain"
    );

    argList::addBoolOption
    (
        "updateHydroFields",
        "bool",
        "Updates the hydraulic heads and piezometric levels "
        "to the original mesh elevation"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    scalar numDomains = args.getOrDefault<scalar>("numDomains",2);
    Info << "Number of mesh domains: " << numDomains << "\n" << endl;

    fileName meshDir;
    if (numDomains == 1)
    {
        meshDir = fvMesh::defaultRegion;
    }
    else if (numDomains == 2)
    {
        meshDir = "surface";
    }
    else
    {
        FatalError
        << "Error. This utility is only applicable to one (surface) or two (surface and subsurface)"
        " domains." << nl
        << exit(FatalError);
    }

    if (args.found("time"))
    {
        if (args["time"] == "constant")
        {
            runTime.setTime(instant(0, "constant"), 0);
        }
        else
        {
            const scalar timeValue = args.get<scalar>("time");
            runTime.setTime(instant(timeValue), 0);
        }
    }

    // Mesh
    fvMesh mesh
    (
       IOobject
       (
           meshDir,
           runTime.timeName(),
           runTime,
           IOobject::MUST_READ
       )
    );

    if (args.found("waterHeightOutlet"))
    {
        uniformDimensionedScalarField hRef
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        );

        label patchID = mesh.boundaryMesh().findPatchID("outlet");

        dimensionedScalar minOutlet
        (
            "minOutlet",dimLength,min(mesh.boundaryMesh()[patchID].localPoints()).component(vector::Z)
        );

        Info << "Minimum outlet elevation in the surface domain: " << minOutlet.value() << " [m]" << endl;

        hRef -= minOutlet;

        Info << "Restored hRef: " << hRef.value() << " [m]\n" << endl;

        hRef.write();
    }

    // Vertical mesh translation (L)
    uniformDimensionedScalarField meshShiftZ
    (
        IOobject
        (
            "meshShiftZ",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (args.found("transformZmesh"))
    {
        scalar mesh_shiftZ = meshShiftZ.value();
        Info << "Vertical mesh translation: " << -1 * mesh_shiftZ << " [m]\n" << endl;

        if (numDomains == 1)
        {
            // Surface mesh points
            pointIOField meshPoints
            (
                IOobject
                (
                    "points",
                    runTime.findInstance(polyMesh::meshSubDir,"points"),
                    polyMesh::meshSubDir,
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

            meshPoints -= mesh_shiftZ * vector(0,0,1);
            meshPoints.write();
        }
        else if (numDomains == 2)
        {
            // Surface mesh points
            pointIOField surfacePoints
            (
                IOobject
                (
                    "points",
                    runTime.findInstance("surface"/polyMesh::meshSubDir,"points"),
                    "surface"/polyMesh::meshSubDir,
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

            // Subsurface mesh points
            pointIOField subSurfacePoints
            (
                IOobject
                (
                    "points",
                    runTime.findInstance("subsurface"/polyMesh::meshSubDir,"points"),
                    "subsurface"/polyMesh::meshSubDir,
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

            surfacePoints -= mesh_shiftZ * vector(0,0,1);
            surfacePoints.write();

            subSurfacePoints -= mesh_shiftZ * vector(0,0,1);
            subSurfacePoints.write();
        }
    }

    if (args.found("updateHydroFields"))
    {
        // Piezometric level (L)
        volScalarField piezoSur
        (
            IOobject
            (
                "piezo",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        piezoSur == piezoSur - meshShiftZ;

        // Hydraulic head (L)
        volScalarField hSur
        (
            IOobject
            (
                "h",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        hSur == hSur - meshShiftZ;

        // Reading transportProperties dictionary
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Whether to clip the air fraction
        bool alphaAirClip
        (
            transportProperties.getOrDefault<Switch>("alphaAirClip",true)
        );

        // Water fraction threshold
        scalar alphaWaterThres
        (
            transportProperties.getOrDefault<scalar>("alphaWaterThres",0.5)
        );

        if (alphaAirClip)
        {
            // Water fraction (-)
            volScalarField alpha1
            (
                IOobject
                (
                    "alpha.water",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            );

            piezoSur *= pos0(alpha1 - alphaWaterThres);
            hSur *= pos0(alpha1 - alphaWaterThres);
        }

        Info << "Piezometric level updated." << endl;
        piezoSur.write();

        Info << "Surface hydraulic head updated.\n" << endl;
        hSur.write();

        if (numDomains == 2)
        {
            // Subsurface mesh
            fvMesh meshSubsurface
            (
               IOobject
               (
                   "subsurface",
                   runTime.timeName(),
                   runTime,
                   IOobject::MUST_READ
               )
            );

            volScalarField hSub
            (
                IOobject
                (
                    "h",
                    runTime.timeName(),
                    meshSubsurface,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                meshSubsurface
            );

            hSub == hSub - meshShiftZ;

            Info << "Subsurface hydraulic head updated.\n" << endl;
            hSub.write();
        }
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
