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
    changeZRefMesh

Description
    Optionally (1) transforms the mesh points in the polyMesh directory according
    to the zRef value provided and, (2) updates the hRef value to the new surface
    mesh location. Applicable to one (surface) or two (surface and subsurface)
    domains.

Usage
    Example usage:
    \verbatim
    changeZRefMesh -numDomains 2 -zRef 0.0 -transformZmesh -waterHeightOutlet
    \endverbatim

Author
    Alvaro Pardo-Alvarez, modified from transformPoints.C, 23.08.2022

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
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
        "zRef",
        "scalar",
        "New mesh zRef. The default value is 0.0"
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
        "Sets water height at the outlet of the domain"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    scalar numDomains = args.getOrDefault<scalar>("numDomains",2);
    Info << "Number of mesh domains: " << numDomains << "\n" << endl;

    scalar zRef = args.getOrDefault<scalar>("zRef",Zero);
    Info << "New reference mesh elevation: " << zRef << " [m]\n" << endl;

    fileName meshDir;
    scalar mesh_shiftZ = Zero;
    if (numDomains == 1)
    {
        // Mesh points
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

        scalar mesh_minZ = min(meshPoints.component(vector::Z));
        //scalar mesh_minZ = meshSurface.bounds().min().component(vector::Z);
        Info << "Minimum surface mesh elevation: " << mesh_minZ << " [m]" << endl;

        mesh_shiftZ = zRef - mesh_minZ;
        Info << "Vertical mesh translation: " << mesh_shiftZ << " [m]\n" << endl;

        if (args.found("transformZmesh"))
        {
            meshPoints += mesh_shiftZ * vector(0,0,1);
            meshPoints.write();
        }

        meshDir = fvMesh::defaultRegion;
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

        scalar meshSur_minZ = min(surfacePoints.component(vector::Z));
        Info << "Minimum surface mesh elevation: " << meshSur_minZ << " [m]" << endl;

        scalar meshSub_minZ = min(subSurfacePoints.component(vector::Z));
        Info << "Minimum subsurface mesh elevation: " << meshSub_minZ << " [m]\n" << endl;

        scalar meshSur_shiftZ = zRef - meshSur_minZ;
        scalar meshSub_shiftZ = zRef - meshSub_minZ;

        scalar min_shiftZ = min(mag(meshSur_shiftZ),mag(meshSub_shiftZ));
        Info << "Minimum mesh shift: " << min_shiftZ << " [m]\n" << endl;

        if (min_shiftZ == mag(meshSur_shiftZ))
        {
            mesh_shiftZ = meshSur_shiftZ;
            Info << "Vertical mesh translation: " << mesh_shiftZ << " [m]\n" << endl;
        }
        else
        {
            mesh_shiftZ = meshSub_shiftZ;
            Info << "Vertical mesh translation: " << mesh_shiftZ << " [m]\n" << endl;
        }

        if (args.found("transformZmesh"))
        {
            surfacePoints += mesh_shiftZ * vector(0,0,1);
            surfacePoints.write();

            subSurfacePoints += mesh_shiftZ * vector(0,0,1);
            subSurfacePoints.write();
        }

        meshDir = "surface";
    }
    else
    {
        FatalError
        << "Error. This utility is only applicable to one (surface) or two (surface and subsurface)"
        " domains." << nl
        << exit(FatalError);
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

    // Vertical mesh translation (L)
    uniformDimensionedScalarField meshShiftZ
    (
        IOobject
        (
            "meshShiftZ",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar("meshShiftZ",dimLength,mesh_shiftZ)
    );

    meshShiftZ.write();

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

        hRef += minOutlet;

        Info << "Updated hRef: " << hRef.value() << " [m]\n" << endl;

        hRef.write();
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
