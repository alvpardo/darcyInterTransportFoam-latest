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
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    checkMeshesMatch

Description
    Checks whether the number of faces and their coordinates match at the 
    interface patch of the surface and subsurface meshes.

Usage
    Example usage:
    \verbatim
    checkMeshesMatch -interfacePatch "riverbed" -distThres 5e-9
    \endverbatim

Author
    Alvaro Pardo-Alvarez, modified from transformPoints.C, 23.09.2022

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "interfacePatch",
        "word",
        "Patch through which the domain descomposition is transfered. "
    );

    argList::addOption
    (
        "distThres",
        "scalar",
        "Upper threshold applied to the 3D distance between the face centres of the surface and "
        "subsurface interface patches. Couples of face centres from both patches whose distance "
        "is equal or smaller than the defined threshold would match. The default value is 0.0"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    word interfacePatch = args.get<word>("interfacePatch");
    Info << "Patch through which the domain descomposition is transfered: " << interfacePatch << "\n" << endl;

    scalar distThres = args.getOrDefault<scalar>("distThres",Zero);
    Info << "Upper threshold applied to the 3D distance between the face centres of both meshes "
        "at the interface patch: " << distThres << "\n" << endl;

    // Surface mesh
    fvMesh mesh
    (
       IOobject
       (
           "surface",
           runTime.timeName(),
           runTime,
           IOobject::MUST_READ
       )
    );

    // Subsurface mesh
    fvMesh meshSub
    (
       IOobject
       (
           "subsurface",
           runTime.timeName(),
           runTime,
           IOobject::MUST_READ
       )
    );

    // Bottom patch of the surface mesh
    label surfaceBottomPatchID = mesh.boundaryMesh().findPatchID(interfacePatch);

    // Top patch of the subsurface mesh
    label subSurfaceTopPatchID = meshSub.boundaryMesh().findPatchID(interfacePatch);

    // Face centre nodes
    vectorField pts_inter_mesh = mesh.C().boundaryField()[surfaceBottomPatchID];
    vectorField pts_inter_meshSub = meshSub.C().boundaryField()[subSurfaceTopPatchID];

    // Face vertices
    //vectorField pts_mesh = mesh.boundaryMesh()[surfaceBottomPatchID].localPoints();
    //vectorField pts_meshSub = meshSub.boundaryMesh()[subSurfaceTopPatchID].localPoints();

    Info << "Matching check of the surface and subsurface meshes at the interface patch:\n" << endl;

    label count = 0;
    while (count < 2)
    {
        if (Pstream::parRun())
        {
            // Parallel mode
            if (Pstream::master())
            {
                if (count == 0)
                {
                    // Master node
                    // Check whether the number of faces of both meshes is the same at the interface patch
                    Info << "    The case is run in parallel mode\n"
                         << "    Both surface and subsurface meshes are splitted into " << Pstream::nProcs() << " subdomains\n" << endl;

                    Info << "    Number of faces of both meshes at the interface patch of each couple of surface-subsurface subdomains:" << endl;
                    Info << "    Processor #" << Pstream::masterNo() << ": surface (" << pts_inter_mesh.size() << ") - subsurface (" << pts_inter_meshSub.size() << ")" << endl;

                    if (pts_inter_mesh.size() != pts_inter_meshSub.size())
                    {
                        FatalError
                        << "The number of faces of both meshes differs at the interface patch of the surface-subsurface subdomain couple #" << Pstream::masterNo()
                        << ": surface (" << pts_inter_mesh.size() << ") - subsurface (" << pts_inter_meshSub.size() << ")\n"
                        << "    This is a necessary condition to achieve a correct output\n"
                        << "    Check blockMeshDict for potential errors and try again" << nl
                        << exit(FatalError);
                    }
                }
                else
                {
                    // Check whether the face coordinates of both meshes match at the interface patch
                    forAll(pts_inter_mesh,i)
                    {
                        label minPointi = 0;
                        scalar minDistSqr = magSqr(pts_inter_meshSub[minPointi] - pts_inter_mesh[i]);

                        for (label j=1;j<pts_inter_meshSub.size();j++)
                        {
                            scalar distSqr = magSqr(pts_inter_meshSub[j] - pts_inter_mesh[i]);

                            if (distSqr < minDistSqr)
                            {
                                minDistSqr = distSqr;
                                minPointi = j;
                            }
                        }

                        if (mag(pts_inter_meshSub[minPointi]-pts_inter_mesh[i]) > distThres)
                        {
                            FatalError
                            << "The face coordinates of both meshes do not match at the interface patch\n"
                            << "    This is a necessary condition to achieve a correct output\n"
                            << "    Check blockMeshDict for potential errors and try again" << nl
                            << exit(FatalError);
                        }
                    }
                }

                // Slave nodes
                for (const int slave : Pstream::subProcs())
                {
                    // Master receiving from slave
                    IPstream fromSlave(Pstream::commsTypes::blocking, slave);
                    fromSlave >> pts_inter_mesh;
                    fromSlave >> pts_inter_meshSub;

                    if (count == 0)
                    {
                        // Check whether the number of faces of both meshes is the same at the interface patch
                        Info << "    Processor #" << slave << ": surface (" << pts_inter_mesh.size() << ") - subsurface (" << pts_inter_meshSub.size() << ")" << endl;

                        if (pts_inter_mesh.size() != pts_inter_meshSub.size())
                        {
                            FatalError
                            << "The number of faces of both meshes differs at the interface patch of the surface-subsurface subdomain couple #" << slave
                            << ": surface (" << pts_inter_mesh.size() << ") - subsurface (" << pts_inter_meshSub.size() << ")\n"
                            << "    This is a necessary condition to achieve a correct output\n"
                            << "    Check blockMeshDict for potential errors and try again" << nl
                            << exit(FatalError);
                        }
                    }
                    else
                    {
                        // Check whether the face coordinates of both meshes match at the interface patch
                        forAll(pts_inter_mesh,i)
                        {
                            label minPointi = 0;
                            scalar minDistSqr = magSqr(pts_inter_meshSub[minPointi] - pts_inter_mesh[i]);

                            for (label j=1;j<pts_inter_meshSub.size();j++)
                            {
                                scalar distSqr = magSqr(pts_inter_meshSub[j] - pts_inter_mesh[i]);

                                if (distSqr < minDistSqr)
                                {
                                    minDistSqr = distSqr;
                                    minPointi = j;
                                }
                            }

                            if (mag(pts_inter_meshSub[minPointi]-pts_inter_mesh[i]) > distThres)
                            {
                                FatalError
                                << "The face coordinates of both meshes do not match at the interface patch\n"
                                << "    This is a necessary condition to achieve a correct output\n"
                                << "    Check blockMeshDict for potential errors and try again" << nl
                                << exit(FatalError);
                            }
                        }
                    }
                }
            }
            else
            {
                // Slave sending to master
                OPstream toMaster(Pstream::commsTypes::blocking, Pstream::masterNo());
                toMaster << pts_inter_mesh;
                toMaster << pts_inter_meshSub;
            }

            if (count == 0)
            {
                Info << nl << "    The number of faces of both meshes is the same at the interface patch of each couple of surface-subsurface subdomains" << endl;
            }
            else
            {
                Info << "    The face coordinates of both meshes also match at the interface patch of each couple of surface-subsurface subdomains\n" << endl;
            }
        }
        else
        {
            // Serial mode
            if (count == 0)
            {
                // Check whether the number of faces of both meshes is the same at the interface patch
                Info << "    The case is run in serial mode\n" << endl;

                Info << "    Number of faces of both meshes at the interface patch:"
                     << " surface (" << pts_inter_mesh.size() << ") - subsurface (" << pts_inter_meshSub.size() << ")" << endl;

                if (pts_inter_mesh.size() != pts_inter_meshSub.size())
                {
                    FatalError
                    << "The number of faces of both meshes differs at the interface patch"
                    << ": surface (" << pts_inter_mesh.size() << ") - subsurface (" << pts_inter_meshSub.size() << ")\n"
                    << "    This is a necessary condition to achieve a correct output\n"
                    << "    Check blockMeshDict for potential errors and try again" << nl
                    << exit(FatalError);
                }
                else
                {
                    Info << nl << "    The number of faces of both meshes is the same at the interface patch" << endl;
                }
            }
            else
            {
                // Check whether the face coordinates of both meshes match at the interface patch
                forAll(pts_inter_mesh,i)
                {
                    label minPointi = 0;
                    scalar minDistSqr = magSqr(pts_inter_meshSub[minPointi] - pts_inter_mesh[i]);

                    for (label j=1;j<pts_inter_meshSub.size();j++)
                    {
                        scalar distSqr = magSqr(pts_inter_meshSub[j] - pts_inter_mesh[i]);

                        if (distSqr < minDistSqr)
                        {
                            minDistSqr = distSqr;
                            minPointi = j;
                        }
                    }

                    if (mag(pts_inter_meshSub[minPointi]-pts_inter_mesh[i]) > distThres)
                    {
                        FatalError
                        << "The face coordinates of both meshes do not match at the interface patch\n"
                        << "    This is a necessary condition to achieve a correct output\n"
                        << "    Check blockMeshDict for potential errors and try again" << nl
                        << exit(FatalError);
                    }
                }

                Info << "    The face coordinates of both meshes also match at the interface patch\n" << endl;
                }
        }

        count++;
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
