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
    transferMeshDecomp

Description
    Transfers the mesh decomposition of one mesh to another by means of the 
    interface patch shared by both meshes.

Usage
    Example usage:
    \verbatim
    transferMeshDecomp -nSubDomains 10 -interfacePatch "riverbed" -transferDir "surToSub" 
    -distThres 0.2 -rbDeltaX 10 -rbDeltaY 10 -rbDeltaZ 10 -overWriteCellDecomp 
    -overWriteCellDecompBase
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
        "nSubDomains",
        "scalar",
        "Number of mesh subdomains"
    );

    argList::addOption
    (
        "interfacePatch",
        "word",
        "Patch through which the domain descomposition is transfered. "
    );

    argList::addOption
    (
        "transferDir",
        "word",
        "Transfering direction of the domain descomposition. "
        "Either surToSub (default) or subToSur."
    );

    argList::addOption
    (
        "distThres",
        "scalar",
        "Upper threshold applied to the 2D distance between the cell centres of both meshes. "
        "The defined threshold sets the limit for the selection of the cells belonging each subdomain "
        "of the target mesh. The default value is 0.0"
    );

    argList::addOption
    (
        "rbDeltaX",
        "scalar",
        "Expansion of the X dimension of the rotated boxes for the full selection of the cells "
        "belonging to each subdomain of the target mesh. Its value should be normally small. "
        "The default value is 0.0"
    );

    argList::addOption
    (
        "rbDeltaY",
        "scalar",
        "Expansion of the Y dimension of the rotated boxes for the full selection of the cells "
        "belonging to each subdomain of the target mesh. Its value should be normally small. "
        "The default value is 0.0"
    );

    argList::addOption
    (
        "rbDeltaZ",
        "scalar",
        "Expansion of the Z dimension of the rotated boxes for the full selection of the cells "
        "belonging to each subdomain of the target mesh. Its value should correspond to the "
        "vertical extension of such mesh. The default value is 0.0"
    );

    argList::addBoolOption
    (
        "overWriteCellDecomp",
        "bool",
        "Overwrite constant/cellDecomposition file based on the cellDist file "
        "of the mesh to which the domain decomposition is transfered"
    );

    argList::addBoolOption
    (
        "overWriteCellDecompBase",
        "bool",
        "Overwrite constant/cellDecomposition file based on the cellDist file "
        "of the mesh from which the domain decomposition is transfered"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    if (args.options().empty())
    {
        FatalErrorInFunction
            << "Missing nSubDomains. Please provide a number."
            << exit(FatalError);
    }

    scalar nSubDomains = args.get<scalar>("nSubDomains");
    Info << "Number of mesh subdomains: " << nSubDomains << "\n" << endl;

    word interfacePatch = args.get<word>("interfacePatch");
    Info << "Patch through which the domain descomposition is transfered: " << interfacePatch << "\n" << endl;

    word transferDir = args.getOrDefault<word>("transferDir","surToSub");
    Info << "Transfering direction of the domain descomposition: " << transferDir << "\n" << endl;

    scalar distThres = args.getOrDefault<scalar>("distThres",Zero);
    Info << "Upper threshold applied to the 2D distance between the cell centres of both meshes: " << distThres << "\n" << endl;

    scalar rbDeltaX = args.getOrDefault<scalar>("rbDeltaX",Zero);
    Info << "Expansion of the X dimension of the rotated boxes for the full selection of the cells "
        "belonging to each subdomain of the target mesh: " << rbDeltaX << endl;

    scalar rbDeltaY = args.getOrDefault<scalar>("rbDeltaY",Zero);
    Info << "Expansion of the Y dimension of the rotated boxes for the full selection of the cells "
        "belonging to each subdomain of the target mesh: " << rbDeltaY << endl;

    scalar rbDeltaZ = args.getOrDefault<scalar>("rbDeltaZ",Zero);
    Info << "Expansion of the Z dimension of the rotated boxes for the full selection of the cells "
        "belonging to each subdomain of the target mesh: " << rbDeltaZ << "\n" << endl;

    bool overWriteCellDecomp = args.found("overWriteCellDecomp");
    bool overWriteCellDecompBase = args.found("overWriteCellDecompBase");

    fileName meshDirFrom;
    fileName meshDirTo;

    // Transfering direction of mesh decomposition: surface -> subsurface
    if (transferDir == "surToSub")
    {
        meshDirFrom = "surface";
        Info << "Overwrite constant/surface/cellDecomposition file based on the cellDist file "
                "of the mesh from which the domain decomposition is transfered: " << overWriteCellDecompBase << endl;

        meshDirTo = "subsurface";
        Info << "Overwrite constant/subsurface/cellDecomposition file based on the cellDist file "
                "of the mesh to which the domain decompostion is transfered: " << overWriteCellDecomp << "\n" << endl;
    }
    // Transfering direction of mesh decomposition: subsurface -> surface
    else if (transferDir == "subToSur")
    {
        meshDirFrom = "subsurface";
        Info << "Overwrite constant/subsurface/cellDecomposition file based on the cellDist file "
                "of the mesh from which the domain decomposition is transfered: " << overWriteCellDecompBase << endl;

        meshDirTo = "surface";
        Info << "Overwrite constant/surface/cellDecomposition file based on the cellDist file "
                "of the mesh to which the domain decompostion is transfered: " << overWriteCellDecomp << "\n" << endl;
    }

    // Mesh from which the domain decomposition is transfered
    fvMesh meshFrom
    (
       IOobject
       (
           meshDirFrom,
           runTime.timeName(),
           runTime,
           IOobject::MUST_READ
       )
    );

    // Mesh to which the domain decompostion is transfered
    fvMesh meshTo
    (
       IOobject
       (
           meshDirTo,
           runTime.timeName(),
           runTime,
           IOobject::MUST_READ
       )
    );

    // Interface patch of both meshes
    label interfacePatchID = meshFrom.boundaryMesh().findPatchID(interfacePatch);

    // Face centres of the interface patch
    fvPatchVectorField faceCtrs = meshFrom.C().boundaryField()[interfacePatchID];

    // Cell IDs of meshFrom
    volScalarField cellDistFrom
    (
        IOobject
        (
            "cellDist",
            runTime.timeName(),
            meshFrom,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        meshFrom
    );

    IOList<label> cellDecompositionFrom
    (
        IOobject
        (
            "cellDecomposition",
            runTime.constant(),
            meshFrom,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Cell IDs of meshTo
    volScalarField cellDistTo
    (
        IOobject
        (
            "cellDist",
            runTime.timeName(),
            meshTo,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        meshTo
    );

    IOList<label> cellDecompositionTo
    (
        IOobject
        (
            "cellDecomposition",
            runTime.constant(),
            meshTo,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    List<scalar> subdomainLimits(nSubDomains+1, Zero);
    subdomainLimits[0] = 0;
    scalar n=0;
    for (label nSub=0; nSub<nSubDomains; nSub++)
    {
        // Number of faces at each subdomain interface patch of meshFrom
        forAll (cellDistFrom.boundaryField()[interfacePatchID],i)
        {
            if (cellDistFrom.boundaryField()[interfacePatchID][i] == nSub)
            {
                n+=1;
            }
            subdomainLimits[nSub+1] = n;
        }

        // Face coordinates of each subdomain interface patch of meshFrom
        vectorField patchFrom (subdomainLimits[nSub+1]-subdomainLimits[nSub],vector::zero);
        scalar fc=0;
        forAll (cellDistFrom.boundaryField()[interfacePatchID],i)
        {
            if (cellDistFrom.boundaryField()[interfacePatchID][i] == nSub)
            {
                patchFrom[fc] = faceCtrs[i];
                fc+=1;
            }
        }

        // Vectors defining the rotated boxes which contain the cell centres of each meshTo subdomain
        vector origin = vector(min(patchFrom.component(vector::X))-rbDeltaX,min(patchFrom.component(vector::Y))-rbDeltaY,min(patchFrom.component(vector::Z))-rbDeltaZ);
        vector i_ = vector(max(patchFrom.component(vector::X))+rbDeltaX-origin.component(vector::X),min(patchFrom.component(vector::Y))-rbDeltaY-origin.component(vector::Y),0.0);
        vector j_ = vector(min(patchFrom.component(vector::X))-rbDeltaX-origin.component(vector::X),max(patchFrom.component(vector::Y))+rbDeltaY-origin.component(vector::Y),0.0);
        vector k_ = vector(0.0,0.0,max(patchFrom.component(vector::Z))+rbDeltaZ-origin.component(vector::Z));

        // Rotated box vertices
        pointField boxPoints(8);
        boxPoints[0] = origin;
        boxPoints[1] = origin + i_;
        boxPoints[2] = origin + i_ + j_;
        boxPoints[3] = origin + j_;
        boxPoints[4] = origin + k_;
        boxPoints[5] = origin + k_ + i_;
        boxPoints[6] = origin + k_ + i_ + j_;
        boxPoints[7] = origin + k_ + j_;

        labelList boxVerts(identity(8));

        const cellModel& hex = cellModel::ref(cellModel::HEX);

        // Get outwards pointing faces
        faceList boxFaces(cellShape(hex, boxVerts).faces());

        // Precalculate normals
        vectorField boxFaceNormals(boxFaces.size());
        forAll(boxFaces,i)
        {
            boxFaceNormals[i] = boxFaces[i].areaNormal(boxPoints);

            //Info << "Face:" << i << " position:" << boxFaces[i].centre(boxPoints)
            //    << " normal:" << boxFaceNormals[i] << endl;
        }

        // Set to zero the third dimension of patchFrom for the correct selection of meshTo cell centres
        forAll (patchFrom,i)
        {
            patchFrom[i].component(vector::Z) *= 0.0;
        }

        // Cell centres of the meshTo
        pointField ctrsTo = meshTo.cellCentres();

        // All the info for nearest
        //HashTable <scalar,point> nearest;

        // Check whether the meshTo cell centres are inside all the faces of the rotated boxes
        forAll(ctrsTo,celli)
        {
            bool inside = true;

            forAll(boxFaces,i)
            {
                const face& f = boxFaces[i];

                if (((ctrsTo[celli] - boxPoints[f[0]]) & boxFaceNormals[i]) > 0)
                {
                    inside = false;
                    break;
                }
            }

            if (inside)
            {
                // Set to zero the third dimension of ctrsTo for the correct selection of meshTo cell centres
                ctrsTo[celli].component(vector::Z) *= 0.0;

                // Selection of the meshTo cells corresponding to each extended meshFrom subdomain
                label minPointi = 0;
                scalar minDistSqr = magSqr(patchFrom[minPointi] - ctrsTo[celli]);

                for (label i=1;i<patchFrom.size();i++)
                {
                    scalar distSqr = magSqr(patchFrom[i] - ctrsTo[celli]);

                    if (distSqr < minDistSqr)
                    {
                        minDistSqr = distSqr;
                        minPointi = i;
                    }
                }

                //nearest.insert(patchFrom[minPointi],mag(patchFrom[minPointi]-ctrsTo[celli]));

                // Subdomain ID assignment to each selected meshTo cell
                if (mag(patchFrom[minPointi]-ctrsTo[celli]) <= distThres)
                {
                    cellDistTo[celli] = nSub;
                    cellDecompositionTo[celli] = nSub;
                }
            }
        }

        //Info << nearest << endl;
    }

    cellDistTo.write();

    // Overwrite cellDecomposition file of meshTo
    if (overWriteCellDecomp)
    {
        cellDecompositionTo.write();
    }

    // Overwrite cellDecomposition file of meshFrom
    if (overWriteCellDecompBase)
    {
        forAll (cellDistFrom,i)
        {
            cellDecompositionFrom[i] = cellDistFrom[i];
        }

        cellDecompositionFrom.write();
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
