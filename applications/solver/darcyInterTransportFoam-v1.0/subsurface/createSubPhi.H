/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Global
    createSubPhi

Description
    Creates and initialises the relative face-flux field phi.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Darcy volumetric flux (L^3/T)
surfaceScalarField phi_q
(
    IOobject
    (
        "phi_q",
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    meshSub,
    dimensionedScalar(dimensionSet(0,3,-1,0,0,0,0),Zero)
);

// Porewater/average linear volumetric flux (L^3/T)
surfaceScalarField phi_USub
(
    IOobject
    (
        "phi_USub",
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    meshSub,
    dimensionedScalar(dimensionSet(0,3,-1,0,0,0,0),Zero)
);

// Solute mass flux (M/T)
surfaceScalarField phiCSub
(
    IOobject
    (
        "phiC",
        runTime.timeName(),
        meshSub,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    meshSub,
    dimensionedScalar(dimensionSet(1,0,-1,0,0,0,0),Zero)
);


// ************************************************************************* //
