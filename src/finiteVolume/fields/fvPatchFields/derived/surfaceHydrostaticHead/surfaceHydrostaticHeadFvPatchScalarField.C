/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "surfaceHydrostaticHeadFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceHydrostaticHeadFvPatchScalarField
::surfaceHydrostaticHeadFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::surfaceHydrostaticHeadFvPatchScalarField
::surfaceHydrostaticHeadFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(patchInternalField());
    }
}


Foam::surfaceHydrostaticHeadFvPatchScalarField
::surfaceHydrostaticHeadFvPatchScalarField
(
    const surfaceHydrostaticHeadFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::surfaceHydrostaticHeadFvPatchScalarField
::surfaceHydrostaticHeadFvPatchScalarField
(
    const surfaceHydrostaticHeadFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf)
{}


Foam::surfaceHydrostaticHeadFvPatchScalarField
::surfaceHydrostaticHeadFvPatchScalarField
(
    const surfaceHydrostaticHeadFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceHydrostaticHeadFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    /*fvMesh mesh
    (
       IOobject
       (
           surDomainName_,
           db().time().timeName(),
           db().time(),
           IOobject::MUST_READ
       )
    );

    IOdictionary inOutStatic_hSur
    (
        IOobject
        (
            "inOutStatic_hSur",
            db().time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalar hPatch
    (
        inOutStatic_hSur.get<scalar>(surPatchName_)
    );*/

    scalarField hPatch = *this;

    forAll (patch(),i)
    {
        if (patch().name() == "inlet")
        {
            hPatch[i] = hInFromSur_;
        }

        if (patch().name() == "outlet")
        {
            hPatch[i] = hOutFromSur_;
        }

        if (patch().name() == "aquifer")
        {
            hPatch[i] = hRiverFromSur_[i];
        }
    }

    operator==(hPatch);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::surfaceHydrostaticHeadFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        surfaceHydrostaticHeadFvPatchScalarField
    );
}

// ************************************************************************* //
