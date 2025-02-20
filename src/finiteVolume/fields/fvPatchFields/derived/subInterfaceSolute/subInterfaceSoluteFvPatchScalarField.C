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

#include "subInterfaceSoluteFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subInterfaceSoluteFvPatchScalarField
::subInterfaceSoluteFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_("phi_USub")
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::subInterfaceSoluteFvPatchScalarField
::subInterfaceSoluteFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_(dict.getOrDefault<word>("phi", "phi_USub"))
{
    refValue() = 0.0;

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

    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::subInterfaceSoluteFvPatchScalarField
::subInterfaceSoluteFvPatchScalarField
(
    const subInterfaceSoluteFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{}


Foam::subInterfaceSoluteFvPatchScalarField
::subInterfaceSoluteFvPatchScalarField
(
    const subInterfaceSoluteFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    phiName_(ptf.phiName_)
{}


Foam::subInterfaceSoluteFvPatchScalarField
::subInterfaceSoluteFvPatchScalarField
(
    const subInterfaceSoluteFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::subInterfaceSoluteFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //fluid flux on this patch
    const fvsPatchScalarField& phi_USubp =
        patch().lookupPatchField<surfaceScalarField,scalar>(phiName_);

    forAll(CFromSur_, facei)
    {
        if(neg(phi_USubp[facei])) //flux surface -> subsurface: fixedValue
        {
            refValue()[facei] = CFromSur_[facei];
            valueFraction()[facei] = 1.0;
        }
        else                      //flux subsurface -> surface: zeroGradient
        {
            refGrad()[facei] = 0.0;
            valueFraction()[facei] = 0.0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::subInterfaceSoluteFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    if (phiName_ != "phi_USub")
    {
        os.writeKeyword("phi_USub") << phiName_ << token::END_STATEMENT << nl;
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        subInterfaceSoluteFvPatchScalarField
    );
}

// ************************************************************************* //
