/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "darcyGradHeadFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::darcyGradHeadFvPatchScalarField::darcyGradHeadFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    qName_("q"),
    KName_("K")
{}


Foam::darcyGradHeadFvPatchScalarField::darcyGradHeadFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    qName_(dict.getOrDefault<word>("q", "q")),
    KName_(dict.getOrDefault<word>("K", "K"))
{
    patchType() = dict.getOrDefault<word>("patchType", word::null);
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::darcyGradHeadFvPatchScalarField::darcyGradHeadFvPatchScalarField
(
    const darcyGradHeadFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    qName_(ptf.qName_),
    KName_(ptf.KName_)
{}


Foam::darcyGradHeadFvPatchScalarField::darcyGradHeadFvPatchScalarField
(
    const darcyGradHeadFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    qName_(wbppsf.qName_),
    KName_(wbppsf.KName_)
{}


Foam::darcyGradHeadFvPatchScalarField::darcyGradHeadFvPatchScalarField
(
    const darcyGradHeadFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    qName_(wbppsf.qName_),
    KName_(wbppsf.KName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::darcyGradHeadFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchVectorField& q =
        patch().lookupPatchField<volVectorField,vector>(qName_);

    const fvPatchTensorField& K =
        patch().lookupPatchField<volTensorField,tensor>(KName_);

    gradient() = -inv(K) & q & patch().nf();
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::darcyGradHeadFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    if (qName_ != "q")
    {
        os.writeKeyword("q") << qName_ << token::END_STATEMENT << nl;
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        darcyGradHeadFvPatchScalarField
    );
}


// ************************************************************************* //
