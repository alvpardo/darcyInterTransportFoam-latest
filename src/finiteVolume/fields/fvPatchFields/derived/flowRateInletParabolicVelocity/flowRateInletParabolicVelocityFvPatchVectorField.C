/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "flowRateInletParabolicVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateInletParabolicVelocityFvPatchVectorField::
flowRateInletParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(),
    volumetric_(false),
    rhoName_("rho"),
    rhoInlet_(0.0),
    Uprofile_("parabolic"),
    geometryInlet_("open"),
    avgUmulti_(1.0),
    depthDir_(0,0,1)
{}


Foam::flowRateInletParabolicVelocityFvPatchVectorField::
flowRateInletParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    rhoInlet_(dict.getOrDefault<scalar>("rhoInlet", -VGREAT)),
    Uprofile_(dict.getOrDefault<word>("Uprofile","parabolic")),
    geometryInlet_(dict.getOrDefault<word>("geometryInlet","open")),
    avgUmulti_(readScalar(dict.lookup("avgUmulti"))),
    depthDir_(dict.get<vector>("depthDir"))
{
    if (dict.found("volumetricFlowRate"))
    {
        volumetric_ = true;
        flowRate_ = Function1<scalar>::New("volumetricFlowRate", dict);
        rhoName_ = "rho";
    }
    else if (dict.found("massFlowRate"))
    {
        volumetric_ = false;
        flowRate_ = Function1<scalar>::New("massFlowRate", dict);
        rhoName_ = dict.getOrDefault<word>("rho", "rho");
    }
    else
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Please supply either 'volumetricFlowRate' or"
            << " 'massFlowRate' and 'rho'" << exit(FatalIOError);
    }

    // Value field require if mass based
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::flowRateInletParabolicVelocityFvPatchVectorField::
flowRateInletParabolicVelocityFvPatchVectorField
(
    const flowRateInletParabolicVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_.clone()),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    Uprofile_(ptf.Uprofile_),
    geometryInlet_(ptf.geometryInlet_),
    avgUmulti_(ptf.avgUmulti_),
    depthDir_(ptf.depthDir_)
{}


Foam::flowRateInletParabolicVelocityFvPatchVectorField::
flowRateInletParabolicVelocityFvPatchVectorField
(
    const flowRateInletParabolicVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_.clone()),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    Uprofile_(ptf.Uprofile_),
    geometryInlet_(ptf.geometryInlet_),
    avgUmulti_(ptf.avgUmulti_),
    depthDir_(ptf.depthDir_)
{}


Foam::flowRateInletParabolicVelocityFvPatchVectorField::
flowRateInletParabolicVelocityFvPatchVectorField
(
    const flowRateInletParabolicVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_.clone()),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    Uprofile_(ptf.Uprofile_),
    geometryInlet_(ptf.geometryInlet_),
    avgUmulti_(ptf.avgUmulti_),
    depthDir_(ptf.depthDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class RhoType>
void Foam::flowRateInletParabolicVelocityFvPatchVectorField::updateValues
(
    const RhoType& rho
)
{
    const scalar t = db().time().timeOutputValue();

    const vectorField n(patch().nf());

    // Get range and orientation
    boundBox bb(patch().patch().localPoints(), true);
    //Info << "Bounding box: " << bb << endl;

    const vectorField& faceCtr = patch().Cf();
    //Info << "Patch faces centres: " << faceCtr << endl;

    // Patch normal vector: weighted mean of patch face normal vectors
    vector n_vector = gSum(patch().magSf()*n)/gSum(patch().magSf());
    //Info << "Patch mean normal: " << n_vector << endl;

    vector transDir = n_vector ^ depthDir_;
    //Info << "Transversal direction: " << transDir << endl;

    scalar pHeight = ((bb.max() - bb.min()) & depthDir_); // Patch height
    scalar pWidth = ((bb.max() - bb.min()) & transDir); // Patch width

    scalarField verCoord(faceCtr.size(),Zero);
    scalarField horCoord(faceCtr.size(),Zero);

    if (geometryInlet_ == "open")
    {
        vector paraCtr = vector(0.5*(bb.max() + bb.min()).x(),0.5*(bb.max() + bb.min()).y(),bb.max().z()); // Top patch centre
        //Info << "Parabola centre: " << paraCtr << endl;

        // Calculate local 1-D coordinate for the parabolic profile
        verCoord = ((faceCtr - paraCtr) & depthDir_)/pHeight;
        horCoord = ((faceCtr - paraCtr) & transDir)/(0.5*pWidth);	// Multiply by 0.5 to get the distance from the centre of the parabola to the patch end
    }
    else if (geometryInlet_ == "closed")
    {
        vector paraCtr = vector(0.5*(bb.max() + bb.min()).x(),0.5*(bb.max() + bb.min()).y(),0.5*(bb.max() + bb.min()).z()); // Patch centre
        //Info << "Parabola centre: " << paraCtr << endl;

        // Calculate local 1-D coordinate for the parabolic profile
        verCoord = ((faceCtr - paraCtr) & depthDir_)/(0.5*pHeight);
        horCoord = ((faceCtr - paraCtr) & transDir)/(0.5*pWidth);	// Multiply by 0.5 to get the distance from the centre of the parabola to the patch end
    }
    else
    {
        FatalErrorInFunction
            << "Please supply either 'geometryInlet_ == open' or"
            << " 'geometryInlet_ == closed'" << exit(FatalError);
    }

    const scalar avgU = -flowRate_->value(t)/gSum(rho*patch().magSf());

    if (Uprofile_ == "parabolic")
    {
        // Parabolic U profile
        operator==(n*avgUmulti_*avgU*(1.0 - sqr(verCoord)));
    }
    else if (Uprofile_ == "paraboloidal")
    {
        // Parabolidal U profile
        operator==(n*avgUmulti_*avgU*(1.0 - sqr(verCoord))*(1.0 - sqr(horCoord)));
        //operator==(n*((16*avgUmulti_*avgU)/(sqr(pHeight*2)*sqr(pWidth)))*(pHeight - ((faceCtr - paraCtr) & depthDir_))*(pHeight + ((faceCtr - paraCtr) & depthDir_))
                 //*(pWidth/2 - ((faceCtr - paraCtr) & transDir))*(pWidth/2 + ((faceCtr - paraCtr) & transDir)));
    }
    else
    {
        FatalErrorInFunction
            << "Please supply either 'Uprofile=parabolic' or"
            << " 'Uprofile=paraboloidal'" << exit(FatalError);
    }
}


void Foam::flowRateInletParabolicVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (volumetric_ || rhoName_ == "none")
    {
        updateValues(one());
    }
    else
    {
        // Mass flow-rate
        if (db().foundObject<volScalarField>(rhoName_))
        {
            const fvPatchField<scalar>& rhop =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);

            updateValues(rhop);
        }
        else
        {
            // Use constant density
            if (rhoInlet_ < 0)
            {
                FatalErrorInFunction
                    << "Did not find registered density field " << rhoName_
                    << " and no constant density 'rhoInlet' specified"
                    << exit(FatalError);
            }

            updateValues(rhoInlet_);
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::flowRateInletParabolicVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    flowRate_->writeData(os);
    if (!volumetric_)
    {
        os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
        os.writeEntryIfDifferent<scalar>("rhoInlet", -VGREAT, rhoInlet_);
    }
    os.writeEntry<word>("Uprofile",Uprofile_);
    os.writeEntry<word>("geometryInlet",geometryInlet_);
    os.writeEntry<scalar>("avgUmulti",avgUmulti_);
    os.writeEntry<vector>("depthDir",depthDir_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       flowRateInletParabolicVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
