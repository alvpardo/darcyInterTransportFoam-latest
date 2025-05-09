/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "VoFmeanVelocityForce.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "uniformDimensionedFields.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(VoFmeanVelocityForce, 0);
    addToRunTimeSelectionTable(option, VoFmeanVelocityForce, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::VoFmeanVelocityForce::writeProps
(
    const scalar gradP
) const
{
    // Only write on output time
    if (mesh_.time().writeTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        );
        propsDict.add("gradient",gradP);
        propsDict.regIOobject::write();
    }

    // write every time step, so that other solver can read gradP
    if(gradPWrite_)
    {
        uniformDimensionedVectorField deltaTgradP
        (
            IOobject
            (
                "deltaTgradP",
                mesh_.time().constant(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            dimensionedVector("deltaTgradP",dimensionSet(1,-2,-2,0,0,0,0),flowDir_*gradP)
        );

        deltaTgradP.write();
    };
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::VoFmeanVelocityForce::VoFmeanVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    Ubar_(coeffs_.get<vector>("Ubar")),
    gradP0_(0.0),
    dGradP_(0.0),
    flowDir_(Ubar_/mag(Ubar_)),
    relaxation_(coeffs_.getOrDefault<scalar>("relaxation", 1)),
    rAPtr_(nullptr),
    gradPWrite_(coeffs_.getOrDefault<bool>("gradPWrite", false))
{
    coeffs_.readEntry("fields", fieldNames_);

    if (fieldNames_.size() != 1)
    {
        FatalErrorIn
        (
            "Foam::fv::VoFmeanVelocityForce::"
            "VoFmeanVelocityForce"
            "("
                "const word&, "
                "const word&, "
                "const dictionary&, "
                "const fvMesh&"
            ")"
        )   << "Source can only be applied to a single field.  Current "
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    fv::option::resetApplied();

    // Read the initial pressure gradient from file if it exists
    IFstream propsFile
    (
        mesh_.time().timePath()/"uniform"/(name_ + "Properties")
    );

    if (propsFile.good())
    {
        Info<< "    Reading pressure gradient from file" << endl;
        dictionary propsDict(propsFile);
        propsDict.readEntry("gradient", gradP0_);
    }

    Info<< "    Initial pressure gradient = " << gradP0_ << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fv::VoFmeanVelocityForce::magUbarAve
(
    const volVectorField& U
) const
{
    // Read IB cell information

    volScalarField& alpha = const_cast<volScalarField&>
        (mesh_.lookupObject<volScalarField>("alpha.water"));

    scalar Vtotal = 0.0;

    scalar magUbarAve = 0.0;

    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label celli = cells_[i];
        scalar volCell = cv[celli];
        Vtotal += volCell*alpha[celli]; 
        magUbarAve += (flowDir_ & U[celli])*volCell*alpha[celli];
    }

    reduce(magUbarAve, sumOp<scalar>());
    reduce(Vtotal, sumOp<scalar>());

    magUbarAve /= Vtotal;

    return magUbarAve;
}


void Foam::fv::VoFmeanVelocityForce::correct(volVectorField& U)
{
    // Read IB cell information

    volScalarField& alpha = const_cast<volScalarField&>
        (mesh_.lookupObject<volScalarField>("alpha.water"));

    scalar Vtotal = 0.0;

    const scalarField& rAU = rAPtr_().internalField();

    // Integrate flow variables over cell set
    scalar rAUave = 0.0;
    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label celli = cells_[i];
        scalar volCell = cv[celli];
        Vtotal += volCell*alpha[celli];
        rAUave += rAU[celli]*volCell*alpha[celli];
    }

    // Collect across all processors
    reduce(rAUave, sumOp<scalar>());
    reduce(Vtotal, sumOp<scalar>());

    // Volume averages
    rAUave /= Vtotal;

    scalar magUbarAve = this->magUbarAve(U);

    // Calculate the pressure gradient increment needed to adjust the average
    // flow-rate to the desired value
    dGradP_ = relaxation_*(mag(Ubar_) - magUbarAve)/rAUave;

    // Apply correction to velocity field
    forAll(cells_, i)
    {
        label celli = cells_[i];
        U[celli] += flowDir_*rAU[celli]*dGradP_;
    }

    U.correctBoundaryConditions();

    scalar gradP = gradP0_ + dGradP_;

    Info<< "Pressure gradient source: uncorrected Ubar = " << magUbarAve
        << ", pressure gradient = " << gradP << endl;

    writeProps(gradP);
}


void Foam::fv::VoFmeanVelocityForce::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    volScalarField& alpha = const_cast<volScalarField&>
        (mesh_.lookupObject<volScalarField>("alpha.water"));

    scalar gradP = gradP0_ + dGradP_;
    scalarField alpha_d (alpha);
    forAll(alpha_d, celli)
    {
        if(alpha[celli]>0.5)
        {
            alpha_d[celli] = 1.0;
        }
        else
        {
            alpha_d[celli] = 0;
        }
    }

    UIndirectList<vector>(Su, cells_) = flowDir_*gradP*alpha_d;

    eqn += Su;
}


void Foam::fv::VoFmeanVelocityForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    this->addSup(eqn, fieldi);
}


void Foam::fv::VoFmeanVelocityForce::constrain
(
    fvMatrix<vector>& eqn,
    const label
)
{
    if (!rAPtr_)
    {
        rAPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_ + ":rA",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                1.0/eqn.A()
            )
        );
    }
    else
    {
        rAPtr_() = 1.0/eqn.A();
    }

    gradP0_ += dGradP_;
    dGradP_ = 0.0;
}


bool Foam::fv::VoFmeanVelocityForce::read(const dictionary& dict)
{
    NotImplemented;

    return false;
}


// ************************************************************************* //
