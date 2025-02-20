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

Class
    Foam::subInterfaceHeatFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    A scalar BC condition for the top patch of the subsurface flow domain.

    This BC relies on the flux direction on subsurface top:

        (1) if flux direction subsurface -> surface (outlet)
            zeroGradient: just let the scalar go out

        (2) if flux direction subsurface <- surface (inlet)
            fixedValue: the value comes from the surface (TFromSur_)

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            subInterfaceHeat;
        phi             phi_USub;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::mixedFvPatchField

SourceFiles
    subInterfaceHeatFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef subInterfaceHeatFvPatchScalarField_H
#define subInterfaceHeatFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class subInterfaceHeatFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class subInterfaceHeatFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- TFromSur: interpolated value of the same field from surface domain
        scalarField TFromSur_;

        //- Name of flux field
        word phiName_;

public:

    //- Runtime scalar information
    TypeName("subInterfaceHeat");


    // Constructors

        //- Construct from patch and internal field
        subInterfaceHeatFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        subInterfaceHeatFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  subInterfaceHeatFvPatchScalarField onto a new patch
        subInterfaceHeatFvPatchScalarField
        (
            const subInterfaceHeatFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        subInterfaceHeatFvPatchScalarField
        (
            const subInterfaceHeatFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new subInterfaceHeatFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        subInterfaceHeatFvPatchScalarField
        (
            const subInterfaceHeatFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new subInterfaceHeatFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        virtual scalarField& TFromSur()
        {
            return TFromSur_;
        }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
