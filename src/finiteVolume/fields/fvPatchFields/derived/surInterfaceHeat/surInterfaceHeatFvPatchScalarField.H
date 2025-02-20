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
    Foam::surInterfaceHeatFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    A scalar BC condition for the bottom patch of the surface flow domain.

    This BC relies on the flux direction on the bottom patch:

        (1) if flux direction subsurface -> surface (inlet)
            fixedValue: the value comes from the subsurface (TFromSub_)

        (2) if flux direction subsurface <- surface (outlet)
            zeroGradient: just let the scalar go out

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            surInterfaceHeat;
        phi             phi;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::mixedFvPatchField

SourceFiles
    surInterfaceHeatFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef surInterfaceHeatFvPatchScalarField_H
#define surInterfaceHeatFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class surInterfaceHeatFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class surInterfaceHeatFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- TFromSub: interpolated value of the same field from subsurface domain
        scalarField TFromSub_;

        //- Name of flux field
        word phiName_;

public:

    //- Runtime scalar information
    TypeName("surInterfaceHeat");


    // Constructors

        //- Construct from patch and internal field
        surInterfaceHeatFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surInterfaceHeatFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  surInterfaceHeatFvPatchScalarField onto a new patch
        surInterfaceHeatFvPatchScalarField
        (
            const surInterfaceHeatFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        surInterfaceHeatFvPatchScalarField
        (
            const surInterfaceHeatFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new surInterfaceHeatFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        surInterfaceHeatFvPatchScalarField
        (
            const surInterfaceHeatFvPatchScalarField&,
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
                new surInterfaceHeatFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        virtual scalarField& TFromSub()
        {
            return TFromSub_;
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
