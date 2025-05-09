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
    Foam::surfaceHydrostaticHeadFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a hydrostatic head distribution at the 
    inlet, outlet or bottom boundary of the subsurface domain based on the head 
    values obtained at the equivalent patches of the surface domain.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            surfaceHydrostaticHead;
        value           uniform 0; // placeholder
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    surfaceHydrostaticHeadFvPatchScalarField.C

Author
    Alvaro Pardo-Alvarez, 26.01.2022

\*---------------------------------------------------------------------------*/

#ifndef surfaceHydrostaticHeadFvPatchScalarField_H
#define surfaceHydrostaticHeadFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class surfaceHydrostaticHeadFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class surfaceHydrostaticHeadFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- hInFromSur: hydrostatic hydraulic head at the inlet patch of the surface domain
        scalar hInFromSur_;

        //- hOutFromSur: hydrostatic hydraulic head at the outlet patch of the surface domain
        scalar hOutFromSur_;

        //- hRiverFromSur: hydrostatic hydraulic head at the bottom patch of the surface domain
        scalarField hRiverFromSur_;

public:

    //- Runtime type information
    TypeName("surfaceHydrostaticHead");


    // Constructors

        //- Construct from patch and internal field
        surfaceHydrostaticHeadFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceHydrostaticHeadFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  surfaceHydrostaticHeadFvPatchScalarField onto a new patch
        surfaceHydrostaticHeadFvPatchScalarField
        (
            const surfaceHydrostaticHeadFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        surfaceHydrostaticHeadFvPatchScalarField
        (
            const surfaceHydrostaticHeadFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new surfaceHydrostaticHeadFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        surfaceHydrostaticHeadFvPatchScalarField
        (
            const surfaceHydrostaticHeadFvPatchScalarField&,
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
                new surfaceHydrostaticHeadFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        scalar& hInFromSur()
        {
            return hInFromSur_;
        }

        scalar& hOutFromSur()
        {
            return hOutFromSur_;
        }

        scalarField& hRiverFromSur()
        {
            return hRiverFromSur_;
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
