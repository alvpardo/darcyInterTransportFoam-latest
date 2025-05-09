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

Class
    Foam::flowRateInletParabolicVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    prescribes a 2D or 3D parabolic velocity field normal to the inlet patch, 
    given a specified flow rate, the type of inlet geometry (open or closed), 
    an average velocity multiplier (to calculate the peak of the parabola/
    paraboloid) and the water depth direction. The resulting velocity profile 
    (either half or a quarter of parabola) will depend on the geometry type of 
    the inlet.

    For a mass-based flux:
    - the flow rate should be provided in kg/s
    - if \c rho is "none" the flow rate is in m3/s
    - otherwise \c rho should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must
      specify the inlet density using the \c rhoInlet entry

    For a volumetric-based flux:
    - the flow rate is in m3/s

Usage
    \table
        Property     | Description             | Required    | Default value
        Uprofile     | velocity profile        | yes         | parabolic
        geometryInlet| open or closed inlet patch geometry | yes | open
        massFlowRate | mass flow rate [kg/s]   | no          |
        volumetricFlowRate | volumetric flow rate [m3/s]| no |
        avgUmulti    | average velocity multiplier to calculate the maximum velocity | yes |
        depthDir     | water depth direction   | yes         |
        rho          | density field name      | no          | rho
        rhoInlet     | inlet density           | no          |
    \endtable

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    <patchName>
    {
        type                flowRateInletParabolicVelocity;
        Uprofile            parabolic;
        geometryInlet       open;
        volumetricFlowRate  0.2;
        avgUmulti           1.3;
        depthDir            (0 0 1);
        value               uniform (0 0 0);
    }
    \endverbatim

    Example of the boundary condition specification for a mass flow rate:
    \verbatim
    <patchName>
    {
        type                flowRateInletParabolicVelocity;
        Uprofile            parabolic;
        geometryInlet       open;
        massFlowRate        0.2;
        avgUmulti           1.3;
        depthDir            (0 0 1);
        rho                 rho;
        rhoInlet            1.0;
        value               uniform (0 0 0);
    }
    \endverbatim

    The \c flowRate entry is a \c Function1 of time, see Foam::Function1Types.

Note
    - \c rhoInlet is required for the case of a mass flow rate, where the
      density field is not available at start-up
    - The value is positive into the domain (as an inlet)
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types
    Foam::flowRateInletVelocityFvPatchVectorField

SourceFiles
    flowRateInletParabolicVelocityFvPatchVectorField.C

Author
    Alvaro Pardo-Alvarez, modified from flowRateInletVelocity BC, 16.12.2022

\*---------------------------------------------------------------------------*/

#ifndef flowRateInletParabolicVelocityFvPatchVectorField_H
#define flowRateInletParabolicVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class flowRateInletParabolicVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class flowRateInletParabolicVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet integral flow rate
        autoPtr<Function1<scalar>> flowRate_;

        //- Is volumetric?
        bool volumetric_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        //- Rho initialisation value (for start; if value not supplied)
        scalar rhoInlet_;

        //- Velocity profile shape
        word Uprofile_;

        //- Inlet patch geometry
        word geometryInlet_;

        //- Average velocity multiplier
        scalar avgUmulti_;

        //- Water depth direction
        vector depthDir_;


    // Private member functions

        //- Update the patch values given the appropriate density type and value
        template<class RhoType>
        void updateValues(const RhoType& rho);


public:

   //- Runtime type information
   TypeName("flowRateInletParabolicVelocity");


   // Constructors

        //- Construct from patch and internal field
        flowRateInletParabolicVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flowRateInletParabolicVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletParabolicVelocityFvPatchVectorField
        //  onto a new patch
        flowRateInletParabolicVelocityFvPatchVectorField
        (
            const flowRateInletParabolicVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        flowRateInletParabolicVelocityFvPatchVectorField
        (
            const flowRateInletParabolicVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new flowRateInletParabolicVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        flowRateInletParabolicVelocityFvPatchVectorField
        (
            const flowRateInletParabolicVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new flowRateInletParabolicVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
