/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1

Description
    This boundary condition provides a turbulentFluxCoupledBaffleMixed condition,
    calculated as:

        \f[
            Q = Q_{0} + Q_{p} + s*Q_{t}
        \f]

    where
    \vartable
        s      | single scalar value [units]
        Q_{0}  | single vector value [units]
        Q_{p}  | vector field across patch [units]
        Q_{t}  | vector function of time [units]
    \endtable

Usage
    \table
        Property   | Description                  | Req'd? | Default
        scalarData | single scalar value          | yes    |
        data       | single vector value          | yes    |
        fieldData  | vector field across patch    | yes    |
        timeVsData | vector function of time      | yes    |
        wordData   | word, eg name of data object | no     | wordDefault
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       turbulentFluxCoupledBaffleMixed;
        scalarData -1;
        data       (1 0 0);
        fieldData  uniform (3 0 0);
        timeVsData table (
                             (0 (0 0 0))
                             (1 (2 0 0))
                          );
        wordName   anotherName;
        value      uniform (4 0 0); // optional initial value
    }
    \endverbatim

SourceFiles
    turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1_H
#define turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1_H

//#include "fixedValueFvPatchFields.H"
//#include "Function1.H"

#include "mixedFvPatchFields.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
/*---------------------------------------------------------------------------*\
                          Class turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1 Declaration
\*---------------------------------------------------------------------------*/

class turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1
:
    //public fixedValueFvPatchVectorField
    public mixedFvPatchVectorField
    //public velocityCoupledBase
{
    // Private Data

        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;

        //- Total contact resistance
        scalar contactRes_;


public:

    //- Runtime type information
    TypeName("turbulentFluxCoupledBaffleMixedOnly1");


    // Constructors

        //- Construct from patch and internal field
        turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1
        (
            const turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );
        /*
        //- Disallow copy without setting internal field reference
        turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1
        (
            const turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1&
        ) = delete;
	*/
        //- Copy constructor setting internal field reference
        turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1
        (
            const turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1&,
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
                new turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions
	/*
        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchVectorField&, const labelList&);


        // Evaluation functions
	*/
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
