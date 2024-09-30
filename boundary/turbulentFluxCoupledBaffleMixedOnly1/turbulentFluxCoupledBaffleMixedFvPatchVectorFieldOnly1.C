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

\*---------------------------------------------------------------------------*/

#include "turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1.H"
#include "zeroGradientFvPatchField.H"
#include "calculatedFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

namespace Foam
{
namespace compressible
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1::
turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    TnbrName_("undefined-Tnbr"),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0)
{
    this->refValue() = vector(0.0, 0.0, 0.0);
    this->refGrad() = vector(0.0, 0.0, 0.0);
    this->valueFraction() = 1.0;
}


turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1::
turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    TnbrName_(dict.lookup("Tnbr")),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0.0)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    if (dict.found("thicknessLayers"))
    {
        dict.lookup("thicknessLayers") >> thicknessLayers_;
        dict.lookup("kappaLayers") >> kappaLayers_;

        if (thicknessLayers_.size() > 0)
        {
            // Calculate effective thermal resistance by harmonic averaging
            forAll(thicknessLayers_, iLayer)
            {
                contactRes_ += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
            }
            contactRes_ = 1.0/contactRes_;
        }
    }

    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = vectorField("refValue", dict, p.size());
        refGrad() = vectorField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = vector(0.0, 0.0, 0.0);
        valueFraction() = 1.0;
    }
}


turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1::
turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1
(
    const turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    contactRes_(ptf.contactRes_)
{}


turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1::
turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1
(
    const turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1& wtcsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(wtcsf, iF),
    TnbrName_(wtcsf.TnbrName_),
    thicknessLayers_(wtcsf.thicknessLayers_),
    kappaLayers_(wtcsf.kappaLayers_),
    contactRes_(wtcsf.contactRes_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    // Calculate the temperature by harmonic averaging
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //typedef turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1 thisType;
    typedef zeroGradientFvPatchField<vector> thisType;
    typedef calculatedFvPatchField<scalar> rhoType;
    
    const fvPatchVectorField& nbrTp =
        nbrPatch.lookupPatchField<volVectorField, vector>(TnbrName_);
        
    const fvPatchScalarField& nbrRhop =
        nbrPatch.lookupPatchField<volScalarField, scalar>("rho");
        
    const fvPatchScalarField& myRhop =
        patch().lookupPatchField<volScalarField, scalar>("rho");
/*
    if (!isA<thisType>(nbrTp))
    {
        FatalErrorInFunction
            << "Patch field for " << internalField().name() << " on "
            << patch().name() << " is of type " << thisType::typeName
            << endl << "The neighbouring patch field " << TnbrName_ << " on "
            << nbrPatch.name() << " is required to be the same, but is "
            << "currently of type " << nbrTp.type() << exit(FatalError);
    }
*/
    const thisType& nbrField = refCast<const thisType>(nbrTp);
    const rhoType& nbrRhoF = refCast<const rhoType>(nbrRhop);
    const rhoType& myRhoF = refCast<const rhoType>(myRhop);

    // Swap to obtain full local values of neighbour internal field
    tmp<vectorField> nbrIntFld(new vectorField(nbrField.size(), vector(0.0, 0.0, 0.0)));
    tmp<scalarField> nbrRho(new scalarField(nbrRhoF.size(), 0.0));
    tmp<scalarField> myRho(new scalarField(myRhoF.size(), 0.0));
    //tmp<vectorField> nbrKDelta(new vectorField(nbrField.size(), vector(0.0, 0.0, 0.0)));

    if (contactRes_ == 0.0)
    {
        nbrIntFld.ref() = nbrField.patchInternalField();
        nbrRho.ref()=nbrRhoF.patchInternalField();
        myRho.ref()=myRhoF.patchInternalField();
        //nbrKDelta.ref() = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();
    }
    else
    {
        nbrIntFld.ref() = nbrField;
        nbrRho.ref()=nbrRhoF;
        myRho.ref()=myRhoF;
        //nbrKDelta.ref() = contactRes_;
    }
    Info << "debug1 " <<endl;
    mpp.distribute(nbrIntFld.ref());
    mpp.distribute(nbrRho.ref());
    //mpp.distribute(myRho.ref());
    //mpp.distribute(nbrKDelta.ref());

    //tmp<vectorField> myKDelta = kappa(*this)*patch().deltaCoeffs();


    // Both sides agree on
    // - temperature : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
    // - gradient    : (temperature-fld)*delta
    // We've got a degree of freedom in how to implement this in a mixed bc.
    // (what gradient, what fixedValue and mixing coefficient)
    // Two reasonable choices:
    // 1. specify above temperature on one side (preferentially the high side)
    //    and above gradient on the other. So this will switch between pure
    //    fixedvalue and pure fixedgradient
    // 2. specify gradient and temperature such that the equations are the
    //    same on both sides. This leads to the choice of
    //    - refGradient = zero gradient
    //    - refValue = neighbour value
    //    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    //Info << refValue() <<endl;
    //Info << nbrIntFld() <<endl;
    //this->refValue() = nbrIntFld()*nbrRho.patchInternalField()/myRho.patchInternalField();
    //Info << "myRho(): " << myRho() <<endl;
    //Info << "nbrRho(): " << nbrRho() <<endl;
    this->refValue() = nbrIntFld()*nbrRho()/(myRho());
    this->refGrad() = vector(0.0, 0.0, 0.0);
    //this->valueFraction() = nbrKDelta()/(nbrKDelta() + myKDelta());
    this->valueFraction() = 1.0;//0.5;
    mixedFvPatchVectorField::updateCoeffs();

    // Restore tag
    UPstream::msgType() = oldTag;
}


void turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1::write
(
    Ostream& os
) const
{
    mixedFvPatchVectorField::write(os);
    writeEntry(os, "Tnbr", TnbrName_);
    writeEntry(os, "thicknessLayers", thicknessLayers_);
    writeEntry(os, "kappaLayers", kappaLayers_);

    //velocityCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    turbulentFluxCoupledBaffleMixedFvPatchVectorFieldOnly1
);
}
}
// ************************************************************************* //
