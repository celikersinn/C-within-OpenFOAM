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

#include "InletPulsatileBCFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
/*Foam::scalar Foam::InletPulsatileBCFvPatchVectorField::t() const
{
    return db().time().timeOutputValue(); //remove this one, we donot need the time. in such cases we need since we need to access the time value from given database
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::InletPulsatileBCFvPatchVectorField::
InletPulsatileBCFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umean_(0.0), // constructor based on defined value. sometimes user donot provide the value for a given the variable in the private data
    period_(0.0)
    
{
}


Foam::InletPulsatileBCFvPatchVectorField::
InletPulsatileBCFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict // based on the dictionary this code ask for the value
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umean_(readScalar(dict.lookup("Umean"))),
    period_(readScalar(dict.lookup("period")))

{


    fixedValueFvPatchVectorField::evaluate();

    /*
    // Initialise with the value entry if evaluation is not possible
    fvPatchVectorField::operator=
    (
        vectorField("value", dict, p.size())
    );
    */
}


Foam::InletPulsatileBCFvPatchVectorField::
InletPulsatileBCFvPatchVectorField
(
    const InletPulsatileBCFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Umean_(ptf.Umean_),
    period_(ptf.period_)
   
{}


Foam::InletPulsatileBCFvPatchVectorField::
InletPulsatileBCFvPatchVectorField
(
    const InletPulsatileBCFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    Umean_(ptf.Umean_),
    period_(ptf.period_)
{}


Foam::InletPulsatileBCFvPatchVectorField::
InletPulsatileBCFvPatchVectorField
(
    const InletPulsatileBCFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    Umean_(ptf.Umean_),
    period_(ptf.period_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::InletPulsatileBCFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    //m(fieldData_, fieldData_); // not using the field data
}


void Foam::InletPulsatileBCFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const InletPulsatileBCFvPatchVectorField& tiptf =
        refCast<const InletPulsatileBCFvPatchVectorField>(ptf);

    //fieldData_.rmap(tiptf.fieldData_, addr);
}


void Foam::InletPulsatileBCFvPatchVectorField::updateCoeffs() //this is the function youre going to do some modification
{
    if (updated())
    {
        return;
    }
	scalar curTime_ = this->db().time().timeOutputValue();

    scalar Umean_ = 0.511;
    scalar period_ = 0.928;
    scalar pi = 3.141592654;

    scalar UTime1 = ((-((Umean_) + (-0.0080467)*cos(1 *2.*pi*curTime_/period_) +
                  (0.086438)*sin(1 *2.*pi*curTime_/period_) +
                  (-0.0057042)*cos(2 *2.*pi*curTime_/period_) +
                  (0.036684)*sin(2 *2.*pi*curTime_/period_) + 
                  (-0.037873)*cos(3 *2.*pi*curTime_/period_) + 
                  (0.023462)*sin(3 *2.*pi*curTime_/period_) +
                  (-0.019703)*cos(4 *2.*pi*curTime_/period_) +
                  (0.0023512)*sin(4 *2.*pi*curTime_/period_) +
                  (-0.019301)*cos(5 *2.*pi*curTime_/period_) +
                  (0.006191)*sin(5 *2.*pi*curTime_/period_) +
                  (-0.015524)*cos(6 *2.*pi*curTime_/period_) +
                  (-0.0082778)*sin(6 *2.*pi*curTime_/period_) +
                  (-0.013376)*cos(7 *2.*pi*curTime_/period_) +
                  (-0.007184)*sin(7 *2.*pi*curTime_/period_) +
                  (8.3728e-005)*cos(8 *2.*pi*curTime_/period_) +
                  (-0.0034172)*sin(8 *2.*pi*curTime_/period_) +
                  (-0.0061556)*cos(9 *2.*pi*curTime_/period_) +
                  (-0.0056501)*sin(9 *2.*pi*curTime_/period_) +
                  (0.0014531)*cos(10*2.*pi*curTime_/period_) +
                  (-0.0035759)*sin(10*2.*pi*curTime_/period_))
           )/ Umean_ * 0.205);
	
	tmp<vectorField> n = patch().nf(); // surface normal to the that particular patch

    vectorField::operator==(n*UTime1);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::InletPulsatileBCFvPatchVectorField::write //if you dont write, after 100s if you restarting it doesnt start
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "Umean", Umean_);
    writeEntry(os, "period", period_);

}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        InletPulsatileBCFvPatchVectorField
    );
}

// ************************************************************************* //
