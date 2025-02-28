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

#include "myTKEBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

/*Foam::scalar Foam::myTKEBCFvPatchScalarField::t() const
{
    return db().time().timeOutputValue(); remove this one, we donot need the time. in such cases we need since we need to access the time value from given database
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myTKEBCFvPatchScalarField::
myTKEBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    a_(0.0), // constructor based on defined value. sometimes user donot provide the value for a given the variable in the private data
    alpha_(0.0),
    Iu_(0.0),
    y_(1,0,0),
    n_(0,1,0)
    
{
}


Foam::myTKEBCFvPatchScalarField::
myTKEBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict // based on the dictionary this code ask for the value
)
:
    fixedValueFvPatchScalarField(p, iF),
    a_(readScalar(dict.lookup("MeanVelocity"))),
    alpha_(readScalar(dict.lookup("TerrRough"))),
    Iu_(readScalar(dict.lookup("TurbIntensity"))),
    y_(pTraits<vector>(dict.lookup("y"))),
    n_(pTraits<vector>(dict.lookup("normal")))
    //Umax_(dict.lookup<scalar>("Umax")),
    //Rm_(dict.lookup<scalar>("Rm")),

{


    fixedValueFvPatchScalarField::evaluate();

    /*
    // Initialise with the value entry if evaluation is not possible
    fvPatchVectorField::operator=
    (
        vectorField("value", dict, p.size())
    );
    */
}


Foam::myTKEBCFvPatchScalarField::
myTKEBCFvPatchScalarField
(
    const myTKEBCFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    a_(ptf.a_),
    alpha_(ptf.alpha_),
    Iu_(ptf.Iu_),
    y_(ptf.y_),
    n_(ptf.n_)
    //Umax_(ptf.Umax_),
    //Rm_(ptf.Rm_),
   
{}


Foam::myTKEBCFvPatchScalarField::
myTKEBCFvPatchScalarField
(
    const myTKEBCFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    a_(ptf.a_),
    alpha_(ptf.alpha_),
    Iu_(ptf.Iu_),
    y_(ptf.y_),
    n_(ptf.n_)
{}


Foam::myTKEBCFvPatchScalarField::
myTKEBCFvPatchScalarField
(
    const myTKEBCFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    a_(ptf.a_),
    alpha_(ptf.alpha_),
    Iu_(ptf.Iu_),
    y_(ptf.y_),
    n_(ptf.n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myTKEBCFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    //m(fieldData_, fieldData_); // not using the field data
}


void Foam::myTKEBCFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const myTKEBCFvPatchScalarField& tiptf =
        refCast<const myTKEBCFvPatchScalarField>(ptf);

    //fieldData_.rmap(tiptf.fieldData_, addr);
}


void Foam::myTKEBCFvPatchScalarField::updateCoeffs() //this is the function youre going to do some modification
{
    if (updated())
    {
        return;
    }
	// U = Umax (1- (r/R)^2), Umax is user defined
	
	const vectorField& c = patch().Cf();// going to give me a vector field y(patch().Cf()) is the face center of the every given face
	
	scalarField coord = c & y_;
	
	scalarField Uinf =  a_*pow(coord, alpha_); 
	scalarField TKE =  1.5*sqr(Uinf*Iu_);
	
	// y_ going to convert the vector field into scalar field by dot product with r and max-min
	//scalarField temp1 = 2*(r & y_)/ ((BB.max() - BB.min()) & y_);  //Calculating the local 1D coordinate
	
	scalarField::operator== (
		TKE 
	);
	
	

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::myTKEBCFvPatchScalarField::write //if you dont write, after 100s if you restarting it doesnt start
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "MeanVelocity", a_);
    writeEntry(os, "TerrRough", alpha_);
    writeEntry(os, "TurbIntensity", Iu_);
    writeEntry(os, "y", y_);
    writeEntry(os, "normal", n_);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        myTKEBCFvPatchScalarField
    );
}

// ************************************************************************* //
