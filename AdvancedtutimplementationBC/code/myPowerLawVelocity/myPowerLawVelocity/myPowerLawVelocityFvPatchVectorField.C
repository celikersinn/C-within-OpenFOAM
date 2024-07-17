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

#include "myPowerLawVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

/*Foam::scalar Foam::myPowerLawVelocityFvPatchVectorField::t() const
{
    return db().time().timeOutputValue(); remove this one, we donot need the time. in such cases we need since we need to access the time value from given database
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myPowerLawVelocityFvPatchVectorField::
myPowerLawVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    a_(0.0), // constructor based on defined value. sometimes user donot provide the value for a given the variable in the private data
    alpha_(0.0),
    y_(1,0,0),
    n_(0,1,0)
    
{
}


Foam::myPowerLawVelocityFvPatchVectorField::
myPowerLawVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict // based on the dictionary this code ask for the value
)
:
    fixedValueFvPatchVectorField(p, iF),
    a_(readScalar(dict.lookup("MeanVelocity"))),
    alpha_(readScalar(dict.lookup("TerrRough"))),
    y_(pTraits<vector>(dict.lookup("y"))),
    n_(pTraits<vector>(dict.lookup("normal")))
    //Umax_(dict.lookup<scalar>("Umax")),
    //Rm_(dict.lookup<scalar>("Rm")),

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


Foam::myPowerLawVelocityFvPatchVectorField::
myPowerLawVelocityFvPatchVectorField
(
    const myPowerLawVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    a_(ptf.a_),
    alpha_(ptf.alpha_),
    y_(ptf.y_),
    n_(ptf.n_)
    //Umax_(ptf.Umax_),
    //Rm_(ptf.Rm_),
   
{}


Foam::myPowerLawVelocityFvPatchVectorField::
myPowerLawVelocityFvPatchVectorField
(
    const myPowerLawVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    a_(ptf.a_),
    alpha_(ptf.alpha_),
    y_(ptf.y_),
    n_(ptf.n_)
{}


Foam::myPowerLawVelocityFvPatchVectorField::
myPowerLawVelocityFvPatchVectorField
(
    const myPowerLawVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    a_(ptf.a_),
    alpha_(ptf.alpha_),
    y_(ptf.y_),
    n_(ptf.n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myPowerLawVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    //m(fieldData_, fieldData_); // not using the field data
}


void Foam::myPowerLawVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const myPowerLawVelocityFvPatchVectorField& tiptf =
        refCast<const myPowerLawVelocityFvPatchVectorField>(ptf);

    //fieldData_.rmap(tiptf.fieldData_, addr);
}


void Foam::myPowerLawVelocityFvPatchVectorField::updateCoeffs() //this is the function youre going to do some modification
{
    if (updated())
    {
        return;
    }
	// U = Umax (1- (r/R)^2), Umax is user defined
	
	//first patch give me access to the given patch, second patch faces of a patch, localpoints access of a vertices of given faces
	boundBox BB (patch().patch().localPoints(),true); //then we can use the constructor boundbox since we providing a pointfield 
	
	// using the bounding box, we can calculate the corners of a given this patch
	
	vector midCord = 0.5* (BB.max()+ BB.min());// center of this given complete patch, giving the length of the radius
	
	const vectorField& c = patch().Cf();// going to give me a vector field y(patch().Cf()) is the face center of the every given face
	
	scalarField coord = c & y_;
	
	vectorField Uinf = n_ * a_*pow(coord, alpha_); // give you the freestream direction velocity
	
	// y_ going to convert the vector field into scalar field by dot product with r and max-min
	//scalarField temp1 = 2*(r & y_)/ ((BB.max() - BB.min()) & y_);  //Calculating the local 1D coordinate
	
	fixedValueFvPatchVectorField::operator== (
		Uinf 
	);
	
	

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::myPowerLawVelocityFvPatchVectorField::write //if you dont write, after 100s if you restarting it doesnt start
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "MeanVelocity", a_);
    writeEntry(os, "TerrRough", alpha_);
    writeEntry(os, "y", y_);
    writeEntry(os, "normal", n_);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        myPowerLawVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
