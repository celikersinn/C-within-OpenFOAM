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

#include "myvelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

/*Foam::scalar Foam::myvelocityFvPatchVectorField::t() const
{
    return db().time().timeOutputValue(); remove this one, we donot need the time. in such cases we need since we need to access the time value from given database
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myvelocityFvPatchVectorField::
myvelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umax_(0.0), // constructor based on defined value. sometimes user donot provide the value for a given the variable in the private data
    y_(1,0,0),
    n_(0,1,0)
    
{
}


Foam::myvelocityFvPatchVectorField::
myvelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict // based on the dictionary this code ask for the value
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umax_(readScalar(dict.lookup("Umax"))),
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


Foam::myvelocityFvPatchVectorField::
myvelocityFvPatchVectorField
(
    const myvelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Umax_(ptf.Umax_),
    y_(ptf.y_),
    n_(ptf.n_)
    //Umax_(ptf.Umax_),
    //Rm_(ptf.Rm_),
   
{}


Foam::myvelocityFvPatchVectorField::
myvelocityFvPatchVectorField
(
    const myvelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    Umax_(ptf.Umax_),
    y_(ptf.y_),
    n_(ptf.n_)
{}


Foam::myvelocityFvPatchVectorField::
myvelocityFvPatchVectorField
(
    const myvelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    Umax_(ptf.Umax_),
    y_(ptf.y_),
    n_(ptf.n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myvelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    //m(fieldData_, fieldData_); // not using the field data
}


void Foam::myvelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const myvelocityFvPatchVectorField& tiptf =
        refCast<const myvelocityFvPatchVectorField>(ptf);

    //fieldData_.rmap(tiptf.fieldData_, addr);
}


void Foam::myvelocityFvPatchVectorField::updateCoeffs() //this is the function youre going to do some modification
{
    if (updated())
    {
        return;
    }
	// U = Umax (1- (r/R)^2), Umax is user defined
	//Bounding range => we will have a patch get to know the info from given patch  => from that we are gonna ask what are the faces this patch has  
	//=> ask for the points known as a localPoints in openfoam  => ask for the min and max value
	
	//first patch give me access to the given patch, second patch faces of a patch, localpoints access of a vertices of given faces
	boundBox BB (patch().patch().localPoints(),true); //then we can use the constructor boundbox since we providing a pointfield 
	
	// using the bounding box, we can calculate the corners of a given this patch
	
	vector midCord = 0.5* (BB.max()+ BB.min());// center of this given complete patch, giving the length of the radius
	
	vectorField r = patch().Cf() - midCord ;// going to give me a vector field y(patch().Cf()) is the face center of the every given face
	
	// y_ going to convert the vector field into scalar field by dot product with r and max-min
	scalarField temp1 = 2*(r & y_)/ ((BB.max() - BB.min()) & y_);  //Calculating the local 1D coordinate
	
	fixedValueFvPatchVectorField::operator== (
		Umax_ * n_* (1- sqr(temp1)) //we need normal since converting the vector field

    /* alternative way 
    const vectorField& c =patch().Cf();
	scalarField coord_r = (c & y_) //Calculaing the local 1D coordinate n terms only y coorditane(0 1 0)
	fixedValueFvPatchVectorField::operator== (
		Umax_ * n_* (1- sqr(coord_r/Rm_)) //we need normal since converting the vector field.n give you flow direction ie x-direction (1 0 0)
		// rm mean radius of specific pipe we specify
    */
	);
	
	

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::myvelocityFvPatchVectorField::write //if you dont write, after 100s if you restarting it doesnt start
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "Umax", Umax_);
    writeEntry(os, "y", y_);
    writeEntry(os, "normal", n_);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        myvelocityFvPatchVectorField
    );
}

// ************************************************************************* //
