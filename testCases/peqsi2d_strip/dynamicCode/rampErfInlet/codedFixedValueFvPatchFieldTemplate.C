/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "codedFixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "read.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ac1b37fca016c2cfe83880ec0e6fe064d68f1e5f
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void rampErfInlet_ac1b37fca016c2cfe83880ec0e6fe064d68f1e5f(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    rampErfInletFixedValueFvPatchVectorField
);


const char* const rampErfInletFixedValueFvPatchVectorField::SHA1sum =
    "ac1b37fca016c2cfe83880ec0e6fe064d68f1e5f";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rampErfInletFixedValueFvPatchVectorField::
rampErfInletFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct rampErfInlet sha1: ac1b37fca016c2cfe83880ec0e6fe064d68f1e5f"
            " from patch/dictionary\n";
    }
}


rampErfInletFixedValueFvPatchVectorField::
rampErfInletFixedValueFvPatchVectorField
(
    const rampErfInletFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct rampErfInlet sha1: ac1b37fca016c2cfe83880ec0e6fe064d68f1e5f"
            " from patch/DimensionedField/mapper\n";
    }
}


rampErfInletFixedValueFvPatchVectorField::
rampErfInletFixedValueFvPatchVectorField
(
    const rampErfInletFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct rampErfInlet sha1: ac1b37fca016c2cfe83880ec0e6fe064d68f1e5f "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rampErfInletFixedValueFvPatchVectorField::
~rampErfInletFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy rampErfInlet sha1: ac1b37fca016c2cfe83880ec0e6fe064d68f1e5f\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void rampErfInletFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs rampErfInlet sha1: ac1b37fca016c2cfe83880ec0e6fe064d68f1e5f\n";
    }

//{{{ begin code
    #line 12 "/home/sunkim/openfoam/RGP-13/test/peqsi2d_strip/0/U/inlet"

            const scalar t = this->db().time().value();
            const scalar ramp = min(t/50.0e-6, 1.0);
            const scalar r = 1.1e-3;
            const scalar sig = 3.0*8.2e-6;
            const scalar Uj = 20.0;
            vectorField& f = *this;
            const vectorField& Cf = patch().Cf();
            forAll(f, i)
            {
                const scalar y = Cf[i].y();
                const scalar prof = 0.5*(1.0 - erf((mag(y) - r)/sig));
                f[i] = vector(Uj*prof*ramp, 0, 0);
            }
        
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

