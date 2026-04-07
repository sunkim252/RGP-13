/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::rfSpecie

Description
    Base class of the real-fluid thermophysical property types.
    Extends Foam::specie with critical properties for real gas models.

SourceFiles
    rfSpecieI.H
    rfSpecie.C

\*---------------------------------------------------------------------------*/

#ifndef rfSpecie_H
#define rfSpecie_H

#include "word.H"
#include "scalar.H"
#include "dictionary.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class rfSpecie;

inline rfSpecie operator+(const rfSpecie&, const rfSpecie&);
inline rfSpecie operator*(const scalar, const rfSpecie&);
inline rfSpecie operator==(const rfSpecie&, const rfSpecie&);

Ostream& operator<<(Ostream&, const rfSpecie&);


/*---------------------------------------------------------------------------*\
                           Class rfSpecie Declaration
\*---------------------------------------------------------------------------*/

class rfSpecie
{
    // Private Data

        //- Name of specie
        word name_;

        //- Mass fraction of this specie in mixture
        scalar Y_;

        //- Molecular weight of specie [kg/kmol]
        scalar molWeight_;

        // Critical data for real gas models
        scalar Tc_;      // Critical temperature [K]
        scalar Pc_;      // Critical pressure    [Pa]
        scalar Vc_;      // Critical volume      [m^3/kmol]
        scalar omega_;   // Acentric factor      [-]
        scalar miui_;    // Dipole moment        [-]
        scalar kappai_;  // Association factor   [-]
        scalar sigmvi_;  // Diffusion volume     [-]


public:

    //- Runtime type information
    ClassName("rfSpecie");


    // Constructors

        //- Construct from components without name
        inline rfSpecie
        (
            const scalar Y,
            const scalar molWeight,
            const scalar Tc,
            const scalar Pc,
            const scalar Vc,
            const scalar omega,
            const scalar miui,
            const scalar kappai,
            const scalar sigmvi
        );

        //- Construct from components with name
        inline rfSpecie
        (
            const word& name,
            const scalar Y,
            const scalar molWeight,
            const scalar Tc,
            const scalar Pc,
            const scalar Vc,
            const scalar omega,
            const scalar miui,
            const scalar kappai,
            const scalar sigmvi
        );

        //- Construct as named copy
        inline rfSpecie(const word& name, const rfSpecie&);

        //- Construct from name and dictionary
        rfSpecie(const word& name, const dictionary& dict);

        //- Copy constructor
        rfSpecie(const rfSpecie&) = default;


    // Member Functions

        // Access

            //- Name
            inline const word& name() const;

            //- Molecular weight [kg/kmol]
            inline scalar W() const;

            //- Mass fraction of this specie in mixture
            inline scalar Y() const;

            //- Gas constant [J/kg/K]
            inline scalar R() const;

            //- Return critical temperature [K]
            inline scalar Tc() const;

            //- Return critical pressure [Pa]
            inline scalar Pc() const;

            //- Return critical volume [m^3/kmol]
            inline scalar Vc() const;

            //- Return acentric factor [-]
            inline scalar omega() const;

            //- Return dipole moment [-]
            inline scalar miui() const;

            //- Return association factor [-]
            inline scalar kappai() const;

            //- Return diffusion volume [-]
            inline scalar sigmvi() const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator=(const rfSpecie&);
        inline void operator+=(const rfSpecie&);
        inline void operator*=(const scalar);


    // Friend operators

        inline friend rfSpecie operator+(const rfSpecie&, const rfSpecie&);
        inline friend rfSpecie operator*(const scalar, const rfSpecie&);
        inline friend rfSpecie operator==(const rfSpecie&, const rfSpecie&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const rfSpecie&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rfSpecieI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
