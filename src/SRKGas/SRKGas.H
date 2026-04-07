/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
    Foam::SRKGas

Description
    Soave-Redlich-Kwong cubic equation of state for gases.

    Validated against NIST data.

Usage
    \table
        Property     | Description
        Tc           | Critical temperature [K]
        Pc           | Critical pressure [Pa]
        omega        | Acentric factor [-]
    \endtable

    The EOS parameters (b, coef1, coef2, coef3) are computed from Tc, Pc
    and omega in the dictionary constructor.

SourceFiles
    SRKGasI.H
    SRKGas.C

\*---------------------------------------------------------------------------*/

#ifndef SRKGas_H
#define SRKGas_H

#include "autoPtr.H"
#include "rfSpecie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class SRKGas;

template<class Specie>
inline SRKGas<Specie> operator+
(
    const SRKGas<Specie>&,
    const SRKGas<Specie>&
);

template<class Specie>
inline SRKGas<Specie> operator*
(
    const scalar,
    const SRKGas<Specie>&
);

template<class Specie>
inline SRKGas<Specie> operator==
(
    const SRKGas<Specie>&,
    const SRKGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const SRKGas<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class SRKGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class SRKGas
:
    public Specie
{
    // Private Data

        //- SRK co-volume parameter [m^3/kmol]
        scalar b_;

        //- Decomposed aAlpha coefficients
        scalar coef1_;
        scalar coef2_;
        scalar coef3_;


    // Private Member Functions

        //- Return Xterm in dHdT equation
        inline scalar Xterm(const scalar p, const scalar T) const;

        //- Return derivative of Z w.r.t. T
        inline scalar dZdT(const scalar p, const scalar T) const;

        //- Return derivative of Z w.r.t. T with pre-computed Z value
        inline scalar dZdT
        (
            const scalar p,
            const scalar T,
            const scalar Zv
        ) const;


public:

    // Constructors

        //- Construct from components
        inline SRKGas
        (
            const Specie& sp,
            const scalar& b,
            const scalar& coef1,
            const scalar& coef2,
            const scalar& coef3
        );

        //- Construct from name and dictionary
        SRKGas(const word& name, const dictionary& dict);

        //- Construct as named copy
        inline SRKGas(const word& name, const SRKGas&);

        //- Construct and return a clone
        inline autoPtr<SRKGas> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "SRKGas<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy contribution [J/kg]
            inline scalar h(const scalar p, const scalar T) const;

            //- Return Cp contribution [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy contribution [J/kg]
            inline scalar e(const scalar p, const scalar T) const;

            //- Return Cv contribution [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy contribution to the integral of Cp/T [J/kg/K]
            inline scalar sp(const scalar p, const scalar T) const;

            //- Return entropy contribution to the integral of Cv/T [J/kg/K]
            inline scalar sv(const scalar p, const scalar T) const;

            //- Return compressibility [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;

            //- Return volumetric coefficient of thermal expansion [1/T]
            inline scalar alphav(const scalar p, const scalar T) const;

            //- Update coefficients for mixture
            void updateEoS
            (
                const scalar& b,
                const scalar& coef1,
                const scalar& coef2,
                const scalar& coef3
            )
            {
                b_ = b;
                coef1_ = coef1;
                coef2_ = coef2;
                coef3_ = coef3;
            }


    // IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member Operators

        inline void operator+=(const SRKGas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend SRKGas operator+ <Specie>
        (
            const SRKGas&,
            const SRKGas&
        );

        friend SRKGas operator* <Specie>
        (
            const scalar s,
            const SRKGas&
        );

        friend SRKGas operator== <Specie>
        (
            const SRKGas&,
            const SRKGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const SRKGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SRKGasI.H"

#ifdef NoRepository
    #include "SRKGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
