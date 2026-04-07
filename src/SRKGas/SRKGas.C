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

\*---------------------------------------------------------------------------*/

#include "SRKGas.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::SRKGas<Specie>::SRKGas
(
    const word& name,
    const dictionary& dict
)
:
    Specie(name, dict)
{
    const scalar Tc = this->Tc();
    const scalar Pc = this->Pc();
    const scalar omega = this->omega();

    b_ = 0.08664*Foam::constant::thermodynamic::RR*Tc/Pc;

    const scalar a =
        0.42747*sqr(Foam::constant::thermodynamic::RR*Tc)/Pc;

    const scalar S = 0.48508 + 1.5517*omega - 0.15613*sqr(omega);

    coef1_ = a*sqr(1.0 + S);

    coef2_ = a*2*S*(1 + S)/sqrt(Tc);

    coef3_ = a*sqr(S)/Tc;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::SRKGas<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const SRKGas<Specie>& srk
)
{
    srk.write(os);
    return os;
}


// ************************************************************************* //
