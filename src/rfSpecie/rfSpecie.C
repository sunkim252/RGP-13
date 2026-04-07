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

\*---------------------------------------------------------------------------*/

#include "rfSpecie.H"
#include "constants.H"

/* * * * * * * * * * * * * * * public constants  * * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(rfSpecie, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rfSpecie::rfSpecie(const word& name, const dictionary& dict)
:
    name_(name),
    Y_(dict.subDict("specie").lookupOrDefault("massFraction", 1.0)),
    molWeight_(dict.subDict("specie").lookup<scalar>("molWeight")),
    Tc_(dict.subDict("rfProperties").lookup<scalar>("Tc")),
    Pc_(dict.subDict("rfProperties").lookup<scalar>("Pc")),
    Vc_(dict.subDict("rfProperties").lookup<scalar>("Vc")),
    omega_(dict.subDict("rfProperties").lookup<scalar>("omega")),
    miui_(dict.subDict("rfProperties").lookup<scalar>("miui")),
    kappai_(dict.subDict("rfProperties").lookup<scalar>("kappai")),
    sigmvi_(dict.subDict("rfProperties").lookup<scalar>("sigmvi"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rfSpecie::write(Ostream& os) const
{
    dictionary dict("specie");
    if (Y_ != 1)
    {
        dict.add("massFraction", Y_);
    }
    dict.add("molWeight", molWeight_);
    os  << indent << dict.dictName() << dict;

    dictionary rfDict("rfProperties");
    rfDict.add("Tc", Tc_);
    rfDict.add("Pc", Pc_);
    rfDict.add("Vc", Vc_);
    rfDict.add("omega", omega_);
    rfDict.add("miui", miui_);
    rfDict.add("kappai", kappai_);
    rfDict.add("sigmvi", sigmvi_);
    os  << indent << rfDict.dictName() << rfDict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const rfSpecie& st)
{
    st.write(os);
    os.check("Ostream& operator<<(Ostream& os, const rfSpecie& st)");
    return os;
}


// ************************************************************************* //
