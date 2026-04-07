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

#include "chungTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::chungTransport<Thermo>::chungTransport
(
    const word& name,
    const dictionary& dict
)
:
    Thermo(name, dict),
    sigmaM_(0.809*pow(this->Vc(), 1.0/3)),
    epsilonkM_(this->Tc()/1.2593),
    MM_(this->W()),
    VcM_(this->Vc()),
    TcM_(this->Tc()),
    omegaM_(this->omega()),
    miuiM_(this->miui()),
    kappaiM_(this->kappai()),
    Ymd_(2),
    Xmd_(2),
    Tcmd_(2),
    Pcmd_(2),
    Mmd_(2),
    sigmd_(2)
{
    // Temporary initialization
    forAll(Ymd_, i)
    {
        Ymd_[i] = this->Y();
    }

    forAll(Xmd_, i)
    {
        Xmd_[i] = this->Y()/this->W();
    }

    List<scalar> TcTemp(2);
    forAll(TcTemp, i)
    {
        TcTemp[i] = this->Tc();
    }
    forAll(Tcmd_, i)
    {
        Tcmd_[i] = TcTemp;
    }

    List<scalar> PcTemp(2);
    forAll(PcTemp, i)
    {
        PcTemp[i] = this->Pc();
    }
    forAll(Pcmd_, i)
    {
        Pcmd_[i] = PcTemp;
    }

    List<scalar> MTemp(2);
    forAll(MTemp, i)
    {
        MTemp[i] = this->W();
    }
    forAll(Mmd_, i)
    {
        Mmd_[i] = MTemp;
    }

    List<scalar> sigTemp(2);
    forAll(sigTemp, i)
    {
        sigTemp[i] = this->sigmvi();
    }
    forAll(sigmd_, i)
    {
        sigmd_[i] = sigTemp;
    }
}


template<class Thermo>
Foam::chungTransport<Thermo>::chungTransport
(
    const Thermo& t,
    const dictionary& dict
)
:
    Thermo(t),
    sigmaM_(0.809*pow(this->Vc(), 1.0/3)),
    epsilonkM_(this->Tc()/1.2593),
    MM_(this->W()),
    VcM_(this->Vc()),
    TcM_(this->Tc()),
    omegaM_(this->omega()),
    miuiM_(this->miui()),
    kappaiM_(this->kappai()),
    Ymd_(2),
    Xmd_(2),
    Tcmd_(2),
    Pcmd_(2),
    Mmd_(2),
    sigmd_(2)
{
    // Temporary initialization
    forAll(Ymd_, i)
    {
        Ymd_[i] = this->Y();
    }

    forAll(Xmd_, i)
    {
        Xmd_[i] = this->Y()/this->W();
    }

    List<scalar> TcTemp(2);
    forAll(TcTemp, i)
    {
        TcTemp[i] = this->Tc();
    }
    forAll(Tcmd_, i)
    {
        Tcmd_[i] = TcTemp;
    }

    List<scalar> PcTemp(2);
    forAll(PcTemp, i)
    {
        PcTemp[i] = this->Pc();
    }
    forAll(Pcmd_, i)
    {
        Pcmd_[i] = PcTemp;
    }

    List<scalar> MTemp(2);
    forAll(MTemp, i)
    {
        MTemp[i] = this->W();
    }
    forAll(Mmd_, i)
    {
        Mmd_[i] = MTemp;
    }

    List<scalar> sigTemp(2);
    forAll(sigTemp, i)
    {
        sigTemp[i] = this->sigmvi();
    }
    forAll(sigmd_, i)
    {
        sigmd_[i] = sigTemp;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::chungTransport<Thermo>::write(Ostream& os) const
{
    os  << this->name() << endl
        << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("sigmaM", sigmaM_);
    dict.add("epsilonkM", epsilonkM_);
    dict.add("MM", MM_);
    dict.add("VcM", VcM_);
    dict.add("TcM", TcM_);
    dict.add("omegaM", omegaM_);
    dict.add("miuiM", miuiM_);
    dict.add("kappaiM", kappaiM_);

    os  << indent << dict.dictName() << dict
        << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const chungTransport<Thermo>& ch
)
{
    ch.write(os);
    return os;
}


// ************************************************************************* //
