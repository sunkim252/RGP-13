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
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::Xterm
(
    const scalar p,
    const scalar T
) const
{
    // Cache Z to avoid redundant cubic solve in dZdT
    const scalar Zv = Z(p, T);
    const scalar B = b_*p/(Foam::constant::thermodynamic::RR*T);
    const scalar sqrtT = sqrt(T);
    const scalar first =
        3*coef2_*log((Zv + B)/Zv)/(4*b_*sqrtT);

    return
    (
        first
      - (
            -coef1_/b_ + 0.5*coef2_*sqrtT/b_
        )
       *(
            b_*p
           /((Zv + B)*Foam::constant::thermodynamic::RR*T*T)
          + (B/(sqr(Zv) + Zv*B))*dZdT(p, T, Zv)
        )
    );
}


template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::dZdT
(
    const scalar p,
    const scalar T
) const
{
    return dZdT(p, T, Z(p, T));
}


template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::dZdT
(
    const scalar p,
    const scalar T,
    const scalar Zv
) const
{
    const scalar aAlpha = coef1_ - coef2_*sqrt(T) + coef3_*T;
    const scalar A =
        aAlpha*p/sqr(Foam::constant::thermodynamic::RR*T);
    const scalar B = b_*p/(Foam::constant::thermodynamic::RR*T);

    const scalar numerator =
       -(A + Zv + 2*B*Zv)*b_*p
       /(Foam::constant::thermodynamic::RR*T*T)
      + ((B - Zv)*(coef3_*T - coef1_)/(2*T))
       *(p/sqr(Foam::constant::thermodynamic::RR*T));

    return numerator/(3*sqr(Zv) - 2*Zv + A - B - sqr(B));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::SRKGas<Specie>::SRKGas
(
    const Specie& sp,
    const scalar& b,
    const scalar& coef1,
    const scalar& coef2,
    const scalar& coef3
)
:
    Specie(sp),
    b_(b),
    coef1_(coef1),
    coef2_(coef2),
    coef3_(coef3)
{}


template<class Specie>
inline Foam::SRKGas<Specie>::SRKGas
(
    const word& name,
    const SRKGas& srk
)
:
    Specie(name, srk),
    b_(srk.b_),
    coef1_(srk.coef1_),
    coef2_(srk.coef2_),
    coef3_(srk.coef3_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::SRKGas<Specie>>
Foam::SRKGas<Specie>::clone() const
{
    return autoPtr<SRKGas<Specie>>
    (
        new SRKGas<Specie>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    const scalar Zv = this->Z(p, T);
    return p/(Zv*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::h
(
    const scalar p,
    const scalar T
) const
{
    scalar B = b_*p/(Foam::constant::thermodynamic::RR*T);
    if (B <= 0)
    {
        B = 1e-16;
    }

    const scalar sqrtT = sqrt(T);
    const scalar aAlpha = coef1_ - coef2_*sqrtT + coef3_*T;
    const scalar daAlpha = -coef2_/(2*sqrtT) + coef3_;

    const scalar Zv = this->Z(p, T);

    if (B == -Zv)
    {
        return
            (Foam::constant::thermodynamic::RR*T*(Zv - 1))
           /this->W();
    }
    else
    {
        return
            (
                Foam::constant::thermodynamic::RR*T*(Zv - 1)
              + (T*daAlpha - aAlpha)/b_*log((Zv + B)/Zv)
            )/this->W();
    }
}


template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::Cp
(
    scalar p,
    scalar T
) const
{
    const scalar sqrtT = sqrt(T);
    const scalar aAlpha = coef1_ - coef2_*sqrtT + coef3_*T;
    const scalar daAlpha = -coef2_/(2*sqrtT) + coef3_;
    const scalar ddaAlpha = coef2_/(4*T*sqrtT);
    const scalar A =
        aAlpha*p/sqr(Foam::constant::thermodynamic::RR*T);
    scalar B = b_*p/(Foam::constant::thermodynamic::RR*T);
    if (B <= 0)
    {
        B = 1e-16;
    }

    const scalar Zv = this->Z(p, T);

    const scalar M = (sqr(Zv) + B*Zv)/(Zv - B);
    const scalar N = daAlpha*B/(b_*Foam::constant::thermodynamic::RR);

    return
    (
        (T/b_)*ddaAlpha*log((Zv + B)/Zv)
      + Foam::constant::thermodynamic::RR*sqr(M - N)
       /(sqr(M) - A*(2*Zv + B))
      - Foam::constant::thermodynamic::RR
    )/this->W();
}


template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::e
(
    const scalar p,
    const scalar T
) const
{
    // e_departure = h_departure without the (Z-1) ideal-gas PV work term
    // Following the PengRobinsonGas pattern: e = h - RR*T*(Z-1)/W
    scalar B = b_*p/(Foam::constant::thermodynamic::RR*T);
    if (B <= 0)
    {
        B = 1e-16;
    }

    const scalar Zv = this->Z(p, T);

    if (B == -Zv)
    {
        // Departure vanishes
        return 0;
    }
    else
    {
        const scalar sqrtT = sqrt(T);
        const scalar aAlpha = coef1_ - coef2_*sqrtT + coef3_*T;
        const scalar daAlpha = -coef2_/(2*sqrtT) + coef3_;

        return
            (
                (T*daAlpha - aAlpha)/b_*log((Zv + B)/Zv)
            )/this->W();
    }
}


template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::Cv
(
    scalar p,
    scalar T
) const
{
    // Cv = Cp - CpMCv, inline both to share a single Z(p,T) solve
    const scalar sqrtT = sqrt(T);
    const scalar aAlpha = coef1_ - coef2_*sqrtT + coef3_*T;
    const scalar daAlpha = -coef2_/(2*sqrtT) + coef3_;
    const scalar ddaAlpha = coef2_/(4*T*sqrtT);
    const scalar A =
        aAlpha*p/sqr(Foam::constant::thermodynamic::RR*T);
    scalar B = b_*p/(Foam::constant::thermodynamic::RR*T);
    if (B <= 0)
    {
        B = 1e-16;
    }

    const scalar Zv = this->Z(p, T);

    const scalar M = (sqr(Zv) + B*Zv)/(Zv - B);
    const scalar N = daAlpha*B/(b_*Foam::constant::thermodynamic::RR);

    const scalar MminusN = M - N;
    const scalar Msq = sqr(M);
    const scalar denom = Msq - A*(2*Zv + B);

    // Cp (departure) contribution
    const scalar CpVal =
    (
        (T/b_)*ddaAlpha*log((Zv + B)/Zv)
      + Foam::constant::thermodynamic::RR*sqr(MminusN)/denom
      - Foam::constant::thermodynamic::RR
    )/this->W();

    // CpMCv contribution
    const scalar CpMCvVal = this->R()*sqr(MminusN)/denom;

    return CpVal - CpMCvVal;
}


template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::sp
(
    const scalar p,
    const scalar T
) const
{
    const scalar daAlpha = -coef2_/(2*sqrt(T)) + coef3_;
    scalar B = b_*p/(Foam::constant::thermodynamic::RR*T);
    if (B <= 0)
    {
        B = 1e-16;
    }

    const scalar Zv = this->Z(p, T);

    return
    (
        Foam::constant::thermodynamic::RR
       *(
            -log(p/Pstd)
          + log(Zv - B)
        )
      - daAlpha/b_*log((Zv + B)/Zv)
    )/this->W();
}


template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::sv
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    const scalar Zv = this->Z(p, T);
    return 1.0/(Zv*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    const scalar aAlpha = coef1_ - coef2_*sqrt(T) + coef3_*T;

    const scalar A =
        aAlpha*p/sqr(Foam::constant::thermodynamic::RR*T);
    const scalar B = b_*p/(Foam::constant::thermodynamic::RR*T);

    const scalar a2 = -1.0;
    const scalar a1 = A - B - sqr(B);
    const scalar a0 = -A*B;

    const scalar Q = (3*a1 - a2*a2)/9.0;
    const scalar Rl = (9*a2*a1 - 27*a0 - 2*a2*a2*a2)/54.0;

    const scalar Q3 = Q*Q*Q;
    const scalar D = Q3 + Rl*Rl;

    scalar root = -1;

    if (D <= 0)
    {
        const scalar th = ::acos(Rl/sqrt(-Q3));
        const scalar qm = 2*sqrt(-Q);
        const scalar r1 = qm*cos(th/3.0) - a2/3.0;
        const scalar r2 =
            qm*cos((th + 2*constant::mathematical::pi)/3.0) - a2/3.0;
        const scalar r3 =
            qm*cos((th + 4*constant::mathematical::pi)/3.0) - a2/3.0;

        root = max(r1, max(r2, r3));
    }
    else
    {
        // One root is real
        const scalar D05 = sqrt(D);
        scalar S = 0;

        if (Rl + D05 < 0)
        {
            S = -cbrt(mag(Rl + D05));
        }
        else
        {
            S = cbrt(Rl + D05);
        }
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -cbrt(mag(Rl - D05));
        }
        else
        {
            Tl = cbrt(Rl - D05);
        }

        root = S + Tl - a2/3.0;
    }

    return root;
}


template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    const scalar aAlpha = coef1_ - coef2_*sqrt(T) + coef3_*T;
    const scalar daAlpha = -coef2_/(2*sqrt(T)) + coef3_;

    const scalar A =
        aAlpha*p/sqr(Foam::constant::thermodynamic::RR*T);
    scalar B = b_*p/(Foam::constant::thermodynamic::RR*T);
    if (B <= 0)
    {
        B = 1e-16;
    }

    const scalar Zv = this->Z(p, T);

    const scalar M = (sqr(Zv) + B*Zv)/(Zv - B);
    const scalar N = daAlpha*B/(b_*Foam::constant::thermodynamic::RR);

    return this->R()*sqr(M - N)/(sqr(M) - A*(2*Zv + B));
}


template<class Specie>
inline Foam::scalar Foam::SRKGas<Specie>::alphav
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::SRKGas<Specie>::operator+=
(
    const SRKGas<Specie>& srk
)
{
    scalar X1 = this->Y()/this->W();
    Specie::operator+=(srk);

    if (mag(this->Y()) > small)
    {
        X1 *= this->W()/this->Y();
        const scalar X2 = this->W()*srk.Y()/(srk.W()*this->Y());

        b_ = X1*b_ + X2*srk.b_;
        coef1_ = X1*coef1_ + X2*srk.coef1_;
        coef2_ = X1*coef2_ + X2*srk.coef2_;
        coef3_ = X1*coef3_ + X2*srk.coef3_;
    }
}


template<class Specie>
inline void Foam::SRKGas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
Foam::SRKGas<Specie> Foam::operator+
(
    const SRKGas<Specie>& srk1,
    const SRKGas<Specie>& srk2
)
{
    Specie sp
    (
        static_cast<const Specie&>(srk1)
      + static_cast<const Specie&>(srk2)
    );

    if (mag(sp.Y()) < small)
    {
        return SRKGas<Specie>
        (
            sp,
            srk1.b_,
            srk1.coef1_,
            srk1.coef2_,
            srk1.coef3_
        );
    }
    else
    {
        const scalar X1 = sp.W()*srk1.Y()/(srk1.W()*sp.Y());
        const scalar X2 = sp.W()*srk2.Y()/(srk2.W()*sp.Y());

        return SRKGas<Specie>
        (
            sp,
            X1*srk1.b_ + X2*srk2.b_,
            X1*srk1.coef1_ + X2*srk2.coef1_,
            X1*srk1.coef2_ + X2*srk2.coef2_,
            X1*srk1.coef3_ + X2*srk2.coef3_
        );
    }
}


template<class Specie>
Foam::SRKGas<Specie> Foam::operator*
(
    const scalar s,
    const SRKGas<Specie>& srk
)
{
    return SRKGas<Specie>
    (
        s*static_cast<const Specie&>(srk),
        srk.b_,
        srk.coef1_,
        srk.coef2_,
        srk.coef3_
    );
}


template<class Specie>
Foam::SRKGas<Specie> Foam::operator==
(
    const SRKGas<Specie>& srk1,
    const SRKGas<Specie>& srk2
)
{
    Specie sp
    (
        static_cast<const Specie&>(srk1)
     == static_cast<const Specie&>(srk2)
    );

    return SRKGas<Specie>
    (
        sp,
        NaN,
        NaN,
        NaN,
        NaN
    );
}


// ************************************************************************* //
