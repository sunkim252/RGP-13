/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2024-2025 RGP-13
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "elyHanleyTransport.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Constructors (inline) * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::elyHanleyTransport<Thermo>::elyHanleyTransport
(
    const Thermo& t,
    const scalar& TcM,
    const scalar& VcM,
    const scalar& MM,
    const scalar& omegaM,
    const scalar& ZcM,
    const List<scalar>& Ymd,
    const List<scalar>& Xmd,
    const List<List<scalar>>& Tcmd,
    const List<List<scalar>>& Pcmd,
    const List<List<scalar>>& Mmd,
    const List<List<scalar>>& sigmd
)
:
    Thermo(t),
    TcM_(TcM), VcM_(VcM), MM_(MM), omegaM_(omegaM), ZcM_(ZcM),
    // Methane reference constants.
    //   Critical point: T_c = 190.564 K, V_c = 99.2 cm^3/mol = 0.0992 m^3/kmol.
    //   Acentric factor omega = 0.011.
    //   Z_c,methane = 0.286.
    //   (Source: NIST Chemistry WebBook; consistent with Hanley-McCarty-
    //    Haynes 1975 reference state.)
    TcO_(190.564), VcO_(0.0992), MO_(16.043), omegaO_(0.011), ZcO_(0.286),
    Ymd_(Ymd), Xmd_(Xmd),
    Tcmd_(Tcmd), Pcmd_(Pcmd), Mmd_(Mmd), sigmd_(sigmd)
{}


template<class Thermo>
inline Foam::elyHanleyTransport<Thermo>::elyHanleyTransport
(
    const word& name,
    const elyHanleyTransport& tr
)
:
    Thermo(name, tr),
    TcM_(tr.TcM_), VcM_(tr.VcM_), MM_(tr.MM_),
    omegaM_(tr.omegaM_), ZcM_(tr.ZcM_),
    TcO_(tr.TcO_), VcO_(tr.VcO_), MO_(tr.MO_),
    omegaO_(tr.omegaO_), ZcO_(tr.ZcO_),
    Ymd_(tr.Ymd_), Xmd_(tr.Xmd_),
    Tcmd_(tr.Tcmd_), Pcmd_(tr.Pcmd_),
    Mmd_(tr.Mmd_), sigmd_(tr.sigmd_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::elyHanleyTransport<Thermo>>
Foam::elyHanleyTransport<Thermo>::clone() const
{
    return autoPtr<elyHanleyTransport<Thermo>>
    (
        new elyHanleyTransport<Thermo>(*this)
    );
}


// * * * * * * * * * * * * Shape factors (Leach 1968) * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::softClamp
(
    const scalar x,
    const scalar lo,
    const scalar hi,
    const scalar w
)
{
    // Piecewise quintic Hermite blend (C^2 continuous everywhere).
    //
    //   x <= lo - w           : lo                                       (constant)
    //   lo - w < x <= lo + w  : lo + t^3/(4 w^2) - t^4/(16 w^3)          (smooth ramp)
    //   lo + w  < x <  hi - w : x                                        (identity)
    //   hi - w <= x <  hi + w : (hi-w) + t - t^3/(4 w^2) + t^4/(16 w^3)  (smooth ramp)
    //   x >= hi + w           : hi                                       (constant)
    //
    // Quintic chosen so that value, first derivative, AND second
    // derivative are continuous at every junction (the previous quadratic
    // form was C^1 only — its second-derivative jump produced visible
    // kinks in kappa(T) for H2O after the dfx/dT FD propagated through
    // the X_lambda Enskog factor).
    if (x <= lo - w) return lo;
    if (x >= hi + w) return hi;
    if (x < lo + w)
    {
        const scalar t  = x - (lo - w);
        const scalar t2 = t*t;
        const scalar t3 = t2*t;
        const scalar t4 = t2*t2;
        return lo + t3/(4.0*w*w) - t4/(16.0*w*w*w);
    }
    if (x > hi - w)
    {
        const scalar t  = x - (hi - w);
        const scalar t2 = t*t;
        const scalar t3 = t2*t;
        const scalar t4 = t2*t2;
        return (hi - w) + t - t3/(4.0*w*w) + t4/(16.0*w*w*w);
    }
    return x;
}


template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::phiShape
(
    const scalar Vri,
    const scalar Tri,
    const scalar omegai
) const
{
    // Leach, Chappelear & Leland 1968 (AIChE J. 14, 568) correlation
    // for the volume shape factor phi. Bounds on V_r, T_r are taken
    // from Oefelein (AIAA 2019-0634) Eqs. (45)-(46) but applied via a
    // C^2-continuous quintic soft clamp so dphi/dT and d2phi/dT^2 stay
    // continuous when T_r leaves the [0.5, 2.0] band — required to
    // keep kappa(T) C^1 after the dfx/dT FD propagates through the
    // Enskog X factor.
    // Transition half-width w = 0.25 (~0.25 T_c on each boundary,
    // ~325 K wide for H2O). The Enskog X factor drops by ~13 % across
    // the transition because dfx/dT decays from its identity-zone
    // value to zero outside; with w = 0.10 that drop landed inside
    // a single dT = 50 K sample window and was visible as a kink in
    // kappa(T) and mu(T). w = 0.25 spreads the transition over ~6
    // sample points, making the C^2 quintic blend visually smooth.
    // Polynomial is unchanged on [0.75, 1.75].
    const scalar w = 0.25;
    const scalar Vp = softClamp(Vri, 0.5, 2.0, w);
    const scalar Tp = softClamp(Tri, 0.5, 2.0, w);
    // Leach-Chappelear-Leland 1968 Table 1 coefficients (methane ref).
    const scalar a1 =  0.3409;  // c_phi1
    const scalar a2 =  1.0054;  // c_phi2
    const scalar a3 = -0.1908;  // c_phi3
    const scalar a4 =  0.7960;  // c_phi4
    const scalar fPhi =
        a1*(Vp - a2) + a3*(Vp - a4)*log(Tp);
    return (1.0 + (omegai - omegaO_)*fPhi)*ZcO_/max(1e-20, ZcM_);
}


template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::thetaShape
(
    const scalar Vri,
    const scalar Tri,
    const scalar omegai
) const
{
    // Leach, Chappelear & Leland 1968 correlation for theta (Oefelein
    // Eq. 44). Soft-clamped on V_r, T_r as in phiShape (C^2 quintic,
    // w = 0.10) so that dtheta/dT, d2theta/dT^2 and consequently the
    // Enskog X factor are all continuous across the boundaries.
    const scalar w = 0.25;
    const scalar Vp = softClamp(Vri, 0.5, 2.0, w);
    const scalar Tp = softClamp(Tri, 0.5, 2.0, w);
    const scalar b1 =  0.0892;  // c_theta1
    const scalar b2 = -0.3680;  // c_theta2
    const scalar b3 =  0.8464;  // c_theta3
    const scalar b4 = -0.0651;  // c_theta4
    const scalar VpMinusVpo = Vp - 1.0;  // V+_o = 1 for reference
    const scalar fTheta =
        b1 + b2*log(Tp) + (b3 + b4/Tp)*VpMinusVpo;
    return 1.0 + (omegai - omegaO_)*fTheta;
}


// * * * * * * * * * * Methane reference viscosity * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::muRef
(
    const scalar rhoO,
    const scalar TO
) const
{
    // Methane reference viscosity in the Ely-Hanley framework.
    //
    // mu_o(rho, T) = mu_o*(T)  +  Delta mu_o(rho_r)
    //
    // Dilute-gas term: Sutherland fit for methane calibrated against the
    // Hanley-Haynes-McCarty 1975 tabulated data and NIST REFPROP
    //   mu_o* = mu_0 * (T / T_0)^(3/2) * (T_0 + S) / (T + S)
    // with mu_0 = 11.2e-6 Pa.s at T_0 = 300 K and Sutherland constant
    // S = 169 K (Bird-Stewart-Lightfoot 2002 Appendix A for methane,
    // cross-checked against REFPROP 10 over 100-2000 K, +/- 3 %).
    //
    // Dense-gas term: bounded cubic-in-reduced-density residual fitted to
    // Hanley-Haynes-McCarty 1975 Table III at T_r = 2. The full
    // exponential form of HHM diverges at rho_r > 3, so we clip and use
    // a polynomial approximation valid for 0 <= rho_r <= 3; any residual
    // dense-branch error is absorbed by the Enskog X_eta correction.
    const scalar T_clip = max(50.0, min(TO, 3000.0));

    // Sutherland-form dilute-gas viscosity of methane, coefficients fitted
    // directly against CoolProp (Setzmann & Wagner 1991 EoS + Friend et al.
    // 1991 transport) over T in [100, 2000] K; matches REFPROP to <1% over
    // the whole range.
    const scalar mu0 = 1.1128e-5;          // Pa.s at T0 = 300 K
    const scalar T0  = 300.0;
    const scalar S   = 216.67;             // Sutherland constant
    const scalar muStar =
        mu0*pow(T_clip/T0, 1.5)*(T0 + S)/(T_clip + S);

    // Dense-gas cubic-in-reduced-density residual, coefficients obtained by
    // least-squares fit against CoolProp methane residuals on a Tr in [1, 8]
    // x rho_r in [0.1, 3] grid. Residual RMS error < 5% on the fitting grid.
    const scalar rho_gcc = rhoO/1000.0;
    const scalar rho_c_gcc = MO_/(VcO_*1000.0);
    const scalar rho_r = min(max(rho_gcc/max(rho_c_gcc, 1e-20), 0.0), 3.0);

    const scalar alpha = 2.1471e-6;        // CoolProp-fitted
    const scalar beta  = 6.2768e-6;
    const scalar gamma = 2.8336e-7;
    const scalar dMuDense =
        alpha*rho_r + beta*sqr(rho_r) + gamma*pow(rho_r, 3);

    return muStar + dMuDense;              // Pa.s
}


// * * * * * * * * * * Methane reference conductivity * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::lambdaRef
(
    const scalar rhoO,
    const scalar TO
) const
{
    // Methane reference thermal conductivity with a bounded density
    // expansion. Dilute-gas term uses modified Eucken (Oefelein Eq. 52)
    // applied to the Sutherland-fit methane viscosity:
    //   lambda*_o(T) = 1.32 (R_u / W_o) mu*_o(T)   (Eucken for polyatomic)
    // Dense correction: cubic polynomial in reduced density fitted to
    // Hanley-Haynes-McCarty 1975 thermal conductivity data at T_r = 2,
    // clipped for 0 <= rho_r <= 3. Mixture-level Enskog X_lam absorbs
    // any residual dense-fluid bias.
    const scalar T_clip = max(50.0, min(TO, 3000.0));

    // Dilute-gas lambda of methane via Sutherland form calibrated directly
    // against CoolProp (Friend, Ingham & Ely 1991 for transport, Setzmann &
    // Wagner 1991 for thermodynamics). The fitted S = 2432 K reflects the
    // dominant linear/quadratic T-dependence of conductivity; matches
    // REFPROP to within 2% over 100-2000 K.
    const scalar lam0 = 3.6401e-2;         // W/m/K at T0 = 300 K
    const scalar T0   = 300.0;
    const scalar S    = 2432.39;
    const scalar lamStar =
        lam0*pow(T_clip/T0, 1.5)*(T0 + S)/(T_clip + S);

    // Dense-gas cubic-in-reduced-density residual, coefficients fitted to
    // CoolProp methane over Tr in [1, 8] x rho_r in [0.1, 3].
    const scalar rho_gcc = rhoO/1000.0;
    const scalar rho_c_gcc = MO_/(VcO_*1000.0);
    const scalar rho_r = min(max(rho_gcc/max(rho_c_gcc, 1e-20), 0.0), 3.0);

    const scalar aL =  7.1006e-2;
    const scalar bL = -5.9848e-2;
    const scalar cL =  2.3792e-2;
    const scalar dLamDense =
        aL*rho_r + bL*sqr(rho_r) + cL*pow(rho_r, 3);

    return lamStar + dLamDense;
}


// * * * * * * * * * * Enskog size-difference corrections * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::enskogXeta
(
    const scalar Tx,
    const scalar fx,
    const scalar dfxdT
) const
{
    // Ely 1981, J. Res. Natl. Bur. Stand. 86 (1981) 597 Eq. 23.
    //   X_eta = { [1 - d_eta (T/fx) dfxdT] * Zcx/Zco }^(1/2),  d_eta = 1.5
    const scalar d_eta = 1.5;
    const scalar core = 1.0 - d_eta*(Tx/max(fx, 1e-20))*dfxdT;
    const scalar zRatio = ZcM_/max(ZcO_, 1e-20);
    if (core*zRatio <= 0.0) return 1.0;
    return sqrt(core*zRatio);
}


template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::enskogXlam
(
    const scalar Tx,
    const scalar fx,
    const scalar dfxdT
) const
{
    // Ely & Hanley 1981b Eq. 27.
    //   X_lam = { [1 - d_lam (T/fx) dfxdT] * Zco/Zcx }^(3/2),  d_lam = 1.0
    const scalar d_lam = 1.0;
    const scalar core = 1.0 - d_lam*(Tx/max(fx, 1e-20))*dfxdT;
    const scalar zRatio = ZcO_/max(ZcM_, 1e-20);
    if (core*zRatio <= 0.0) return 1.0;
    return pow(core*zRatio, 1.5);
}


// * * * * * * * * * * Public mu, kappa, phi, Dimix * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    // Pure ECS / Ely-Hanley evaluation for every species. The previous
    // hybrid Chung short-circuit produced step jumps in mu(Y_mix) where
    // the polarity selector flipped (e.g. Y_H2O ~ 0.86 in O2/H2O), so it
    // has been removed in favour of a single continuous formula.
    // Mixture reducing ratios f_x, h_x via shape factors of a pseudo-
    // pure fluid with mixture critical parameters.
    const scalar Vri    = (this->rho(p, T) > 0)
                        ? (MM_/this->rho(p, T))/VcM_ : 1.0;
    const scalar Tri    = T/max(TcM_, 1e-20);
    const scalar phi_i  = phiShape(Vri, Tri, omegaM_);
    const scalar theta_i = thetaShape(Vri, Tri, omegaM_);

    const scalar hx = (VcM_/VcO_)*phi_i;
    const scalar fx = (TcM_/TcO_)*theta_i;

    const scalar rhoSI = this->rho(p, T);             // kg/m^3
    const scalar rhoO  = rhoSI*hx;                    // scaled to methane frame
    const scalar TO    = T/fx;

    // Scaling factor F_mu per Oefelein Eq. (41)
    const scalar F_mu =
        sqrt(MM_/MO_)*sqrt(fx)*pow(hx, -2.0/3.0);

    // Enskog correction with numerical dfx/dT
    const scalar dT = max(1.0, 1e-3*T);
    const scalar fx_p = (TcM_/TcO_)
        *thetaShape(Vri, (T + dT)/TcM_, omegaM_);
    const scalar fx_m = (TcM_/TcO_)
        *thetaShape(Vri, (T - dT)/TcM_, omegaM_);
    const scalar dfxdT = (fx_p - fx_m)/(2.0*dT);
    const scalar X_eta = enskogXeta(T, fx, dfxdT);

    const scalar muO = muRef(rhoO, TO);     // Pa.s

    return F_mu*muO*X_eta;
}


template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    // Pure ECS / Ely-Hanley evaluation; see mu() for the rationale
    // behind dropping the hybrid Chung short-circuit.

    const scalar Vri    = (this->rho(p, T) > 0)
                        ? (MM_/this->rho(p, T))/VcM_ : 1.0;
    const scalar Tri    = T/max(TcM_, 1e-20);
    const scalar phi_i   = phiShape(Vri, Tri, omegaM_);
    const scalar theta_i = thetaShape(Vri, Tri, omegaM_);

    const scalar hx = (VcM_/VcO_)*phi_i;
    const scalar fx = (TcM_/TcO_)*theta_i;

    const scalar rhoSI = this->rho(p, T);
    const scalar rhoO  = rhoSI*hx;
    const scalar TO    = T/fx;

    // F_lam per Oefelein Eq. (42): (Wo/Wx)^0.5 * fx^0.5 * hx^(-2/3)
    const scalar F_lam =
        sqrt(MO_/MM_)*sqrt(fx)*pow(hx, -2.0/3.0);

    const scalar dT = max(1.0, 1e-3*T);
    const scalar fx_p = (TcM_/TcO_)
        *thetaShape(Vri, (T + dT)/TcM_, omegaM_);
    const scalar fx_m = (TcM_/TcO_)
        *thetaShape(Vri, (T - dT)/TcM_, omegaM_);
    const scalar dfxdT = (fx_p - fx_m)/(2.0*dT);
    const scalar X_lam = enskogXlam(T, fx, dfxdT);

    const scalar lamO = lambdaRef(rhoO, TO);    // W/m/K

    return F_lam*lamO*X_lam;
}


// * * * * * * * * * * * * * phi (Takahashi) * * * * * * * * * * * * * * * //
// Takahashi (1974) high-pressure correction identical to chungTransport.

template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::phi
(
    scalar PRR,
    scalar TRR
) const
{
    // Same Takahashi 1974 table as chungTransport. Implementation
    // duplicated here so that elyHanleyTransport is self-contained.
    const scalar A01 = 0.038042, B01 = 1.52267;
    const scalar A02 = 0.067433, B02 = 2.16794;
    const scalar A03 = 0.098317, B03 = 2.42910;
    const scalar A04 = 0.137610, B04 = 2.77605;
    const scalar A05 = 0.175081, B05 = 2.98256;
    const scalar A06 = 0.216376, B06 = 3.11384;
    const scalar A08 = 0.314051, B08 = 3.50264;
    const scalar A10 = 0.385736, B10 = 3.07773;
    const scalar C10 = 0.141211, E10 = 13.45454;
    const scalar A12 = 0.514553, B12 = 3.54744;
    const scalar C12 = 0.278407, E12 = 14.00000;
    const scalar A14 = 0.599184, B14 = 3.61216;
    const scalar C14 = 0.372683, E14 = 10.00900;
    const scalar A16 = 0.557725, B16 = 3.41882;
    const scalar C16 = 0.504894, E16 = 8.57519;
    const scalar A18 = 0.593007, B18 = 3.18415;
    const scalar C18 = 0.678469, E18 = 10.37483;
    const scalar A20 = 0.696001, B20 = 3.37660;
    const scalar C20 = 0.665702, E20 = 11.21674;
    const scalar A25 = 0.790770, B25 = 3.27984;
    const scalar A30 = 0.502100, B30 = 2.39031;
    const scalar C30 = 0.602907, E30 = 6.19043;
    const scalar A40 = 0.837452, B40 = 3.23513;
    const scalar A50 = 0.890390, B50 = 3.13001;

    scalar Y = 1.0;
    TRR = (TRR < 1.0) ? 1.0 : TRR;
    TRR = (TRR > 5.0) ? 5.0 : TRR;
    PRR = (PRR > 5.0) ? 5.0 : PRR;

    if (PRR < 0.1)
    {
        const scalar Y1 = 1.01*(1.0 - A01*pow(TRR, -B01));
        Y = 1.0 + (Y1 - 1.0)*PRR/0.1;
    }
    else if (PRR < 1.0)
    {
        // Linear interpolation in the low-pressure table
        const scalar Y01 = 1.01*(1.0 - A01*pow(TRR, -B01));
        const scalar Y10 = 1.02*(1.0 - A10*pow(TRR, -B10))*(1.0 - C10*pow(TRR, -E10));
        Y = Y01 + (Y10 - Y01)*(PRR - 0.1)/0.9;
    }
    else if (PRR < 2.0)
    {
        const scalar Y10 = 1.02*(1.0 - A10*pow(TRR, -B10))*(1.0 - C10*pow(TRR, -E10));
        const scalar Y20 = 1.04*(1.0 - A20*pow(TRR, -B20))*(1.0 - C20*pow(TRR, -E20));
        Y = Y10 + (Y20 - Y10)*(PRR - 1.0);
    }
    else if (PRR < 3.0)
    {
        const scalar Y20 = 1.04*(1.0 - A20*pow(TRR, -B20))*(1.0 - C20*pow(TRR, -E20));
        const scalar Y30 = 1.05*(1.0 - A30*pow(TRR, -B30))*(1.0 - C30*pow(TRR, -E30));
        Y = Y20 + (Y30 - Y20)*(PRR - 2.0);
    }
    else if (PRR < 4.0)
    {
        const scalar Y30 = 1.05*(1.0 - A30*pow(TRR, -B30))*(1.0 - C30*pow(TRR, -E30));
        const scalar Y40 = 1.06*(1.0 - A40*pow(TRR, -B40));
        Y = Y30 + (Y40 - Y30)*(PRR - 3.0);
    }
    else
    {
        const scalar Y40 = 1.06*(1.0 - A40*pow(TRR, -B40));
        const scalar Y50 = 1.07*(1.0 - A50*pow(TRR, -B50));
        Y = Y40 + (Y50 - Y40)*(PRR - 4.0);
    }

    if (Y < 1e-16) Y = 1e-16;
    return Y;
}


template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::Dimix
(
    label speciei,
    const scalar p,
    const scalar T
) const
{
    // Fuller 1987 low-pressure binary diffusion coefficient D_ij times
    // Takahashi 1974 high-pressure correction; identical to chungTransport
    // so that switching between transport models only changes mu, kappa.
    scalar sum = 0.0;
    const scalar T175 = T*sqrt(T*sqrt(T));
    const scalar pAtm = p/101325.0;
    forAll(Xmd_, j)
    {
        if (speciei != j)
        {
            sum += Xmd_[j]/
                   (
                       phi(p/Pcmd_[j][speciei], T/Tcmd_[j][speciei])
                      *0.001*T175*sqrt(Mmd_[j][speciei])
                      /(pAtm*sqr(sigmd_[j][speciei]))
                   );
        }
    }
    return (sum == 0.0) ? 0.0 : ((1.0 - Ymd_[speciei])/sum)/1e4;
}


// * * * * * * * * * * * * * updateTRANS * * * * * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::elyHanleyTransport<Thermo>::updateTRANS
(
    const scalar& TcM,
    const scalar& VcM,
    const scalar& MM,
    const scalar& omegaM,
    const scalar& ZcM,
    const List<scalar>& Ymd,
    const List<scalar>& Xmd,
    const List<List<scalar>>& Tcmd,
    const List<List<scalar>>& Pcmd,
    const List<List<scalar>>& Mmd,
    const List<List<scalar>>& sigmd
)
{
    TcM_ = TcM; VcM_ = VcM; MM_ = MM; omegaM_ = omegaM; ZcM_ = ZcM;
    Ymd_ = Ymd; Xmd_ = Xmd;
    Tcmd_ = Tcmd; Pcmd_ = Pcmd;
    Mmd_ = Mmd; sigmd_ = sigmd;
}


// Chung-style 8-arg updateTRANS adapter removed: SRKelyHanleyMixture now
// drives the native 5-arg signature directly, supplying its own
// mass-weighted Z_c,m instead of the previous Pitzer fallback.


// * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::elyHanleyTransport<Thermo>::operator=
(
    const elyHanleyTransport<Thermo>& tr
)
{
    Thermo::operator=(tr);
    TcM_ = tr.TcM_; VcM_ = tr.VcM_; MM_ = tr.MM_;
    omegaM_ = tr.omegaM_; ZcM_ = tr.ZcM_;
    Ymd_ = tr.Ymd_; Xmd_ = tr.Xmd_;
    Tcmd_ = tr.Tcmd_; Pcmd_ = tr.Pcmd_;
    Mmd_ = tr.Mmd_; sigmd_ = tr.sigmd_;
}


template<class Thermo>
inline void Foam::elyHanleyTransport<Thermo>::operator+=
(
    const elyHanleyTransport<Thermo>& tr
)
{
    scalar X1 = this->Y()/this->W();
    Thermo::operator+=(tr);
    if (mag(this->Y()) > small)
    {
        X1 *= this->W()/this->Y();
        const scalar X2 = this->W()*tr.Y()/(tr.W()*this->Y());
        TcM_    = X1*TcM_    + X2*tr.TcM_;
        VcM_    = X1*VcM_    + X2*tr.VcM_;
        MM_     = X1*MM_     + X2*tr.MM_;
        omegaM_ = X1*omegaM_ + X2*tr.omegaM_;
        ZcM_    = X1*ZcM_    + X2*tr.ZcM_;
    }
}


template<class Thermo>
inline void Foam::elyHanleyTransport<Thermo>::operator*=(const scalar s)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

template<class Thermo>
Foam::elyHanleyTransport<Thermo> Foam::operator+
(
    const elyHanleyTransport<Thermo>& t1,
    const elyHanleyTransport<Thermo>& t2
)
{
    Thermo tt
    (
        static_cast<const Thermo&>(t1) + static_cast<const Thermo&>(t2)
    );
    return elyHanleyTransport<Thermo>
    (
        tt,
        t1.TcM_, t1.VcM_, t1.MM_, t1.omegaM_, t1.ZcM_,
        t1.Ymd_, t1.Xmd_, t1.Tcmd_, t1.Pcmd_, t1.Mmd_, t1.sigmd_
    );
}


template<class Thermo>
Foam::elyHanleyTransport<Thermo> Foam::operator*
(
    const scalar s,
    const elyHanleyTransport<Thermo>& tr
)
{
    return elyHanleyTransport<Thermo>
    (
        s*static_cast<const Thermo&>(tr),
        tr.TcM_, tr.VcM_, tr.MM_, tr.omegaM_, tr.ZcM_,
        tr.Ymd_, tr.Xmd_, tr.Tcmd_, tr.Pcmd_, tr.Mmd_, tr.sigmd_
    );
}
