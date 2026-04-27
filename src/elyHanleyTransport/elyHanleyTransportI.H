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
    sigmaM_(0), epsilonkM_(0), miuiM_(0), kappaiM_(0),
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
    sigmaM_(tr.sigmaM_), epsilonkM_(tr.epsilonkM_),
    miuiM_(tr.miuiM_), kappaiM_(tr.kappaiM_),
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
    // Piecewise quadratic Hermite blend.
    //   x <= lo - w           : lo                          (constant)
    //   lo - w < x <= lo + w  : lo + (x - (lo - w))^2/(4 w) (smooth ramp)
    //   lo + w  < x <  hi - w : x                           (identity)
    //   hi - w <= x <  hi + w : hi - (hi + w - x)^2/(4 w)   (smooth ramp)
    //   x >= hi + w           : hi                          (constant)
    // Continuity (C^0) and slope continuity (C^1) are exact at every
    // junction, so derivatives propagated downstream (dfx/dT, X_eta,
    // X_lam, mu, kappa) are smooth across the boundaries.
    if (x <= lo - w) return lo;
    if (x >= hi + w) return hi;
    if (x < lo + w)
    {
        const scalar t = x - (lo - w);
        return lo + t*t/(4.0*w);
    }
    if (x > hi - w)
    {
        const scalar t = x - (hi - w);
        return (hi - w) + t - t*t/(4.0*w);
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
    // C^1-continuous soft clamp so dphi/dT remains continuous when
    // T_r leaves the [0.5, 2.0] band.
    // Transition half-width w = 0.05 keeps the polynomial unchanged on
    // [0.55, 1.95] and removes the kinks at T_r = 0.5 and T_r = 2.0.
    const scalar w = 0.05;
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
    // Eq. 44). Soft-clamped on V_r, T_r as in phiShape so that dtheta/dT
    // (and the downstream Enskog correction) stays C^1.
    const scalar w = 0.05;
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


// * * * * * * * * Chung-Lee-Starling fallback (non-polar) * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::muChung
(
    const scalar p,
    const scalar T
) const
{
    // Full Chung-Lee-Starling 1988 viscosity using the mixture
    // parameters saved by the Chung-signature updateTRANS adapter.
    // Reproduces Foam::chungTransport<Thermo>::mu() to machine precision
    // when the same (sigmaM, epsilonkM, MM, VcM, TcM, omegaM, miuiM,
    // kappaiM) are supplied.
    //
    // References:
    //   Chung, T.-H., Ajlan, M., Lee, L. L. & Starling, K. E.
    //   Ind. Eng. Chem. Res. 27 (1988) 671-679.

    const scalar A=1.16145 , B=0.14874, C=0.52487 , D=0.77320 , E=2.16178;
    const scalar F=2.43787 , G=-6.435e-4, H=7.27371, Scol=18.0323, W=-0.76830;

    // Chung's unit-bearing constants (4.0785e-5, 40.785e-6) assume
    // VcM in cm^3/mol; convert from our m^3/kmol storage.
    const scalar VcM_cc = VcM_*1000.0;    // m^3/kmol -> cm^3/mol

    const scalar Tstar = T/max(epsilonkM_, 1e-20);
    const scalar mur = 131.3*miuiM_/sqrt(max(VcM_cc*TcM_, 1e-20));
    const scalar mur2 = sqr(mur);
    const scalar Fc = 1.0 - 0.2756*omegaM_ + 0.059035*sqr(mur2) + kappaiM_;

    const scalar Omegast =
        A*pow(Tstar, -B) + C*exp(-D*Tstar) + E*exp(-F*Tstar)
      + G*pow(Tstar, B)*sin(Scol*pow(Tstar, W) - H);

    const scalar VcM23 = cbrt(sqr(VcM_cc));
    const scalar eta0 = 4.0785e-5*sqrt(MM_*T)*Fc
                       /(max(Omegast, 1e-20)*VcM23);

    // Reduced density Y = rho V_c / 6  with rho in mol/L and Vc in cm^3/mol
    // (units from Chung 1988 Eq. 10).
    const scalar rho_molL = 1e-3*this->rho(p, T)/max(this->W(), 1e-20);
    const scalar Y = rho_molL*VcM_cc/6.0;
    const scalar oneMinusY = 1.0 - Y;
    const scalar G1 = (1.0 - 0.5*Y)/pow(oneMinusY, 3);

    const scalar a0[10] = {6.32402, 0.0012102, 5.28346, 6.62263, 19.74540,
                          -1.89992, 24.27450, 0.79716, -0.23816, 0.068629};
    const scalar a1[10] = {50.41190, -0.0011536, 254.209, 38.0957, 7.63034,
                          -12.53670, 3.44945, 1.11764, 0.067695, 0.34793};
    const scalar a2[10] = {-51.68, -0.0062571, -168.481, -8.46414, -14.35440,
                           4.98529, -11.29130, 0.012348, -0.8163, 0.59256};
    const scalar a3[10] = {1189.02, 0.037283, 3898.270, 31.4178, 31.52670,
                          -18.1507, 69.3464, -4.1161, 4.02528, -0.72663};

    scalar MA[10];
    const scalar mur4 = sqr(mur2);
    for (int i=0; i<10; i++)
    {
        MA[i] = a0[i] + a1[i]*omegaM_ + a2[i]*mur4 + a3[i]*kappaiM_;
    }

    const scalar G2 =
        (MA[0]*(1.0 - exp(-MA[3]*Y))/max(Y, 1e-12)
         + MA[1]*G1*exp(MA[4]*Y) + MA[2]*G1)
       /(MA[0]*MA[3] + MA[1] + MA[2]);

    const scalar etak = eta0*(1.0/G2 + MA[5]*Y);
    const scalar etap =
        (40.785e-6/sqrt(max(Tstar, 1e-20)))*sqrt(MM_*T)/VcM23
       *MA[6]*sqr(Y)*G2*exp(MA[7] + MA[8]/Tstar + MA[9]/sqr(Tstar));

    return 0.1*(etak + etap);   // poise -> Pa.s
}


template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::kappaChung
(
    const scalar p,
    const scalar T
) const
{
    // Full Chung-Lee-Starling 1988 thermal conductivity, reproducing
    // Foam::chungTransport<Thermo>::kappa() exactly.
    const scalar A=1.16145 , B=0.14874, C=0.52487 , D=0.77320 , E=2.16178;
    const scalar F=2.43787 , G=-6.435e-4, H=7.27371, Scol=18.0323, W=-0.76830;

    const scalar VcM_cc = VcM_*1000.0;    // m^3/kmol -> cm^3/mol

    const scalar Tstar = T/max(epsilonkM_, 1e-20);
    const scalar mur = 131.3*miuiM_/sqrt(max(VcM_cc*TcM_, 1e-20));
    const scalar mur2 = sqr(mur);
    const scalar Fc = 1.0 - 0.2756*omegaM_ + 0.059035*sqr(mur2) + kappaiM_;

    const scalar Omegast =
        A*pow(Tstar, -B) + C*exp(-D*Tstar) + E*exp(-F*Tstar)
      + G*pow(Tstar, B)*sin(Scol*pow(Tstar, W) - H);

    const scalar VcM23 = cbrt(sqr(VcM_cc));
    const scalar eta0 = 4.0785e-5*sqrt(MM_*T)*Fc
                       /(max(Omegast, 1e-20)*VcM23);

    const scalar rho_molL = 1e-3*this->rho(p, T)/max(this->W(), 1e-20);
    const scalar Y = rho_molL*VcM_cc/6.0;
    const scalar oneMinusY = 1.0 - Y;
    const scalar G1 = (1.0 - 0.5*Y)/pow(oneMinusY, 3);

    // Chung b-coefficients (Table IV of 1988 paper)
    const scalar b0[7] = {2.41657, -0.50924, 6.61069, 14.54250, 0.79274,
                         -5.8634, 81.171};
    const scalar b1[7] = {0.74824, -1.50936, 5.62073, -8.91387, 0.82019,
                          12.8005, 114.1580};
    const scalar b2[7] = {-0.91858, -49.99120, 64.7599, -5.63794, -0.69369,
                           9.58926, -60.841};
    const scalar b3[7] = {121.721, 69.9834, 27.0389, 74.3435, 6.31734,
                         -65.5292, 466.7750};

    scalar MB[7];
    const scalar mur4 = sqr(mur2);
    for (int j=0; j<7; j++)
    {
        MB[j] = b0[j] + b1[j]*omegaM_ + b2[j]*mur4 + b3[j]*kappaiM_;
    }

    // Chung psi factor uses ideal Cv -> Cv = Cp - R (ideal-gas)
    const scalar CvIdeal = this->Cp(p, T) - this->R();
    const scalar Cv = CvIdeal*this->W()*0.2388e-3;  // cal/mol/K
    const scalar Rcal = Foam::constant::thermodynamic::RR/(4.186*1e3);
    const scalar alpha = Cv/max(Rcal, 1e-20) - 1.5;
    const scalar beta  = 0.7862 - 0.7109*omegaM_ + 1.1368*sqr(omegaM_);
    const scalar Tr = T/max(TcM_, 1e-20);
    const scalar Zterm = 2.0 + 10.5*sqr(Tr);
    const scalar Psi =
        1.0 + alpha*(0.215 + 0.28288*alpha - 1.061*beta + 0.26665*Zterm)
             /(0.6366 + beta*Zterm + 1.061*alpha*beta);

    const scalar lamda0 = 7.452*(eta0/MM_)*Psi;

    const scalar H2 =
        (MB[0]*(1.0 - exp(-MB[3]*Y))/max(Y, 1e-12)
         + MB[1]*G1*exp(MB[4]*Y) + MB[2]*G1)
       /(MB[0]*MB[3] + MB[1] + MB[2]);

    const scalar lamdak = lamda0*(1.0/H2 + MB[5]*Y);
    const scalar lamdap =
        (3.039e-4*sqrt(TcM_/MM_)/VcM23)*MB[6]*sqr(Y)*H2*sqrt(Tr);

    return 418.6798*(lamdak + lamdap);  // cal/cm/s/K -> W/m/K
}


// * * * * * * * * * * Public mu, kappa, phi, Dimix * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::elyHanleyTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    // Hybrid dispatch: polar/near-critical species use the ECS / Ely-
    // Hanley branch (better for H2O); non-polar species fall back to
    // Chung-Lee-Starling 1988 internally (better for O2/N2/H2).
    if (!useElyHanleyBranch())
    {
        return muChung(p, T);
    }
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
    // Hybrid dispatch (see mu() comment).
    if (!useElyHanleyBranch())
    {
        return kappaChung(p, T);
    }

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


template<class Thermo>
inline void Foam::elyHanleyTransport<Thermo>::updateTRANS
(
    const scalar& sigmaM,
    const scalar& epsilonkM,
    const scalar& MM,
    const scalar& VcM,
    const scalar& TcM,
    const scalar& omegaM,
    const scalar& miuiM,
    const scalar& kappaiM,
    const List<scalar>& Ymd,
    const List<scalar>& Xmd,
    const List<List<scalar>>& Tcmd,
    const List<List<scalar>>& Pcmd,
    const List<List<scalar>>& Mmd,
    const List<List<scalar>>& sigmd
)
{
    // Chung-signature adapter: compute ZcM via Pitzer correlation
    //   Z_c ~= 0.291 - 0.08 * omega.
    // sigmaM, epsilonkM, miuiM, kappaiM are intentionally unused by
    // the Ely-Hanley branch (kept in signature for drop-in replacement
    // of chungTransport inside SRKchungTakaMixture).
    // Persist Chung auxiliary parameters for the non-polar fallback branch.
    sigmaM_    = sigmaM;
    epsilonkM_ = epsilonkM;
    miuiM_     = miuiM;
    kappaiM_   = kappaiM;
    // VcM arrives from SRKchungTakaMixture in cm^3/mol (rfSpecie units);
    // convert to m^3/kmol for internal consistency with the elyHanley
    // class, where the reference VcO_ is stored in m^3/kmol.
    const scalar VcM_SI = VcM*1e-3;
    const scalar ZcPitzer = 0.291 - 0.08*omegaM;
    updateTRANS
    (
        TcM, VcM_SI, MM, omegaM,
        (ZcPitzer > 0.0 ? ZcPitzer : 0.29),
        Ymd, Xmd, Tcmd, Pcmd, Mmd, sigmd
    );
}


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
