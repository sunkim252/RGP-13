/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2024-2025 RGP-13
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

Class
    Foam::elyHanleyTransport

Description
    Real-gas transport model for viscosity and thermal conductivity using
    Ely & Hanley's extended corresponding states (ECS) framework with
    methane as the reference fluid.

    Theory:
      - Ely, J.F. & Hanley, H.J.M., "Prediction of Transport Properties.
        1. Viscosity of Fluids and Mixtures", Ind. Eng. Chem. Fundam.
        20 (1981) 323-332.
      - Ely, J.F. & Hanley, H.J.M., "Prediction of Transport Properties.
        2. Thermal Conductivity of Pure Fluids and Mixtures", Ind. Eng.
        Chem. Fundam. 22 (1981) 90-97.
      - Hanley, H.J.M., McCarty, R.D. & Haynes, W.M., "The viscosity and
        thermal conductivity coefficients of dilute and dense fluid
        methane", J. Phys. Chem. Ref. Data 4 (1975) 135.
      - Leach, J.W., Chappelear, P.S. & Leland, T.W., "Use of molecular
        shape factors in vapor-liquid equilibrium calculations with the
        corresponding states principle", AIChE J. 14 (1968) 568.

    Core equations (Oefelein AIAA 2019-0634 Eqs. 39-54):
      mu(rho, T)     = mu_o(rho_o, T_o) * F_mu * X_eta
      lambda(rho, T) = lam_o(rho_o, T_o) * F_lam * X_lam
    where
      rho_o = rho * h_x,  T_o = T / f_x
      F_mu  = (Wx/Wo)^0.5 * f_x^0.5 * h_x^(-2/3)
      F_lam = (Wo/Wx)^0.5 * f_x^0.5 * h_x^(-2/3)
      h_x, f_x built from Leach-Chappelear-Leland shape factors.

    Mass diffusivity is left to the Fuller + Takahashi pipeline already
    present in the repository (shared with chungTransport) so that only
    mu and kappa change when this transport model is selected.

    Drop-in replacement for chungTransport in the same template slot.

SourceFiles
    elyHanleyTransportI.H
    elyHanleyTransport.C

\*---------------------------------------------------------------------------*/

#ifndef elyHanleyTransport_H
#define elyHanleyTransport_H

#include "PtrList.H"

namespace Foam
{

template<class Thermo> class elyHanleyTransport;

template<class Thermo>
inline elyHanleyTransport<Thermo> operator+
(
    const elyHanleyTransport<Thermo>&,
    const elyHanleyTransport<Thermo>&
);

template<class Thermo>
inline elyHanleyTransport<Thermo> operator*
(
    const scalar,
    const elyHanleyTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const elyHanleyTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                      Class elyHanleyTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class elyHanleyTransport
:
    public Thermo
{

protected:

    // Protected data (mixture-averaged quantities cached per cell)

        //- Mixture critical temperature [K]
        scalar TcM_;

        //- Mixture critical volume [m^3/kmol]
        scalar VcM_;

        //- Mixture molecular weight [kg/kmol]
        scalar MM_;

        //- Mixture acentric factor [-]
        scalar omegaM_;

        //- Mixture critical compressibility factor Zc = pc Vc / (R Tc) [-]
        scalar ZcM_;

        //- Reference fluid (methane) critical temperature [K]
        scalar TcO_;

        //- Reference fluid (methane) critical molar volume [m^3/kmol]
        scalar VcO_;

        //- Reference fluid (methane) molecular weight [kg/kmol]
        scalar MO_;

        //- Reference fluid (methane) acentric factor [-]
        scalar omegaO_;

        //- Reference fluid (methane) critical compressibility [-]
        scalar ZcO_;

        //- Chung auxiliary parameters stored for the non-polar fallback
        //  branch (see mu(), kappa() below).
        scalar sigmaM_;      // effective LJ diameter [Angstrom]
        scalar epsilonkM_;   // effective LJ epsilon/k [K]
        scalar miuiM_;       // reduced dipole moment
        scalar kappaiM_;     // association factor

        //- Members for diffusivity evaluation (Fuller + Takahashi).
        //  Kept for API compatibility with chungTransport.

            List<scalar> Ymd_;
            List<scalar> Xmd_;
            List<List<scalar>> Tcmd_;
            List<List<scalar>> Pcmd_;
            List<List<scalar>> Mmd_;
            List<List<scalar>> sigmd_;


    // Protected Member Functions

        //- Leach-Chappelear-Leland shape factor phi(V_r, T_r, omega).
        //  Returns the factor such that h_i = (Vc_i/Vc_o) * phi.
        inline scalar phiShape
        (
            const scalar Vri,
            const scalar Tri,
            const scalar omegai
        ) const;

        //- Leach-Chappelear-Leland shape factor theta(V_r, T_r, omega).
        //  Returns theta such that f_i = (Tc_i/Tc_o) * theta.
        inline scalar thetaShape
        (
            const scalar Vri,
            const scalar Tri,
            const scalar omegai
        ) const;

        //- Methane reference viscosity mu_o(rho_o, T_o) [Pa.s].
        //  Three-term decomposition per Hanley-McCarty-Haynes 1975:
        //    mu_o = mu_o*(T_o) + mu_o**(T_o) * rho_o + Delta_mu_o(rho_o, T_o)
        inline scalar muRef
        (
            const scalar rhoO,
            const scalar TO
        ) const;

        //- Methane reference thermal conductivity [W/m/K].
        inline scalar lambdaRef
        (
            const scalar rhoO,
            const scalar TO
        ) const;

        //- Enskog size-difference correction X_eta (Ely 1981 Eq. 53)
        inline scalar enskogXeta
        (
            const scalar Tx,
            const scalar fx,
            const scalar dfxdT
        ) const;

        //- Enskog size-difference correction X_lam (Ely 1981 Eq. 54)
        inline scalar enskogXlam
        (
            const scalar Tx,
            const scalar fx,
            const scalar dfxdT
        ) const;

        //- Chung-Lee-Starling viscosity (fallback for non-polar species).
        //  Self-contained copy of the Chung formula using the stored
        //  sigmaM, epsilonkM, omegaM, miuiM, kappaiM coefficients.
        inline scalar muChung
        (
            const scalar p,
            const scalar T
        ) const;

        //- Chung-Lee-Starling thermal conductivity fallback.
        inline scalar kappaChung
        (
            const scalar p,
            const scalar T
        ) const;

        //- Polarity test. Returns true if ECS/Ely-Hanley should be used;
        //  false to fall back to Chung-Lee-Starling.
        //
        //  Criterion (tightened for multi-component mixtures):
        //      |omega_mix| >= 0.30  AND  miuiM_mix > 0
        //
        //  Rationale based on measured mixture RMSE:
        //    - Pure H2O (omega=0.344, miui=1.844):       use Ely-Hanley
        //    - Pure O2  (omega=0.025):                    use Chung
        //    - Pure H2  (omega=-0.217, miui=0):           use Chung
        //    - O2/H2O mix with Y_H2O ~ 0.5 (omega_mix ~ 0.23):
        //          use Chung (mid-comp ECS causes mu regression)
        //    - O2/H2O mix with Y_H2O ~ 0.9 (omega_mix ~ 0.32):
        //          use Ely-Hanley (delivers -100 pp kappa improvement)
        //
        //  The 0.30 threshold isolates the regime where the mixture
        //  really is dominated by a strongly polar species.
        inline bool useElyHanleyBranch() const
        {
            return (mag(omegaM_) >= 0.30) && (miuiM_ > 1e-12);
        }


public:

    // Constructors

        inline elyHanleyTransport
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
        );

        inline elyHanleyTransport(const word&, const elyHanleyTransport&);

        elyHanleyTransport(const word& name, const dictionary& dict);

        elyHanleyTransport(const Thermo& t, const dictionary& dict);

        inline autoPtr<elyHanleyTransport> clone() const;


    // Member Functions

        static word typeName()
        {
            return "elyHanley<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Takahashi correction coefficient for binary diffusion.
        //  Implemented identically to chungTransport::phi(PRR, TRR) so that
        //  Dimix() remains drop-in compatible.
        inline scalar phi(scalar PRR, scalar TRR) const;

        //- Mass diffusivity of species i (Fuller + Takahashi) [m^2/s]
        inline scalar Dimix
        (
            label speciei,
            const scalar p,
            const scalar T
        ) const;

        //- Update mixture-level quantities (native Ely-Hanley signature)
        inline void updateTRANS
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
        );

        //- Update with Chung-style signature so that SRKchungTakaMixture
        //  can drive either transport backend without changes.
        //  sigmaM, epsilonkM, miuiM, kappaiM are accepted but unused.
        //  ZcM is estimated via Pitzer correlation Zc ~ 0.291 - 0.08*omega.
        inline void updateTRANS
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
        );

        //- Access functions
        inline scalar TcM()    const { return TcM_; }
        inline scalar VcM()    const { return VcM_; }
        inline scalar MM()     const { return MM_;  }
        inline scalar omegaM() const { return omegaM_; }
        inline scalar ZcM()    const { return ZcM_; }
        inline List<scalar> Ymd() const { return Ymd_; }
        inline List<scalar> Xmd() const { return Xmd_; }
        inline List<List<scalar>> Tcmd() const { return Tcmd_; }
        inline List<List<scalar>> Pcmd() const { return Pcmd_; }
        inline List<List<scalar>> Mmd()  const { return Mmd_;  }
        inline List<List<scalar>> sigmd() const { return sigmd_; }

        void write(Ostream& os) const;


    // Member Operators

        inline void operator=(const elyHanleyTransport&);
        inline void operator+=(const elyHanleyTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend elyHanleyTransport operator+ <Thermo>
        (
            const elyHanleyTransport&,
            const elyHanleyTransport&
        );

        friend elyHanleyTransport operator* <Thermo>
        (
            const scalar,
            const elyHanleyTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const elyHanleyTransport&
        );
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "elyHanleyTransportI.H"

#ifdef NoRepository
    #include "elyHanleyTransport.C"
#endif

#endif
