/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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
    Foam::solvers::peqsiFluid

Description
    Pressure-Equilibrium Semi-Implicit (PEQSI) solver module for transient
    real-fluid flows, implementing the characteristic fractional-step scheme
    of

    \verbatim
        Wada, S., Kai, R., Pillai, A. L., Yamada, T. & Kurose, R. (2024).
        Pressure-equilibrium semi-implicit solver for real fluids.
        Physics of Fluids 36 (11), 116104.  (PEQSI)

        Wada, S., Kai, R. & Kurose, R. (2024).
        A new semi-implicit pressure-based solver considering real gas
        effect.  Journal of Computational Physics 501, 112782.  (WKK)
    \endverbatim

    Each time step consists of
      1. an explicit advective substep (SSP-RK3, non-conservative advective
         form) for rho*, (rho u)*, h* (WKK Eqs. 15-17) and p* (PEQSI Eq. 10,
         replacing WKK's EOS inversion for p*);
      2. an implicit acoustic substep: the modified Helmholtz equation for
         dp = p^{n+1} - p* in the pressure-equilibrium-consistent form
         (PEQSI Eq. 19), solved with a stabilised bi-conjugate-gradient
         method;
      3. the coupled updates of rho, rho u, rho h and p with dp
         (WKK Eqs. 22-24, PEQSI Eqs. 11-13).

    The real-gas coefficients (WKK Eqs. 9, 11, 12)
        xi    = (dh/dT)_v
        alpha = (1/(rho xi)) (dp/dT)_v
        beta  = (1/rho) [ (dp/dv)_T - (1/xi)(dp/dT)_v (dh/dv)_T ]
    satisfy beta/(1-alpha) == -rho c^2 (WKK App. D), which is used as a
    run-time consistency check.

    Mass is conserved by construction (the continuity equation closes the
    Helmholtz derivation, WKK Eq. 27); the conservation defect of the family
    is parked in the energy equation (< 5% in the reference publications)
    and is reported by the run-time audit.

    Transient (Euler fractional-step) only: no LTS, no mesh motion.

SourceFiles
    peqsiFluid.C
    advectSubstep.C
    correctPressurePEQSI.C
    thermoClosure.C

See also
    Foam::solvers::isothermalFluid
    Foam::solvers::fgmFluid

\*---------------------------------------------------------------------------*/

#ifndef peqsiFluid_H
#define peqsiFluid_H

#include "isothermalFluid.H"
#include "fluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class peqsiFluid Declaration
\*---------------------------------------------------------------------------*/

class peqsiFluid
:
    public isothermalFluid
{

protected:

    // Thermophysical properties

        //- Multicomponent thermo (SRK real-fluid instantiations from
        //  libRGP13realFluid; single-species N2 for the inert validation)
        fluidMulticomponentThermo& thermo_;


    // Transported fields

        //- Specific total enthalpy [J/kg] (transported; the thermo he()
        //  field is NOT used as storage so that no thermo call can
        //  silently overwrite the transported value)
        volScalarField h_;


    // PEQSI real-gas coefficient fields (evaluated in thermoClosure)

        //- alpha = (1/(rho xi)) (dp/dT)_v   [-]
        volScalarField alpha_;

        //- beta = (1/rho)[(dp/dv)_T - (1/xi)(dp/dT)_v (dh/dv)_T]  [m2/s2]
        //  Note beta < 0; beta/(1-alpha) == -rho c^2
        volScalarField beta_;


    // Substep state (set by advectSubstep, consumed by
    // correctPressurePEQSI, cleared afterwards)

        //- rho at step start (n-state)
        autoPtr<volScalarField> rhoN_;

        //- U at step start (n-state)
        autoPtr<volVectorField> UN_;

        //- p at step start (n-state)
        autoPtr<volScalarField> pN_;

        //- h at step start (n-state)
        autoPtr<volScalarField> hN_;

        //- The advective substep's own compression bookkeeping
        //  sComp = (1/6 rho^n + 1/6 rho^(1) + 2/3 rho^(2)) div(phiv):
        //  the SSP-RK3 effective quadrature of rho div(u).  Using this as
        //  the Helmholtz source makes int(coef dp) cancel the substep's
        //  mass change identically (see correctPressurePEQSI).
        autoPtr<volScalarField> sComp_;

        //- Guard so the acoustic substep runs once per time step
        label acousticTimeIndex_;

        //- Previous-step cp (saved by updateCoefficients; used by the
        //  LAD explicit-stability cap)
        mutable scalarField cpPrev_;

        //- Cached SGS-active flag (simulationType != laminar),
        //  resolved once
        mutable label sgsActive_ = -1;

        //- Cached LAD directional geometry (static mesh): direction
        //  face weights (n.e_l)^2 and per-cell directional spacings
        mutable autoPtr<PtrList<surfaceScalarField>> ladWDir_;
        mutable autoPtr<PtrList<volScalarField>> ladDeltaDir_;

        //- Per-cell smallest INTERNAL-face spacing (great where a
        //  direction has no internal faces, e.g. the empty direction):
        //  sets the LAD explicit-stability cap
        mutable autoPtr<scalarField> ladDeltaMin_;

        //- Per-phase cumulative CPU-time accounting (peqsiTimers switch):
        //  setup, sgs, lad, rk3, helmAsm, helmSolve, updates, newton,
        //  coeffs
        mutable scalar tPhase_[9] = {0,0,0,0,0,0,0,0,0};
        mutable label tPhaseSteps_ = 0;

        //- Initial domain mass, set on the first acoustic substep;
        //  reference for the cumulative conservation audit (the paper's
        //  metric: relative drift of int(rho) dV)
        scalar initialMass_;

        //- Initial domain integral of rho*h (energy-audit reference)
        scalar initialRhoH_;


    // Protected Member Functions

        //- Build (once) the cached directional geometry shared by the
        //  LAD and the explicit filter: face weights (n.e_l)^2 and
        //  per-cell directional spacings
        void ensureDirGeometry() const;

        //- Explicit 8th-order selective low-pass filter (TK 2012
        //  Sec. 2.7 family; Visbal-Gaitonde explicit variant) applied
        //  once per direction to p and U after the acoustic update:
        //  q -= sigma/256 (delta_l^2 d2/dx_l^2)^4 q.  Nyquist gain
        //  1-sigma, smooth modes ~k^8 untouched; exact identity on a
        //  uniform field (pressure equilibrium preserved).  rho and
        //  rho*h are NOT filtered: the mass and energy ledgers stay
        //  untouched (documented deviation from TK's variable set).
        void applyFilter(const scalar sigma, const bool filterRho);

        //- Evaluate alpha_, beta_ from the current (T, v, Y) state and
        //  verify beta/(1-alpha) against -rho c^2
        void updateCoefficients();

        //- Invert h(T, v, Y) for T at fixed specific volume with Newton
        //  iterations using xi = (dh/dT)_v (WKK Fig. 3; NOT cp), then
        //  refresh transport properties at the new state
        void invertTemperature();


public:

    //- Runtime type information
    TypeName("peqsiFluid");


    // Constructors

        //- Construct from region mesh
        peqsiFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        peqsiFluid(const peqsiFluid&) = delete;


    //- Destructor
    virtual ~peqsiFluid();


    // Member Functions

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- No-op override: the base class solves a continuity predictor
        //  for rho here, which would advance the density outside the
        //  fractional step and corrupt the scheme's mass bookkeeping
        virtual void prePredictor();

        //- Advective substep: SSP-RK3 explicit update of rho*, (rho u)*,
        //  p*, h* with the frozen n-state velocity (characteristic split)
        virtual void momentumPredictor();

        //- No-op: the thermodynamic closure runs after the acoustic
        //  substep (see correctPressurePEQSI)
        virtual void thermophysicalPredictor();

        //- Acoustic substep: modified Helmholtz for dp (PEQSI Eq. 19),
        //  coupled updates (WKK Eqs. 22-24) and thermodynamic closure
        virtual void pressureCorrector();

        //- Static-mesh no-op (mesh motion is not supported)
        virtual void motionCorrector();

        //- No-op override: the base post-corrector re-syncs state the
        //  fractional step already finalised.  TODO(V3): call the
        //  momentum transport correct() here when LES is enabled.
        virtual void postCorrector();

        //- Override WITHOUT the base's rho_ = thermo.rho() re-sync: the
        //  transported (continuity) density is the scheme's mass ledger
        //  and must never be replaced by the EOS density (that stomp was
        //  measured as +0.8% mass/step on the 1-D case A); their
        //  difference is the family's energy-side drift, reported by the
        //  thermo closure
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const peqsiFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
