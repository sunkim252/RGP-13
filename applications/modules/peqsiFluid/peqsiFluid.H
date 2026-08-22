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
#include "FGMTable.H"

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


    // FGM (flamelet-manifold) coupling -- OFF unless peqsiFGM is set.
    //
    // Rationale and equations: PEQSI Sec. IV states the single-species
    // limitation and the extension route (h = h(rho, p, Y); the species
    // diffusion terms enter the pressure equation).  The multi-species
    // pressure-evolution equation itself is published in Terashima &
    // Koshi 2012 Sec. 2.4:
    //   Dp/Dt = (dp/drho)_{e,Y} Drho/Dt + (dp/de)_{rho,Y} De/Dt
    //         + sum_i (dp/dY_i)_{rho,e} DY_i/Dt
    // and the manifold closure is the project's existing validated
    // FGM/FPV table (fgmFluid, FGMTable).  Nothing here is a new
    // scheme: PEQSI's semi-implicit split + TK's multi-species pressure
    // evolution + the existing manifold.
    //
    // Stage 1 (this commit) declares the state only; the transport,
    // manifold closure and the composition source term follow in the
    // staged work (wiki "PEQSI Method and Discretisation" Sec. 12).

        //- Master switch (PIMPLE dict key peqsiFGM, default off).  With
        //  it off every field below stays unallocated and the solver is
        //  bit-identical to the validated inert path.
        Switch fgmActive_;

        //- Transported mixture fraction (mean)
        autoPtr<volScalarField> Z_;

        //- Transported mixture-fraction variance
        autoPtr<volScalarField> Zvar_;

        //- Transported progress variable (unnormalised Yc, as in
        //  fgmFluid: removes the Yc_eq cross terms by construction)
        autoPtr<volScalarField> Yc_;

        //- Manifold source of the progress variable (table lookup)
        autoPtr<volScalarField> sourceYc_;

        //- The manifold (stage 2a): constructed only when peqsiFGM is
        //  set, from constant/fgmProperties -- the baked table with the
        //  PEQSI coefficient blocks and W (verified against md5 on
        //  delivery; the reader FatalErrors on any missing block).
        autoPtr<FGMTable> fgmTable_;

        //- n-state copies of the transported manifold coordinates for
        //  the SSP-RK3 combination (held per advective substep, like
        //  rhoN_/UN_)
        autoPtr<volScalarField> ZN_;
        autoPtr<volScalarField> ZvarN_;
        autoPtr<volScalarField> YcN_;

        //- Manifold-supplied real-gas coefficients and the cells that
        //  may use them (stage 2a+).  fgmClosure fills these from the
        //  table with the phase rule applied; updateCoefficients uses
        //  them in place of the runtime SRK derivatives where the mask
        //  allows, and falls back to the runtime evaluation elsewhere.
        autoPtr<scalarField> alphaTab_;
        autoPtr<scalarField> betaTab_;
        autoPtr<boolList> tabUsable_;

        //- Manifold temperature guess for the Newton inversion, with
        //  the table's own (dT/dv)_h correcting it to the cell's
        //  specific volume.  Written by fgmClosure, consumed by
        //  invertTemperature.
        autoPtr<scalarField> Tguess_;

        //- Opt-1 stencil coordinates.  The mixture locates the 16 corner
        //  node mixtures from these, so they must be the SAME values
        //  fgmClosure feeds makeStencil: c normalised by Cnorm(Z), and
        //  dh measured from the mixing line and the dhRef slice.
        autoPtr<scalarField> cNormF_;
        autoPtr<scalarField> dhF_;

        //- Opt-1 armed (base-thermo blend read off the manifold nodes
        //  instead of summed over the mixture species)
        bool baseBlendArmed_ = false;

        //- Algebraic mixture-fraction variance closure.
        //  Zvar was being ADVECTED with no source at all -- neither of
        //  the two published families, just an incomplete transport
        //  equation whose segregation axis never got fed.  This is the
        //  gradient / scale-similarity closure
        //      Z''^2 = Cv L^2 |grad Z|^2,   gZ = Z''^2 / (Z (1 - Z))
        //  (Cook & Riley, Phys. Fluids 6:2868, 1994; filter-width
        //  mixing length after Pierce & Moin, JFM 504:73, 2004), which
        //  is what fgmFluid already runs in production.
        enum class varModel { les, laminar };
        varModel varModel_ = varModel::laminar;
        scalar Cv_ = 0.1;
        bool varModelReported_ = false;

        //- Tier-2 real-gas coefficients per cell, in coeffNames() order.
        //  Opt-1 removes the base-thermo species sum; this removes the
        //  O(n^2) calculateRealGas pair sum behind it (11236 pairs per
        //  cell at 106 species).  Filled by fgmClosure on the stencil it
        //  already builds; read by the mixture in place of the sum.
        List<scalarField> RGfields_;
        bool tier2Armed_ = false;

        //- Mesh state the manifold lookups were armed for.  Opt-1 and
        //  Tier-2 both recover a cell index from the ADDRESS RANGE of
        //  cached field pointers, so a redistribution (DLB) or a
        //  refinement reallocates those fields and the recovery then
        //  yields a wrong IN-RANGE index -- an out-of-bounds read, not a
        //  clean failure.  The cell count alone misses an equal-count
        //  redistribution, which moves the data pointer without changing
        //  the size, so the species-0 data pointer is tracked too.
        //- FPV composition source for the pressure substep, S_Y.
        //  alpha and beta are derived at FIXED composition (WKK Eq. 9,
        //  11, 12), so nothing in the pressure equation knows that the
        //  manifold moves Y.  Measured consequence: heat release produces
        //  no dilatation at all.  Filled by updateCompositionSource(),
        //  consumed by the p* stage.
        autoPtr<volScalarField> sourceP_;

        //- (drho/dc) and (dh/dc) at fixed (p, T), differenced along the
        //  MANIFOLD rather than in time -- an Eulerian time difference
        //  would carry u.grad(Y) into what must be a material rate.
        autoPtr<scalarField> dYdcRho_;
        autoPtr<scalarField> dYdcH_;

        //- Cnorm(Z), the progress-variable normalisation.  Distinct from
        //  cNormF_, which holds the already-normalised c -- S_Y needs the
        //  normalisation itself to turn dYc/dt into Dc/Dt.
        autoPtr<scalarField> cnormZ_;

        //- Table density at the base state, filled on the stencil
        //  fgmClosure already builds
        autoPtr<scalarField> rhoTabF_;

        //- Reference internal-energy field for peqsiEnergyAudit
        autoPtr<scalarField> reInit_;
        scalar E0_ = 0;

        label armedNCells_ = -1;
        const scalar* armedYData_ = nullptr;


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

        //- LAD explicit-diffusion time-step limit, refreshed by the
        //  advective substep from the UNCAPPED artificial-diffusivity
        //  demand.  The 141520 post-mortem: the demand is set by the
        //  flow (1.15e-3 m^2/s at the jet-head interface) while the
        //  explicit stability cap is 0.45 d^2/dt -- raise dt past
        //  0.45 d^2/D and the cap clips the dissipation the interface
        //  needs; dp then grew 2.3-7x per step and the run died twice
        //  at the same dt.  This limit keeps dt below that crossing
        //  instead of finding it by crashing.
        scalar ladDtLimit_;

        // Single-species SRK replica for the constant-v closure: only
        // p(v, T) is replicated (15 lines of cubic algebra, molar);
        // enthalpy comes from the library's he(p, T).  Cross-checked
        // against the library state on first use.  The FGM stage
        // replaces this with manifold-baked coefficients.
        bool srkReplicaValid_;
        mutable bool srkChecked_;
        scalar srkW_;
        scalar srkB_;
        scalar srkCoef1_;
        scalar srkCoef2_;
        scalar srkCoef3_;
        scalar srkC_;
        scalar srkCq0_;
        scalar srkCq1_;
        scalar srkCq2_;
        scalar srkCTlo_;
        scalar srkCThi_;

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

        //- Sensor-gated conservative 2nd-order shock-capturing filter
        //  (Bogey, de Cacqueray & Bailly, JCP 228 (2009) 1447)
        void applySCFilter(const scalar sigmaMax, const scalar rTh);

        //- Stage 2a: one SSP-RK3 stage of the manifold coordinates
        //  (advective, non-conservative form -- composition must not
        //  respond to the acoustic projection at a contact)
        void advectManifoldStage
        (
            const surfaceScalarField& phiv,
            const volScalarField& divPhiv,
            const scalar cOld,
            const scalar cNew,
            const scalar dtv
        );

        //- Stage 2a: manifold lookup -> composition + diagnostics
        void fgmClosure();

        //- Fill the segregation factor gZ (stored in Zvar_) from the
        //  algebraic variance closure
        void updateSegregation();

        //- Assemble S_Y from the manifold-direction composition
        //  derivatives and the reaction rate
        void updateCompositionSource();

        //- Bounding artificial diffusivity for a [0,1] scalar, after
        //  Terashima & Koshi 2012 Eq. (27):
        //
        //    D* = C_Y c bar[(Y-1)H(Y-1) - Y(1-H(Y))] Delta_Y
        //
        //  with C_Y = 100 (their standard value), c the local sound
        //  speed, Delta_Y = (dx dy dz)^(1/3) and the overbar the same
        //  truncated-Gaussian filter used by the LAD terms.  The bracket
        //  vanishes inside [0,1], so the term is EXACTLY zero wherever
        //  the scalar is bounded and grows with the violation -- it
        //  removes over/undershoots (and with them the temperature
        //  oscillations they drive) without touching resolved fronts.
        //  Used for the FGM manifold coordinates (Z, normalised c),
        //  which are bounded by construction.
        tmp<volScalarField> boundingArtDiffusivity
        (
            const volScalarField& Y,
            const volScalarField& c,
            const scalar CY
        ) const;

        //- Face mask: 0 on every face of a boundary-adjacent cell
        //  (filter identity in the boundary cell layer), exchanged
        //  across coupled patches
        tmp<surfaceScalarField> filterBoundaryMask() const;

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

        //- Return the maximum time-step: the base Courant limit
        //  combined with the LAD explicit-diffusion limit above
        virtual scalar maxDeltaT() const;

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
