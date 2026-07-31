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

Description
    fgmFluid pressure corrector -- overrides isothermalFluid::pressureCorrector
    so the pressure equation can be replaced by a pressure-equilibrium-
    preserving (pressure-evolution) formulation. STEP 1 (this file, first
    build): correctPressurePEP() is a FAITHFUL COPY of the base
    isothermalFluid::correctPressure (non-buoyant, SIMPLErho path), so the 1-D
    advecting-interface benchmark must reproduce the baseline spurious pressure
    spike unchanged -- confirming the override mechanism is isolated before the
    PEP modification is introduced. Reference for the PEP target:
    Terashima & Koshi, J. Comput. Phys. 231 (2012) 6907; Kai/Kurose PEQSI,
    Phys. Fluids 36 (2024) 116104.

\*---------------------------------------------------------------------------*/

#include "fgmFluid.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "fvcMeshPhi.H"
#include "fvcFlux.H"
#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "safeReconstruct.H"
#include "fvcSmooth.H"
#include "fvcVolumeIntegrate.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "fvcAverage.H"
#include "fvmDdt.H"
#include "fvmSup.H"
#include "upwind.H"
#include "zeroGradientFvPatchFields.H"
#include "tabulatedRealGasMixture.H"

namespace
{

// Root-selection hysteresis (2026-07-24, extended): for cells with a
// genuine 3-real-root EOS ambiguity, pick whichever ACTUAL cubic root
// (never an interpolated non-root value) gives a density closer to the
// cell's PREVIOUS-timestep density -- "sticky" branch selection driven by
// history rather than an instantaneous, near-the-crossing noise-sensitive
// fugacity comparison (see PEP and LAD Spike Suppression wiki secs 15-18).
//
// CRITICAL: rho and psi must be re-derived from the SAME chosen branch. An
// earlier version corrected rho alone and left psi (read live via
// thermo.psi() at the psis-construction call site) on the ORIGINAL hard-
// argmin branch -- rho from branch B, psi from branch A at the same cell,
// which corrupts the diagonal of the globally-coupled PEP pressure
// equation and was found to blow the domain pressure up to the pMaxPaDome
// clamp (wiki sec 18b). Passing a non-null psi here keeps both consistent.
Foam::label correctRootHysteresis
(
    const Foam::tabulatedRealGasMixture& hook,
    const Foam::PtrList<Foam::volScalarField>& Y,
    const Foam::volScalarField& p,
    const Foam::volScalarField& T,
    const Foam::scalarField& rhoOld,
    Foam::volScalarField& rho,
    Foam::volScalarField* psi,
    // Near-critical (T,p) band pre-filter (2026-07-24c): a genuine 3-real-
    // root ambiguity is only possible near the mixture's critical point --
    // cells outside this band are GUARANTEED D>0 (single real root) and
    // would always be rejected by the EOS evaluation below anyway, but
    // only after paying for it to find that out. Skips ~40-50% of the
    // domain (measured via the band_cells diagnostic scan). Defaults match
    // the established diagnostic near-critical scan band (O2
    // Tc=154.6K/Pc=50.4bar).
    const Foam::scalar TMin = 140,
    const Foam::scalar TMax = 200,
    const Foam::scalar pMin = 4.0e6,
    const Foam::scalar pMax = 5.5e6,
    // Tabulated-coefficient fast path (2026-07-24d): when non-null (and
    // sized for the 5 SRK coefficients bM/coef1-3/cM, coeffNames() order),
    // read bM/coef1-3/cM directly from this per-cell TABLE -- the SAME
    // table fgmFluid::RGcoeffFields_ that thermo_.correct() already
    // populated this timestep via enableCoeffTabulation() -- an O(1)
    // lookup, instead of recomputing them from composition via
    // calculateRealGas()'s O(n^2) species-pair mixing loop (~103 species
    // in thermo.wang2011). That redundant re-derivation, not the O(1)
    // psi/Cp/h work trimmed earlier, was the actual dominant cost (PEP and
    // LAD Spike Suppression wiki sec 18f). Null falls back to the live
    // composition-based eosRootBranchProperties() (e.g. tabulation
    // disabled/not yet armed).
    const Foam::PtrList<Foam::volScalarField>* RGcoeffFields = nullptr,
    // Mean molecular weight field, needed only by the tabulated path (the
    // live composition-based fallback derives Wm internally from Y). The
    // caller passes thermo_.W() -- computed once outside the per-cell loop
    // (it returns a tmp<volScalarField>, so calling it per-cell would
    // itself allocate a whole field just to read one value).
    const Foam::volScalarField* Wfield = nullptr,
    // Diagnostic only (2026-07-24e, drift-hypothesis follow-up): when
    // non-null, set to 1 for cells where the sticky pick DIFFERS from the
    // just-computed hard-argmin value, 0 otherwise (overwritten fresh each
    // call -- a snapshot of "which cells are hysteresis overriding RIGHT
    // NOW", not an accumulator). Lets a caller check whether overridden
    // cells sit at a density-jump boundary with their (uncorrected)
    // neighbours -- see PEP and LAD Spike Suppression wiki sec 20.
    Foam::volScalarField* flagField = nullptr,
    // Neighbour-consistency cap (2026-07-24f): the (T,p) band pre-filter
    // contains the domain-wide pressure blow-up (wiki sec 21) by letting a
    // blown-up cell fall out of hysteresis control and relax, but a
    // residual LOCAL spike survives -- confirmed to correlate with
    // corrected cells sitting at a density-jump boundary with their
    // uncorrected neighbours (corrected-uncorrected pairs show a ~5x
    // higher jump rate than the general population). Narrowing the (T,p)
    // band further is NOT viable for a combustion case (cells sweep a wide
    // T range continuously as they heat -- a narrow band only catches a
    // transient slice and reopens real bistable-pair coverage gaps; wiki
    // sec 21 follow-up). Instead: when rhoNeighbourAvg is non-null and
    // capRatio<GREAT, clip the STICKY pick itself to
    // [rhoNeighbourAvg/capRatio, rhoNeighbourAvg*capRatio] before it's
    // written -- mirrors the existing psisCapRatio pattern (face-averaged
    // neighbour value via fvc::average(fvc::interpolate(.))) but applied
    // to rho at the exact cells hysteresis touches, computed from the
    // PRE-correction field so it reflects genuine spatial consensus, not
    // hysteresis's own output.
    const Foam::volScalarField* rhoNeighbourAvg = nullptr,
    const Foam::scalar capRatio = Foam::GREAT
)
{
    const Foam::scalarField& pc = p.primitiveField();
    const Foam::scalarField& Tc = T.primitiveField();
    Foam::scalarField& rhoc = rho.primitiveFieldRef();
    Foam::scalarField* psic = psi ? &psi->primitiveFieldRef() : nullptr;
    Foam::scalarField* flagc =
        flagField ? &flagField->primitiveFieldRef() : nullptr;
    if (flagc)
    {
        *flagc = 0;
    }
    const bool tabulated =
        RGcoeffFields && RGcoeffFields->size() >= 5 && Wfield;

    Foam::List<Foam::scalar> Yc(Y.size());
    Foam::label nFlipped = 0;

    forAll(rhoc, celli)
    {
        if
        (
            Tc[celli] < TMin || Tc[celli] > TMax
         || pc[celli] < pMin || pc[celli] > pMax
        )
        {
            continue;
        }

        Foam::scalar rMin, rMax, psiMin, psiMax;
        bool ok;

        if (tabulated)
        {
            ok = hook.eosRootBranchPropertiesFromCoeffs
            (
                (*RGcoeffFields)[0].primitiveField()[celli],  // bM
                (*RGcoeffFields)[1].primitiveField()[celli],  // coef1
                (*RGcoeffFields)[2].primitiveField()[celli],  // coef2
                (*RGcoeffFields)[3].primitiveField()[celli],  // coef3
                (*RGcoeffFields)[4].primitiveField()[celli],  // cM
                Wfield->primitiveField()[celli],              // Wm
                pc[celli], Tc[celli],
                rMin, rMax, psiMin, psiMax,
                psic != nullptr
            );
        }
        else
        {
            forAll(Y, i)
            {
                Yc[i] = Y[i].primitiveField()[celli];
            }

            Foam::scalar CpMin, CpMax, hMin, hMax;
            ok = hook.eosRootBranchProperties
            (
                Yc, pc[celli], Tc[celli],
                rMin, rMax, psiMin, psiMax, CpMin, CpMax, hMin, hMax,
                psic != nullptr,  // needPsi: only the site-A (rho,psi) call
                false             // needCpH: never consumed in this module
            );
        }

        if (ok)
        {
            bool pickMin =
                Foam::mag(rMin - rhoOld[celli])
             <= Foam::mag(rMax - rhoOld[celli]);

            // Neighbour-consistency tie-break (2026-07-24f, corrected): if
            // the history-preferred branch is badly inconsistent with the
            // cell's (pre-correction) face-averaged neighbourhood, and the
            // OTHER branch would be more consistent, switch to it instead.
            // CRITICAL: this must pick between the two REAL cubic roots,
            // never clamp to an arbitrary in-between value -- an earlier
            // version clamped rHyst but left psi keyed to the ORIGINAL
            // pickMin, silently reintroducing the exact rho/psi branch
            // mismatch this whole mechanism exists to prevent (confirmed:
            // it neither closed the residual spike nor kept bistable pairs
            // at 0). Switching pickMin itself keeps rho AND psi consistent
            // with the SAME (possibly reconsidered) branch.
            if (rhoNeighbourAvg && capRatio < Foam::GREAT)
            {
                const Foam::scalar rAvg =
                    rhoNeighbourAvg->primitiveField()[celli];
                const Foam::scalar rPick = pickMin ? rMin : rMax;
                const Foam::scalar ratio =
                    Foam::max(rPick, rAvg)
                   /Foam::max(Foam::min(rPick, rAvg), Foam::small);
                if (ratio > capRatio)
                {
                    pickMin =
                        Foam::mag(rMin - rAvg) <= Foam::mag(rMax - rAvg);
                }
            }

            const Foam::scalar rHyst = pickMin ? rMin : rMax;

            // Relative tolerance, not ==: rhoc comes from the BLENDED
            // two-root EOS path (and may additionally have been relaxed), so
            // it is never bit-equal to a pure cubic root even where
            // hysteresis changed nothing. An exact test flagged every in-band
            // cell and made the diagnostic useless.
            if
            (
                Foam::mag(rHyst - rhoc[celli])
              > 1e-3*Foam::mag(rhoc[celli])
            )
            {
                nFlipped++;
                if (flagc)
                {
                    (*flagc)[celli] = 1;
                }
            }
            rhoc[celli] = rHyst;

            if (psic)
            {
                (*psic)[celli] = pickMin ? psiMin : psiMax;
            }
        }
    }

    rho.correctBoundaryConditions();
    if (psi)
    {
        psi->correctBoundaryConditions();
    }
    if (flagField)
    {
        flagField->correctBoundaryConditions();
    }

    // Global count: the callers report this through Info<<, which prints on
    // the master only. Without the reduction a run whose in-band cells all
    // live off-master reports zero and the diagnostic is blind in parallel.
    Foam::reduce(nFlipped, Foam::sumOp<Foam::label>());

    return nFlipped;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::fgmFluid::correctPressurePEP()
{
    volScalarField& rho(rho_);
    volScalarField& p(p_);
    volVectorField& U(U_);
    surfaceScalarField& phi(phi_);

    // --- RANK 1: refresh the real-fluid rho, psi at the CURRENT pressure each
    // pressure corrector (not once per outer). A stale SRK compressibility
    // psi = (drho/dp)_T at the stiff cold-LOX injector cells lets the pressure
    // corrector run against a near-singular diagonal and generate the injector
    // pressure spike -- the pressure-velocity / EOS-stiffness ill-conditioning
    // that neither PEP nor LAD addresses (Ma, Lv & Ihme, J. Comput. Phys. 340
    // (2017) 330). Property-update-per-corrector is the realFluidFoam recipe
    // (Nguyen & Yoo, Comput. Phys. Commun. 312 (2025) 109600). thermo_.correct()
    // inverts the manifold-seeded he to T and refreshes rho/psi/mu at (p,T,Y).
    // (b) he->T-drift-free variant: re-seed he from the manifold at the CURRENT
    // pressure (updateManifold) BEFORE correct(), so the he->T inversion returns
    // T_table exactly (no drift) while rho/psi still refresh at the new p. This
    // fixes the RANK-1 side effect where bare thermo_.correct() inverted a stale
    // he against the new p and drifted T, moving the spike into the chamber.
    //
    // OPT: this manifold-reseed + EOS refresh is ~93% of the pressure-corrector
    // cost and, at nOuter x nCorr, is called 4x/step (the dominant s/step term).
    // thermoPerCorrector=false hoists it to once-per-outer (standard PIMPLE),
    // done in pressureCorrector() before the corrector loop -- cheaper, but the
    // psi/rho used by the 2nd+ correctors is then the outer's start-of-step
    // state (stale), which is exactly the ill-conditioning RANK 1 was added to
    // cure, so it must be validated against the injector spike before use.
    const Switch thermoPerCorrector
    (
        pimple.dict().lookupOrDefault<Switch>("thermoPerCorrector", true)
    );
    if (thermoPerCorrector)
    {
        updateManifold();
        thermo_.correct();
    }

    // Per-corrector transported-density re-sync (modified-PIMPLE, cf.
    // realFluidFoam/Jarczyk-Pfitzner). The transported rho_ is otherwise
    // advanced only by the correctRho(psi*dp) increments; at a flame-zone
    // pressure spike those increments are huge and (with the pMinPa clamp)
    // inconsistent, so rho_ drifts from the EOS state and accumulates error
    // -- observed: rho_ down to -1175 kg/m3 at the ox tangential holes while
    // the EOS density stayed positive. Snapping rho_ to the just-corrected
    // thermo state each corrector removes the drift (rho.oldTime() is
    // untouched, so the ddt history stays consistent).
    //
    // (2026-07-26) This used to be done twice: an `rho_ = thermo.rho();
    // rho_.correctBoundaryConditions();` pair here, immediately followed by
    // the `rho = thermo.rho()` below -- and `rho` is a reference to `rho_`
    // (see the alias at the top of this function), so the first assignment
    // was unconditionally overwritten by the second and had no effect. The
    // dead pair is removed; the re-sync itself is unchanged and now happens
    // exactly once, at the assignment below.
    volScalarField psi(thermo.psi());
    rho = thermo.rho();
    rho.relax();

    // Root-selection hysteresis: correct rho AND psi together, from the
    // SAME chosen branch, before psis (below) is built from them. See
    // correctRootHysteresis() at the top of this file.
    if (pimple.dict().lookupOrDefault<Switch>("rootHysteresis", false))
    {
        const tabulatedRealGasMixture* hook =
            dynamic_cast<const tabulatedRealGasMixture*>(&thermo_);

        if (hook)
        {
            // Route through the tabulated bM/coef1-3/cM (RGcoeffFields_)
            // when armed -- avoids re-deriving them via calculateRealGas()'s
            // O(n^2) mixing (see correctRootHysteresis() comment above).
            tmp<volScalarField> tW;
            const volScalarField* Wptr = nullptr;
            if (tabRealGasCoeffs_)
            {
                tW = thermo_.W();
                Wptr = &tW();
            }

            // Neighbour-consistency cap (2026-07-24f): pre-correction
            // face-averaged rho, so the cap reflects genuine spatial
            // consensus (not hysteresis's own output). Only computed when
            // the cap is actually requested (fvc::average/interpolate cost
            // avoided otherwise).
            const scalar rootHysteresisCapRatio
            (
                pimple.dict()
                   .lookupOrDefault<scalar>("rootHysteresisCapRatio", GREAT)
            );
            tmp<volScalarField> tRhoAvg;
            const volScalarField* rhoAvgPtr = nullptr;
            if (rootHysteresisCapRatio < GREAT)
            {
                tRhoAvg = fvc::average(fvc::interpolate(rho));
                rhoAvgPtr = &tRhoAvg();
            }

            const label nFlipped = correctRootHysteresis
            (
                *hook, Y_, thermo_.p(), thermo_.T(),
                rho_.oldTime().primitiveField(), rho, &psi,
                pimple.dict().lookupOrDefault<scalar>("rootHysteresisTMin", 140),
                pimple.dict().lookupOrDefault<scalar>("rootHysteresisTMax", 200),
                pimple.dict().lookupOrDefault<scalar>("rootHysteresisPMin", 4.0e6),
                pimple.dict().lookupOrDefault<scalar>("rootHysteresisPMax", 5.5e6),
                tabRealGasCoeffs_ ? &RGcoeffFields_ : nullptr,
                Wptr,
                nullptr,
                rhoAvgPtr,
                rootHysteresisCapRatio
            );

            if (nFlipped > 0)
            {
                Info<< "rootHysteresis: " << nFlipped
                    << " cell(s) switched to the history-consistent branch"
                    << " (rho, psi)" << endl;
            }
        }
    }

    fvVectorMatrix& UEqn = tUEqn.ref();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi*p);

    // Face density used to (a) convert the RC transient term to volumetric and
    // (b) RECONSTRUCT the mass flux from the solved volumetric flux (below).
    //
    // 'rhofUpwind' (default off, read each step): use the UPWIND (donor-cell)
    // density instead of the linear interpolate. Rationale, ported from the
    // real-fluid literature: the spurious pressure/velocity at a transcritical
    // contact is attributed to FLUX INCONSISTENCY at the face -- the PEP
    // condition of the conservative schemes is precisely an algebraic
    // consistency requirement between the numerical energy flux and the
    // numerical MASS flux (dF_rhoe = alpha_hat*dF_rho; EPEP-RG, arXiv:2605.03617),
    // and the double-flux family removes the oscillation by evaluating the
    // face state TWICE from the two donor cells rather than once from an
    // averaged state (Ma, Lv & Ihme, JCP 340 (2017) 330).
    //
    // The same inconsistency exists here in discrete form: the convection
    // operators transport with the UPWIND density while phi is rebuilt with a
    // LINEARLY interpolated rhof. Across the LOX/gas contact (rho 25-60x) the
    // two disagree by O(rho_ratio) -- the identical defect already measured and
    // fixed on the Courant number (see setRDeltaT.C 'densityConsistentCourant',
    // maxCo_eff ~0.001 vs the target 0.2).
    const Switch rhofUpwind
    (
        pimple.dict().lookupOrDefault<Switch>("rhofUpwind", false)
    );
    const surfaceScalarField rhof
    (
        "rhof",
        rhofUpwind
      ? upwind<scalar>(mesh, phi).interpolate(rho)
      : fvc::interpolate(rho)
    );

    const volScalarField rAU("rAU", 1.0/UEqn.A());
    const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    tmp<volScalarField> rAtU
    (
        pimple.consistent()
      ? volScalarField::New("rAtU", 1.0/(1.0/rAU - UEqn.H1()))
      : tmp<volScalarField>(nullptr)
    );

    tmp<surfaceScalarField> rhorAtUf
    (
        pimple.consistent()
      ? surfaceScalarField::New("rhoRAtUf", fvc::interpolate(rho*rAtU()))
      : tmp<surfaceScalarField>(nullptr)
    );

    const volScalarField& rAAtU = pimple.consistent() ? rAtU() : rAU;
    const surfaceScalarField& rhorAAtUf =
        pimple.consistent() ? rhorAtUf() : rhorAUf;

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    if (pimple.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    // --- PEP (pressure-evolution) pressure equation -------------------------
    // Replace the base continuity pEqn's fvc::ddt(rho) -- whose Eulerian
    // thermal density change radiates the spurious contact pressure spike --
    // with the pressure-evolution form  psis*dp/dt + div(u) = 0, where psis =
    // 1/(rho c^2) is the isentropic compressibility and div(u) is the
    // VOLUMETRIC velocity divergence (NOT the mass-flux divergence). For a
    // uniform-p, uniform-u contact div(u) = 0 and dp/dt = 0 by construction, so
    // the density jump advects without radiating pressure (Terashima & Koshi,
    // J. Comput. Phys. 231 (2012) 6907; Kai/Kurose PEQSI, Phys. Fluids 36
    // (2024) 116104). ATTEMPT 1: ideal-gas isentropic compressibility
    // 1/(gamma p) -- validates the STRUCTURE (spike removal); refined to the
    // real SRK sound speed next. constrainPressure / MRF-mass / consistent
    // branches dropped for this first cyclic-benchmark attempt.
    // C1a: HARMONIC-mean face mobility (Rhie-Chow across the stiff density
    // jump). rAU = 1/A() ~ dt/rho jumps ~60x cell-to-cell at the recess-tip
    // LOx/gas interface; the arithmetic fvc::interpolate(rAU) over-weights the
    // low-density (large rAU) side and over-stiffens the pEqn face coefficient
    // by ~O(30), so a small dp drives a huge face flux -> the spurious pressure
    // spike with |U| pinned at the limiter (co-located max(p)/|U|-cap). The
    // harmonic (series-resistance-correct) mean is the density-jump-robust
    // choice for a face mobility (Ferziger & Peric; Rhie & Chow, AIAA J. 21
    // (1983) 1525). 1/rAU = A() > 0 (momentum diagonal) so no floor needed.
    const surfaceScalarField rAUf("rAUf", 1.0/fvc::interpolate(1.0/rAU));
    // ATTEMPT 2: real SRK isothermal compressibility psis = (drho/dp)_T/rho =
    // kappa_T (thermo.psi() now returns the real (drho/dp)_T -- see SRKGasI.H).
    // This supplies the true dense-fluid stiffness that the ideal-gas 1/(gamma
    // p) lacked (~100x too soft for liquid LOX).
    // C3: pressure-based double-flux analogue -- bound the cell-to-cell
    // CONTRAST of the pEqn diagonal compressibility psis = kappa_T across the
    // stiff-LOX / soft-gas interface. psis jumps enormously there (dense LOX
    // stiff, warm gas soft); a large diagonal-ratio next to the (now harmonic)
    // off-diagonal still lets a residual dp blow up. Floor each cell's psis at
    // (face-averaged psis)/psisCapRatio so the interface diagonal ratio is
    // capped -- the honest pressure-based cousin of freezing a local effective
    // gamma* (Ma, Lv & Ihme, JCP 340 (2017) 330). psisCapRatio default GREAT
    // (off); read each step. Conservative at convergence (psis*ddt(p) -> 0).
    const scalar psisCapRatio
    (
        pimple.dict().lookupOrDefault<scalar>("psisCapRatio", GREAT)
    );
    // psisIsentropic: use the ISENTROPIC compressibility kappa_s = kappa_T/gamma
    // = 1/(rho c^2) instead of the isothermal kappa_T = thermo.psi()/rho. For a
    // pressure-EVOLUTION (PEP) equation the acoustically-correct diagonal is the
    // isentropic one; the isothermal psi() over-stiffens the dense cold LOX
    // (Terashima-Koshi; design-agent refinement). Default off.
    const Switch psisIsentropic
    (
        pimple.dict().lookupOrDefault<Switch>("psisIsentropic", false)
    );
    // pepFull (2026-07-18): COMPLETE Terashima-Koshi pressure-evolution
    // equation. The half-PEP (isentropic diagonal alone) measurably WORSENED
    // the interface trap metrics (A/B: dome capped cells 23.6k -> 76.8k,
    // downstream 3.9k -> 28.1k over one equal segment) because the full
    // formulation balances the acoustically-correct diagonal with TWO terms
    // the lite version lacks:
    //   psis*(dp/dt + u.grad(p)) + div(u) = (gamma-1)*psis*div(Deff grad h)
    // (Terashima & Koshi JCP 231 (2012) 6907, volumetric form; the RHS is the
    // thermal/species-diffusion expansion source written with the SAME
    // effective enthalpy diffusivity as the h equation for discrete
    // consistency; reaction expansion enters through the per-corrector
    // rho(T) refresh). pepFull implies the isentropic diagonal. Default off.
    const Switch pepFull
    (
        pimple.dict().lookupOrDefault<Switch>("pepFull", false)
    );
    // v2 term-isolation switches for the 1-D benchmark: pepAdvect enables the
    // pressure-advection term alone, pepRHS the diffusive-expansion RHS alone.
    // pepFull implies both. Each read every step.
    const Switch pepAdvect
    (
        pepFull
     || pimple.dict().lookupOrDefault<Switch>("pepAdvect", false)
    );
    const Switch pepRHS
    (
        pepFull
     || pimple.dict().lookupOrDefault<Switch>("pepRHS", false)
    );
    // psisFreezeOuter (2026-07-24): double-flux/RFQC-literal "freeze the
    // acoustic coefficient over the step" -- psis is computed ONCE per OUTER
    // iteration (in pressureCorrector(), before the corrector loop, stored in
    // the registry as "psisFrozen") instead of being re-evaluated from the
    // current (highly nonlinear near the transcritical interface) EOS state
    // every corrector. This is a more literal implementation of the double-
    // flux idea than psisAdvect (a separate material-transport PDE for psis,
    // tested 2026-07-23 and found to DESTABILISE further) -- no new PDE, just
    // hoisting psis out of the corrector loop the same way thermoPerCorrector
    // =false already hoists the manifold/EOS refresh. Default off.
    const Switch psisFreezeOuter
    (
        pimple.dict().lookupOrDefault<Switch>("psisFreezeOuter", false)
    );
    // psisTabulated (2026-07-24): use the per-cell FGM:psisTab field (offline
    // SRK+JANAF re-derivation of psis, gamma-clipped/smoothed at manifold-
    // build time -- build_psis_table_v2.py) instead of the pointwise runtime
    // EOS evaluation. Highest priority in the chain below: the whole point is
    // to bypass BOTH the raw pointwise nonlinearity AND the ad-hoc
    // psisCapRatio patch that reacts to it. No-op (falls through) if the
    // table wasn't built with a psisTab field (tabPsis_ false). Default off.
    const Switch psisTabulated
    (
        pimple.dict().lookupOrDefault<Switch>("psisTabulated", false)
    );
    volScalarField psis
    (
        "psis",
        (psisTabulated && tabPsis_)
      ? volScalarField(psisTabField_())
      : (psisFreezeOuter && mesh.foundObject<volScalarField>("psisFrozen"))
      ? volScalarField(mesh.lookupObject<volScalarField>("psisFrozen"))
      : (psisIsentropic || pepFull)
      ? volScalarField(psi/(rho*thermo.gamma()))
      : volScalarField(psi/rho)
    );
    // psisSmooth (2026-07-24): replace psisCapRatio's crude single-neighbour-
    // average floor with the SAME established iterative Laplacian-type
    // smoothing (fvc::smooth) this solver already uses for rDeltaT
    // (rDeltaTSmoothingCoeff, setRDeltaT.C) -- a more principled spatial
    // regularisation than an ad-hoc cap, reusing a well-tested OpenFOAM
    // operator instead of inventing one. Mutually exclusive with psisCapRatio
    // (takes precedence if on). Default off.
    const Switch psisSmooth
    (
        pimple.dict().lookupOrDefault<Switch>("psisSmooth", false)
    );
    const scalar psisSmoothCoeff
    (
        pimple.dict().lookupOrDefault<scalar>("psisSmoothCoeff", 0.1)
    );
    // Skip both regularisers when psisTabulated is active: the tabulated
    // field is already smoothed offline (gamma-clipped at build time), so a
    // runtime patch on top would just re-introduce the mesh-local noise the
    // table was built to avoid.
    if (psisTabulated && tabPsis_)
    {
        // no-op
    }
    else if (psisSmooth)
    {
        fvc::smooth(psis, psisSmoothCoeff);
    }
    else if (psisCapRatio < GREAT)
    {
        const volScalarField psisSm(fvc::average(fvc::interpolate(psis)));
        // Diagnostic (2026-07-23): count how many cells the neighbour-average
        // floor actually overrides each corrector, to test whether pepFull's
        // long-time drift correlates with a growing psisCapRatio-affected
        // (near-interface, nonlinear-EOS) cell population.
        const scalarField& psisIf = psis.primitiveField();
        const scalarField psisFloor(psisSm.primitiveField()/psisCapRatio);
        label nCapped = 0;
        forAll(psisIf, celli)
        {
            if (psisIf[celli] < psisFloor[celli])
            {
                nCapped++;
            }
        }
        Info<< "psisCapRatio: " << nCapped << " of " << psisIf.size()
            << " cells capped (" << (100.0*nCapped/psisIf.size()) << "%)"
            << endl;
        psis = max(psis, psisSm/psisCapRatio);
        psis.correctBoundaryConditions();
    }

    // --- psisAdvect: ADVECTED acoustic coefficient (RFQC port) ---------------
    // The quasi-conservative real-fluid family does NOT re-evaluate the
    // thermodynamic coefficient pointwise from the cubic EOS every step at a
    // contact: the double-flux model FREEZES gamma* = rho c^2/p and e*_0 over
    // the step (Ma, Lv & Ihme, JCP 340 (2017) 330), and its generalisation RFQC
    // ADVECTS the Grueneisen-related coefficient xi = h/c^2 and the remainder
    // E0 with the material, d_t xi + u.grad(xi) = 0 (Bai et al., JCP 564 (2026)
    // 115156), recovering p = (rho e - E0)/xi. Pointwise re-evaluation is what
    // injects the spurious pressure, because the coefficient jumps
    // discontinuously as the interface sweeps a cell.
    //
    // In a PRESSURE-BASED solver the acoustic coefficient is not an energy-flux
    // weight but the pEqn DIAGONAL psis = kappa (d rho/dp)/rho. No published
    // pressure-based double-flux implementation exists (all are density-based;
    // the open question in the literature is exactly WHAT the freeze should
    // modify in a segregated algorithm), so the port here is: transport psis as
    // its own scalar with a slow relaxation back to the EOS value, and use the
    // transported field as the pEqn diagonal.
    //
    //     d(psisStar)/dt + div(phiv, psisStar) = (psisEOS - psisStar)/tau
    //
    // tau is a PHYSICAL relaxation time, so the formulation is well defined
    // under LTS (where "frozen over the time step" is ambiguous -- the local
    // dt differs per cell). tau -> 0 recovers the present pointwise EOS value;
    // tau -> large is a pure freeze. The relaxation is the continuous analogue
    // of the RFQC re-projection (local O(dt^2), global O(dt) there).
    //
    // psisAdvect (default off) + psisTau [s] (default 1e-5), read each step.
    const Switch psisAdvect
    (
        pimple.dict().lookupOrDefault<Switch>("psisAdvect", false)
    );
    if (psisAdvect)
    {
        const dimensionedScalar psisTau
        (
            "psisTau",
            dimTime,
            pimple.dict().lookupOrDefault<scalar>("psisTau", 1e-5)
        );

        if (!mesh.foundObject<volScalarField>("psisStar"))
        {
            // zeroGradient BCs so the implicit psisEqn can form patch
            // coefficients (psis itself carries 'calculated' BCs, which have no
            // matrix contribution and abort in valueInternalCoeffs).
            mesh.objectRegistry::store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "psisStar",
                        runTime.name(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    psis,
                    wordList
                    (
                        mesh.boundary().size(),
                        zeroGradientFvPatchField<scalar>::typeName
                    )
                )
            );
        }

        volScalarField& psisStar =
            mesh.lookupObjectRef<volScalarField>("psisStar");

        const surfaceScalarField phivPsis("phivPsis", phi/rhof);

        fvScalarMatrix psisEqn
        (
            fvm::ddt(psisStar)
          + fvm::div(phivPsis, psisStar)
          + fvm::Sp(1.0/psisTau, psisStar)
         ==
            psis/psisTau
        );
        psisEqn.solve();

        // Positivity guard: the pEqn diagonal must stay strictly positive.
        psisStar.max(SMALL*dimensionedScalar(psis.dimensions(), 1));
        psisStar.correctBoundaryConditions();

        psis = psisStar;
    }

    // Volumetric predicted face flux with the RHO-CONSISTENT transient
    // Rhie-Chow correction (C2 proper). The base compressible form is the MASS
    // flux  rhof*fvc::flux(HbyA) + rhorAUf*fvc::ddtCorr(rho,U,phi,rhoUf)
    // (isothermalFluid::correctPressure); the volumetric PEP is its /rhof. The
    // rho-aware ddtCorr(rho,U,phi,rhoUf) (rhoUf is null on a static mesh, where
    // the overload falls back to the rho-weighted form) keeps the transient RC
    // CONSISTENT with the momentum flux fvm::div(phi,U) across the ~60x rho jump
    // -- unlike the earlier ddtCorr(U, phi/rhof), whose arithmetic rhof in the
    // denominator injected a spurious face velocity (the recess-tip/injector
    // spike). This replaces the rcDdtScale 0 workaround, which suppressed the
    // spurious RC flux but left the injector fine cells to CHECKERBOARD (a new
    // spike re-formed there). rcDdtScale (default 1) still gates the term for
    // A/B testing. Read each step (runTimeModifiable).
    const scalar rcDdtScale
    (
        pimple.dict().lookupOrDefault<scalar>("rcDdtScale", scalar(1))
    );
    surfaceScalarField phiHbyAv
    (
        "phiHbyAv",
        fvc::flux(HbyA)
      + rcDdtScale*rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf)/rhof  // rho-consistent RC
    );
    MRF.makeRelative(phiHbyAv);

    // Update the pressure BCs for flux consistency (3D: waveTransmissive
    // outlet, fixedFluxPressure walls). Volumetric-flux form: pass the
    // volumetric predicted flux and the rAUf (velocity-level) coefficient,
    // matching the pEqn's laplacian(rAUf, p). The rho overload is for MASS
    // flux; on a flow-carrying fixedFluxPressure patch it would subtract
    // rho_b*(Sf&U_b) from the volumetric phiHbyAv and skew snGrad(p) by
    // ~rho_b (LOX ~1000x). (backport RGP-13-GPU 2005bc9)
    constrainPressure(p, U, phiHbyAv, rAUf, MRF);

    // --- RANK 4: Artificial Mass Diffusivity (AMD) on DENSITY -----------------
    // Kawai, Terashima & Negishi, J. Comput. Phys. 300 (2015) 116: at a large
    // density ratio, diffusing the SCALARS/temperature (as the manifold LAD does
    // on Z,C,h) drives spurious p/u oscillations THROUGH the nonlinear EOS,
    // whereas diffusing MASS/DENSITY is the consistent choice. Smooth the steep
    // transcritical density interface at source by adding a density-gradient-
    // sensed diffusion of rho to the (volumetric) continuity, WITHOUT touching
    // the manifold scalars. Kinematic AMD coefficient [m^2/s], sized to the
    // local cell and active only where rho varies:
    //   Dr = LADrhoCoeff * V^(1/3) * |U| * s,  s = min(|grad rho| V^(1/3)/rho, 1)
    // where s is a DIMENSIONLESS density-gradient sensor CAPPED at 1: the raw
    // |grad rho| blows up at the spurious spike itself, so an uncapped Dr
    // violates the diffusive CFL (Dr ~ 164 m^2/s -> dt collapse). The cap ties
    // Dr_max = LADrhoCoeff*V^(1/3)*|U| to the convective scale, keeping the
    // diffusive step >= the convective one (Olson & Lele, JCP 246 (2013) 207).
    // The continuity then gains -(1/rho) div(Dr grad rho) [1/s], matching the
    // volumetric pressure-evolution form. LADrhoCoeff read each step (default 0).
    const scalar LADrhoCoeff
    (
        pimple.dict().lookupOrDefault<scalar>("LADrhoCoeff", scalar(0))
    );
    // ladOddEven: gate the AMD by a Jameson-type odd-even (checkerboard)
    // detector on PRESSURE, so the mass diffusion fires ONLY on the spurious
    // interface cells and not on the physical LOx/gas contact or the physical
    // high-speed jet. The plain |grad rho| sensor cannot distinguish the two:
    // BOTH the physical contact and the spurious spike carry a steep density
    // gradient, so boosting LADrhoCoeff globally smears the physical jet
    // (measured: the physical cold-flow ceiling is ~112 m/s Bernoulli, yet the
    // spurious tail reaches the limitU clamp of 800). The distinguishing
    // signature is PRESSURE: across a real contact the pressure is smooth
    // (mechanical equilibrium) while the spurious velocity is driven by a
    // cell-to-cell pressure CHECKERBOARD. The detector
    //   theta = |lap(p)|*dx^2 / (|lap(p)|*dx^2 + |grad(p)|*dx + eps)  in [0,1)
    // is the ratio of the undivided 2nd difference to the 1st difference: ~1 at
    // an odd-even oscillation (2nd diff >> 1st diff), ~0 for a smooth or
    // monotone pressure field, and ~0 in quiescent cells (eps dominates). With
    // the gate on, LADrhoCoeff can be pushed hard (2-4) to dissipate the
    // spurious spike WITHOUT touching the 95% physical bulk -- unlike the
    // earlier uniform LAD boost (Seg A), which smeared everything and worsened
    // dt. ladOddEven default off; read each step. epsPa is the quiescent-cell
    // pressure floor (default 1e3 Pa).
    const Switch ladOddEven
    (
        pimple.dict().lookupOrDefault<Switch>("ladOddEven", false)
    );
    const scalar ladEpsPa
    (
        pimple.dict().lookupOrDefault<scalar>("ladEpsPa", scalar(1e3))
    );
    volScalarField Dr
    (
        IOobject
        (
            "Dr_amd",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimArea/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    );
    if (LADrhoCoeff > 0)
    {
        const scalarField V13(pow(scalarField(mesh.V()), 1.0/3.0));
        // Dimensionless density-gradient sensor in [0,1]: the relative rho
        // change across a cell, CAPPED at 1 so the spurious spike itself (huge
        // |grad rho|) cannot drive the coefficient past the diffusive-CFL limit.
        scalarField sensor
        (
            min
            (
                mag(fvc::grad(rho))().primitiveField()*V13/rho.primitiveField(),
                scalar(1)
            )
        );
        if (ladOddEven)
        {
            const scalarField V23(V13*V13);
            const scalarField d2p
            (
                mag(fvc::laplacian(p))().primitiveField()*V23     // 2nd diff [Pa]
            );
            const scalarField d1p
            (
                mag(fvc::grad(p))().primitiveField()*V13          // 1st diff [Pa]
            );
            const scalarField theta(d2p/(d2p + d1p + ladEpsPa));
            sensor *= theta;
            Info<< "LAD-rho oddEven: theta max = " << gMax(theta) << endl;
        }
        Dr.primitiveFieldRef() =
            LADrhoCoeff*V13*mag(U)().primitiveField()*sensor;
        Dr.correctBoundaryConditions();
        Info<< "LAD-rho: Dr_amd max = " << gMax(Dr.primitiveField())
            << " m^2/s" << endl;
    }

    fvScalarMatrix pDDtEqn
    (
        psis*fvm::ddt(p)
      + fvc::div(phiHbyAv)
      - fvc::laplacian(Dr, rho)/rho
    );

    // --- pepFull: pressure ADVECTION + energy-diffusion RHS ----------------
    // (1) psis*(u.grad p) in exact flux form: div(F p) - p div(F) = F.grad(p)
    //     with F = psisf * (volumetric flux). Uses the CURRENT mass flux
    //     phi/rhof (last corrector) as the advecting velocity.
    // (2) RHS (gamma-1)*psis*laplacian(Dh, h): the diffusive expansion source
    //     of the full pressure-evolution equation, discretised with the SAME
    //     Deff("h") as the h transport for consistency.
    if (pepAdvect)
    {
        const surfaceScalarField phivAdv("phivAdv", phi/rhof);
        // psis face value: HARMONIC mean (v2) -- the arithmetic interpolate
        // over-weights the soft-gas side across the ~250x psis jump, leaving
        // a discrete residual in the div(F p) - p div(F) identity that acts
        // as a spurious source at the interface (v1 suspect #3).
        const surfaceScalarField Fpsis
        (
            "Fpsis", (1.0/fvc::interpolate(1.0/psis))*phivAdv
        );
        // v2.2: EXPLICIT advection source. Implicit forms (fvm::div + Sp or
        // SuSp pairs, v1/v2.1) both blew up on the 1-D benchmark -- the
        // CONCEPTUAL flaw is that the PIMPLE pressure equation is a
        // CORRECTION-form solve, and advecting the FULL p implicitly inside
        // it mixes correction and total quantities (Terashima-Koshi advance p
        // directly, outside a correction structure). Explicit treatment keeps
        // the matrix symmetric (PCG-compatible) and defers the advection to
        // the outer-corrector lag, which nOuter >= 2 absorbs.
        pDDtEqn +=
            fvc::div(Fpsis*fvc::interpolate(p)) - p*fvc::div(Fpsis);
    }
    if (pepRHS)
    {
        if (fgmTable_.useEnthalpy())
        {
            const volScalarField& h = hPtr_();
            const volScalarField DhP("DhP", Deff("h"));
            pDDtEqn -=
                (thermo.gamma() - scalar(1))*psis
               *fvc::laplacian(DhP, h);
        }
    }

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rAUf, p));

        pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        fvConstraints().constrain(pEqn);

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            // Reconstruct the mass flux from the corrected volumetric flux
            phi = rhof*(phiHbyAv + pEqn.flux());
        }
    }

    // --- POSITIVITY GUARD: floor the SOLVED pressure to a physical minimum ---
    // ROOT ISSUE (2026-07-02): the swirl low-pressure core + stiff cold-LOX SRK
    // drive the pEqn to undershoot the SOLVED pressure NEGATIVE (min p ~ -1e7
    // observed even after the EOS-internal density floor). Absolute p<0 is
    // thermodynamically impossible; via SRK rho = p/(Z R T) it makes rho<0 and,
    // through -grad(p), pumps |U| overshoots that collapse dt. The EOS-internal
    // floor only protects rho -- the solved p field itself must also be bounded.
    // Floor p here (right after the solve, before correctRho/relax/U) so every
    // downstream step sees a positive pressure. pMinPa read each step (default 0
    // = off). Positivity guard confined to sub-physical cells, paired with the
    // EOS floor; a physically-exact version would use p_sat(T) instead.
    const scalar pMinPa
    (
        pimple.dict().lookupOrDefault<scalar>("pMinPa", scalar(0))
    );
    if (pMinPa > 0)
    {
        p = max(p, dimensionedScalar("pMinPa", p.dimensions(), pMinPa));
        p.correctBoundaryConditions();
    }

    // Upper positivity-guard twin of pMinPa: cap the solved pressure. Used as
    // a TEMPORARY surgical bound while the ignition-phase flame-zone spike
    // dissolves (fvConstraints limitPressure proved inert here even after
    // constrain(p) -- root cause not chased; this knob is on the proven
    // pMinPa path). pMaxPa read each step (default 0 = off).
    const scalar pMaxPa
    (
        pimple.dict().lookupOrDefault<scalar>("pMaxPa", scalar(0))
    );
    // pMaxXmin (2026-07-18): x-gate the upper pressure cap. The LOX dome
    // tangential-jet impingement stagnation is PHYSICAL up to
    // rho U^2/2 ~ 1100*(240)^2/2 ~ 3.2e7 Pa; a global pMaxPa sized for the
    // annulus-interface trap (~1e7) CLIPS that physics (observed: 77% of the
    // capped cells sitting in the dome at the injection radius). Applying the
    // cap only for cells with x > pMaxXmin (default -GREAT = everywhere)
    // leaves the dome unclipped while still guarding the downstream
    // liquid-film/interface region. Read each step (runTimeModifiable).
    const scalar pMaxXmin
    (
        pimple.dict().lookupOrDefault<scalar>("pMaxXmin", -great)
    );
    // pMaxPaDome (2026-07-18): two-tier cap. Releasing the dome entirely
    // (pMaxXmin gating alone) let the SAME trapped-pressure instability
    // re-inflate THERE (observed 54,127 bar at the LOX dome, 170x the
    // physical jet-stagnation ~317 bar) -- the dome liquid shares the
    // stiff-liquid pEqn trap mechanism with the annulus film. The dome
    // therefore needs its own, LOOSER cap: high enough for physical
    // impingement stagnation (rho U^2/2 ~ 3.2e7), low enough to break the
    // runaway loop. Default 0 = dome uncapped. Read each step.
    const scalar pMaxPaDome
    (
        pimple.dict().lookupOrDefault<scalar>("pMaxPaDome", scalar(0))
    );
    if (pMaxPa > 0)
    {
        if (pMaxXmin > -great/2)
        {
            const vectorField& cc = mesh.C().primitiveField();
            scalarField& pf = p.primitiveFieldRef();
            forAll(pf, celli)
            {
                if (cc[celli].x() > pMaxXmin)
                {
                    if (pf[celli] > pMaxPa) pf[celli] = pMaxPa;
                }
                else if (pMaxPaDome > 0 && pf[celli] > pMaxPaDome)
                {
                    pf[celli] = pMaxPaDome;
                }
            }
        }
        else
        {
            p = min(p, dimensionedScalar("pMaxPa", p.dimensions(), pMaxPa));
        }
        p.correctBoundaryConditions();
    }

    // Thermodynamic density update: the stock SIMPLErho increment
    // correctRho(psi*dp) is DISABLED in the PEP path. At a flame-zone
    // pressure spike the increment is huge (psi_gas ~1e-5 x dp ~ -4e8 ->
    // drho ~ -4000) and drove the THERMO density to -4869 kg/m3 in written
    // states -- and the per-corrector updateManifold()+thermo_.correct()+
    // rho_ re-sync at the top of this function already provides the full
    // EOS-consistent density update (one-corrector lag, nOuter >= 3).

    // Continuity diagnostics (base isothermalFluid::continuityErrors wraps
    // this grandparent call)
    fluidSolver::continuityErrors(rho, thermo.rho(), phi);

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Apply the field-level fvConstraints (limitPressure min/max). The stock
    // correctPressure applies them after the solve; the PEP override had
    // omitted this call, leaving constant/fvConstraints limitP INERT (both
    // its min and the surgical max) -- only the pMinPa floor was active.
    fvConstraints().constrain(p);

    // faceGradP: balanced-force face-consistent pressure gradient for the CELL
    // velocity update. The Green-Gauss fvc::grad(p) rings cell-to-cell at the
    // stiff liquid/gas interface (annulus rho 870<->30) while the FACE snGrad
    // driving phi stays smooth -- the observed cell-U checkerboard with a
    // coherent phi. Reconstructing the gradient FROM the face snGrad makes the
    // cell velocity feel exactly the face pressure differences that drive the
    // flux, removing the decoupling mechanism at the same operation count
    // (Francois et al., balanced-force, JCP 213 (2006) 141; interFoam-family
    // standard). Read each step (runTimeModifiable); default off.
    const Switch faceGradP
    (
        pimple.dict().lookupOrDefault<Switch>("faceGradP", false)
    );
    if (faceGradP)
    {
        // safeReconstruct: see safeReconstruct.H (AMR hanging-node
        // singular-tensor SIGFPE guard, 2026-07-23). Falls back to the same
        // fvc::grad(p) the non-faceGradP branch below uses.
        U = HbyA - rAAtU*safeReconstruct(fvc::snGrad(p)*mesh.magSf(), fvc::grad(p));
    }
    else
    {
        U = HbyA - rAAtU*fvc::grad(p);
    }
    U.correctBoundaryConditions();
    fvConstraints().constrain(U);
    K = 0.5*magSqr(U);

    // SIMPLErho: density from the equation of state
    if (pimple.simpleRho())
    {
        rho = thermo.rho();
        rho.relax();

        // Root-selection hysteresis (2026-07-24): near the O2 near-critical
        // band the hard minimum-fugacity switch in SRKGasI.H::Z() can flip
        // liquid/vapour branch between adjacent cells (or between
        // iterations) for p perturbations as small as 0.05-3%, jumping rho
        // by 2-5x -- see PEP and LAD Spike Suppression wiki secs 6-9. The
        // rootBlendTol logistic blend (secs 10-14) fixed that but corrupted
        // the cold-start pressure field (2-3x nominal overshoot, root cause
        // not isolated after two redesign attempts -- sec 14-16). This is a
        // different, non-blending fix: for cells with a genuine 3-real-root
        // ambiguity, pick whichever ACTUAL root (never an interpolated
        // non-root value) gives a density closer to this cell's PREVIOUS-
        // timestep density -- "sticky" branch selection driven by the
        // cell's own history rather than an instantaneous (and, near the
        // crossing, essentially noise-sensitive) fugacity comparison.
        // MUST run here too, after this LAST unconditional rho =
        // thermo.rho() in this function -- this is the value that reaches
        // fvc::correctRhoUf and next-timestep's oldTime() (an earlier
        // attempt placed the ONLY correction right after the FIRST such
        // assignment near the top of correctPressurePEP() and was silently
        // overwritten by this one -- diagnosed 2026-07-24 by a matched
        // hysteresis-on/off A/B that came back bit-identical). psi is NOT
        // re-corrected here: it was already made rho-consistent, earlier in
        // this function, at the point psis is actually built from it (see
        // correctRootHysteresis() call above and at the top of this file);
        // nothing downstream of this point reads psi again. Off by default;
        // enable with 'rootHysteresis true;' in PIMPLE.
        if (pimple.dict().lookupOrDefault<Switch>("rootHysteresis", false))
        {
            const tabulatedRealGasMixture* hook =
                dynamic_cast<const tabulatedRealGasMixture*>(&thermo_);

            if (hook)
            {
                tmp<volScalarField> tW;
                const volScalarField* Wptr = nullptr;
                if (tabRealGasCoeffs_)
                {
                    tW = thermo_.W();
                    Wptr = &tW();
                }

                // Diagnostic (2026-07-24e): 'rootHysteresisDiag true;' writes
                // a per-cell "FGM:hystCorrected" flag (1 = this cell's final
                // rho was hysteresis-overridden this corrector) so it can be
                // checked post-hoc for spatial correlation with neighbour
                // density jumps -- see PEP and LAD Spike Suppression wiki
                // sec 20.
                volScalarField* flagPtr = nullptr;
                if
                (
                    pimple.dict()
                       .lookupOrDefault<Switch>("rootHysteresisDiag", false)
                )
                {
                    if (mesh.foundObject<volScalarField>("FGM:hystCorrected"))
                    {
                        flagPtr = &mesh.lookupObjectRef<volScalarField>
                        (
                            "FGM:hystCorrected"
                        );
                    }
                    else
                    {
                        flagPtr = new volScalarField
                        (
                            IOobject
                            (
                                "FGM:hystCorrected",
                                mesh.time().name(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            mesh,
                            dimensionedScalar(dimless, 0)
                        );
                        mesh.objectRegistry::store(flagPtr);
                    }
                }

                const scalar rootHysteresisCapRatio
                (
                    pimple.dict().lookupOrDefault<scalar>
                    (
                        "rootHysteresisCapRatio", GREAT
                    )
                );
                tmp<volScalarField> tRhoAvg;
                const volScalarField* rhoAvgPtr = nullptr;
                if (rootHysteresisCapRatio < GREAT)
                {
                    tRhoAvg = fvc::average(fvc::interpolate(rho_));
                    rhoAvgPtr = &tRhoAvg();
                }

                const label nFlipped = correctRootHysteresis
                (
                    *hook, Y_, thermo_.p(), thermo_.T(),
                    rho_.oldTime().primitiveField(), rho_, nullptr,
                    pimple.dict().lookupOrDefault<scalar>("rootHysteresisTMin", 140),
                    pimple.dict().lookupOrDefault<scalar>("rootHysteresisTMax", 200),
                    pimple.dict().lookupOrDefault<scalar>("rootHysteresisPMin", 4.0e6),
                    pimple.dict().lookupOrDefault<scalar>("rootHysteresisPMax", 5.5e6),
                    tabRealGasCoeffs_ ? &RGcoeffFields_ : nullptr,
                    Wptr,
                    flagPtr,
                    rhoAvgPtr,
                    rootHysteresisCapRatio
                );

                if (nFlipped > 0)
                {
                    Info<< "rootHysteresis: " << nFlipped
                        << " cell(s) switched to the history-consistent"
                        << " branch (final rho)" << endl;
                }
            }
        }
    }

    // Correct rhoUf if the mesh is moving
    fvc::correctRhoUf(rhoUf, rho, U, phi, MRF);

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);
    }
}


void Foam::solvers::fgmFluid::pressureCorrector()
{
    if (buoyancy.valid())
    {
        FatalErrorInFunction
            << "fgmFluid PEP pressure corrector does not support buoyant "
            << "(p_rgh) cases." << exit(FatalError);
    }

    // OPT: when the per-corrector FGM/EOS refresh is disabled, do it once here
    // (standard PIMPLE: thermo updated once per outer, before the correctors).
    const Switch thermoPerCorrector
    (
        pimple.dict().lookupOrDefault<Switch>("thermoPerCorrector", true)
    );
    if (!thermoPerCorrector)
    {
        updateManifold();
        thermo_.correct();
    }

    // psisFreezeOuter: compute psis ONCE here (double-flux/RFQC-literal
    // freeze, see the psis block in correctPressurePEP() for the full
    // rationale) using the JUST-refreshed (or previous-outer-final, if
    // thermoPerCorrector=true) rho/thermo state, store it in the registry so
    // every corrector this outer reuses the SAME value instead of
    // recomputing pointwise from the current (possibly still-settling)
    // thermodynamic state.
    const Switch psisFreezeOuter
    (
        pimple.dict().lookupOrDefault<Switch>("psisFreezeOuter", false)
    );
    if (psisFreezeOuter)
    {
        const Switch psisIsentropic
        (
            pimple.dict().lookupOrDefault<Switch>("psisIsentropic", false)
        );
        const Switch pepFullSw
        (
            pimple.dict().lookupOrDefault<Switch>("pepFull", false)
        );

        // Root-selection hysteresis: keep this frozen-psis snapshot
        // rho/psi-consistent too (see correctRootHysteresis() at the top of
        // this file / the identical correction in correctPressurePEP()).
        volScalarField psiNow(thermo.psi());
        if (pimple.dict().lookupOrDefault<Switch>("rootHysteresis", false))
        {
            const tabulatedRealGasMixture* hook =
                dynamic_cast<const tabulatedRealGasMixture*>(&thermo_);

            if (hook)
            {
                tmp<volScalarField> tW;
                const volScalarField* Wptr = nullptr;
                if (tabRealGasCoeffs_)
                {
                    tW = thermo_.W();
                    Wptr = &tW();
                }

                const scalar rootHysteresisCapRatio
                (
                    pimple.dict().lookupOrDefault<scalar>
                    (
                        "rootHysteresisCapRatio", GREAT
                    )
                );
                tmp<volScalarField> tRhoAvg;
                const volScalarField* rhoAvgPtr = nullptr;
                if (rootHysteresisCapRatio < GREAT)
                {
                    tRhoAvg = fvc::average(fvc::interpolate(rho_));
                    rhoAvgPtr = &tRhoAvg();
                }

                correctRootHysteresis
                (
                    *hook, Y_, thermo_.p(), thermo_.T(),
                    rho_.oldTime().primitiveField(), rho_, &psiNow,
                    pimple.dict().lookupOrDefault<scalar>("rootHysteresisTMin", 140),
                    pimple.dict().lookupOrDefault<scalar>("rootHysteresisTMax", 200),
                    pimple.dict().lookupOrDefault<scalar>("rootHysteresisPMin", 4.0e6),
                    pimple.dict().lookupOrDefault<scalar>("rootHysteresisPMax", 5.5e6),
                    tabRealGasCoeffs_ ? &RGcoeffFields_ : nullptr,
                    Wptr,
                    nullptr,
                    rhoAvgPtr,
                    rootHysteresisCapRatio
                );
            }
        }

        const volScalarField psisNow
        (
            (psisIsentropic || pepFullSw)
          ? psiNow/(rho_*thermo.gamma())
          : psiNow/rho_
        );
        if (!mesh.foundObject<volScalarField>("psisFrozen"))
        {
            mesh.objectRegistry::store
            (
                new volScalarField("psisFrozen", psisNow)
            );
        }
        else
        {
            mesh.lookupObjectRef<volScalarField>("psisFrozen") = psisNow;
        }
    }

    while (pimple.correct())
    {
        correctPressurePEP();
    }

    tUEqn.clear();
}


// ************************************************************************* //
