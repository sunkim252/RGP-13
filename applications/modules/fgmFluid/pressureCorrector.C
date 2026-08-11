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
#include "fvcCurl.H"
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

    // rhoResyncCoeff (2026-08-07): how much of the CONTINUITY density is
    // discarded in favour of the EOS density each corrector.
    //
    //     rho <- (1 - a)*rho_continuity + a*thermo.rho()
    //
    // a = 1 is the historical behaviour: the pressure equation's own density
    // (the one it just solved the continuity equation FOR, incremented by
    // psi*dp) is overwritten outright. That is why `mcpCorrectRho` measured
    // EXACTLY zero effect -- the increment is wiped on the next corrector's
    // first line. It is also why the base solver framework passes BOTH
    // densities to fluidSolver::continuityErrors(rho, thermo.rho(), phi):
    // they are meant to be allowed to differ, and the gap is the diagnostic.
    //
    // Measured 2026-08-07, and the reason this knob exists: the mass error
    // tracks the RE-SYNC FREQUENCY, not the pressure equation form.
    //   massConservativeP, thermoPerCorrector false (1 re-sync/outer): 0.91%
    //   massConservativeP, thermoPerCorrector true  (every corrector):  4.78%
    // Re-syncing more often throws the continuity density away more often.
    //
    // a < 1 lets the continuity density survive. The classic risk is unbounded
    // drift between the two densities, which is exactly what a = 1 was
    // suppressing -- so this is a trade, not a free win, and the sweep is the
    // experiment. Default 1 = unchanged behaviour.
    const scalar rhoResyncCoeff
    (
        pimple.dict().lookupOrDefault<scalar>("rhoResyncCoeff", scalar(1))
    );
    // rhoTransport (2026-08-07): Kawai Eq. (31) as a SEPARATE transported
    // continuity equation, so rho is no longer an EOS lookup at all. When it
    // is on, the EOS re-sync below must not run -- the transported density IS
    // the density. See the update block after the pressure correctors.
    const Switch rhoTransport
    (
        pimple.dict().lookupOrDefault<Switch>("rhoTransport", false)
    );
    if (rhoTransport)
    {
        // rho carries over from the previous corrector's Eq. (31) update.
    }
    else if (rhoResyncCoeff >= 1)
    {
        rho = thermo.rho();
    }
    else
    {
        rho = (1 - rhoResyncCoeff)*rho + rhoResyncCoeff*thermo.rho();
        Info<< "rhoResyncCoeff " << rhoResyncCoeff
            << ": |rho - rhoEOS| max = "
            << gMax(mag(rho - thermo.rho())().primitiveField())
            << " kg/m^3" << endl;
    }
    if (!rhoTransport)
    {
        rho.relax();
    }

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
    // rhofScheme (2026-08-07): make the FACE DENSITY consistent with the face
    // reconstruction that actually produced rho.
    //
    // At a passive contact (uniform p, uniform u) the exact answer is that rho
    // advects while p and u are untouched. The mass-conserving pressure
    // equation radiates nothing only if its two sources cancel discretely:
    //     fvc::ddt(rho)  +  fvc::div(phiHbyAm)  ==  0
    // But rho is NOT transported here -- it is slaved to (p, h, Z, C) through
    // the FGM manifold, and h/Z/C are convected with limitedLinear/
    // limitedLinear01, whereas rhof below defaults to PLAIN LINEAR. The same
    // density field therefore gets two different face reconstructions, the
    // cancellation fails by O(limiter difference x rho jump), and the residual
    // is a pressure source sitting exactly on the contact. That is the spike.
    //
    // The PEP form dodges this by deleting fvc::ddt(rho) altogether -- at the
    // cost of not conserving mass (measured: ~9% of throughput, 2026-08-06).
    // This switch attacks the cancellation instead, so mass conservation and
    // spike suppression stop being mutually exclusive.
    //
    // Set to the name of an entry in fvSchemes/interpolationSchemes, e.g.
    //     interpolationSchemes { rhoc  limitedLinear01 1; }
    //     PIMPLE            { rhofScheme rhoc; }
    // Empty (default) keeps the historical plain linear interpolate.
    // Supersedes rhofUpwind (which is the crude 1st-order member of this same
    // family, and which Stage-1 screening showed behaves as a smearer).
    const word rhofScheme
    (
        pimple.dict().lookupOrDefault<word>("rhofScheme", word::null)
    );
    tmp<surfaceScalarField> trhof;
    if (rhofUpwind)
    {
        trhof = upwind<scalar>(mesh, phi).interpolate(rho);
    }
    else if (rhofScheme != word::null)
    {
        trhof = fvc::interpolate(rho, phi, rhofScheme);
    }
    else
    {
        trhof = fvc::interpolate(rho);
    }
    const surfaceScalarField rhof("rhof", trhof);

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
    // rAUfDensityWeighted (2026-08-07): the momentum-weighted-interpolation
    // (MWI) recipe from the pressure-based LARGE-DENSITY-RATIO multiphase
    // literature, which routinely runs rho ratios of 1e6 -- three decades past
    // our 85x -- with the interfacial imbalance driven down to solver
    // tolerance (Bartholomew, Denner et al., JCP 375 (2018) 177; Denner & van
    // Wachem, Numer. Heat Transfer B 65 (2014) 218).
    //
    // The point: rAU = 1/A() ~ dt/rho JUMPS with the density, so ANY direct
    // face average of it (arithmetic or harmonic) is a choice between two bad
    // options at the contact. But the PRODUCT rho*rAU ~ dt is SMOOTH across
    // the same face. So interpolate the smooth product and divide by a
    // consistent face density:
    //     rAUf = interp(rho*rAU)/rhof = rhorAUf/rhof
    // which is dimensionally the same velocity-level mobility the pEqn
    // laplacian wants, but built from a quantity that has no jump to smear.
    // rhorAUf is already assembled above for the mass-level flux, so this
    // costs one division. Default off (the harmonic mean stays the baseline).
    const Switch rAUfDensityWeighted
    (
        pimple.dict().lookupOrDefault<Switch>("rAUfDensityWeighted", false)
    );
    const surfaceScalarField rAUf
    (
        "rAUf",
        rAUfDensityWeighted
      ? (rhorAUf/rhof)()
      : (1.0/fvc::interpolate(1.0/rAU))()
    );
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
    // pepTauS (2026-08-06): add the VISCOUS DISSIPATION tau:S that the
    // Terashima-Koshi pressure-evolution equation omits but Kawai Eq. (19)
    // carries (see the term itself below). NOT implied by pepFull -- kept an
    // independent switch so the A/B against the current production stack is
    // clean. Read every step (runTimeModifiable).
    const Switch pepTauS
    (
        pimple.dict().lookupOrDefault<Switch>("pepTauS", false)
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

    // --- massConservativeP: restore the STOCK mass-conserving pressure eqn ---
    // This module always replaces isothermalFluid's continuity-based pressure
    // equation with the volumetric pressure-evolution (PEP) form -- there was
    // no way back, and `pepFull` only ADDS terms on top of the replacement.
    // Measured consequence (2026-08-06): the domain leaks ~9% of its mass
    // throughput, systematically and without grid convergence, because no
    // discrete continuity equation is solved anywhere. Setting this switch
    // rebuilds the original
    //     fvc::ddt(rho) + psi*correction(fvm::ddt(p)) + fvc::div(phiHbyAm) = 0
    // with the MASS flux phiHbyAm = rhof*phiHbyAv and the mass-level face
    // mobility rhorAUf, giving a discretely mass-conserving solve. Expect the
    // transcritical contact to radiate its spurious pressure again -- that
    // trade is exactly what this switch is for measuring. Default off.
    const Switch massConservativeP
    (
        pimple.dict().lookupOrDefault<Switch>("massConservativeP", false)
    );
    const surfaceScalarField phiHbyAm("phiHbyAm", rhof*phiHbyAv);

    // Update the pressure BCs for flux consistency (3D: waveTransmissive
    // outlet, fixedFluxPressure walls). Volumetric-flux form: pass the
    // volumetric predicted flux and the rAUf (velocity-level) coefficient,
    // matching the pEqn's laplacian(rAUf, p). The rho overload is for MASS
    // flux; on a flow-carrying fixedFluxPressure patch it would subtract
    // rho_b*(Sf&U_b) from the volumetric phiHbyAv and skew snGrad(p) by
    // ~rho_b (LOX ~1000x). (backport RGP-13-GPU 2005bc9)
    if (massConservativeP)
    {
        constrainPressure(p, U, phiHbyAm, rhorAUf, MRF);
    }
    else
    {
        constrainPressure(p, U, phiHbyAv, rAUf, MRF);
    }

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
        // ladDriftSensor (2026-08-08): gate the artificial mass diffusion on
        // the TRANSPORTED-vs-EOS density discrepancy instead of on |grad rho|.
        //
        // Every LAD sensor tried so far fails for the same reason: it cannot
        // tell the PHYSICAL density jump from the SPURIOUS one. |grad rho| is
        // maximal at both. The `ladOddEven` pressure detector measured theta
        // max = 0.998 on this testbed -- it passes essentially every cell, so
        // LADrhoCoeff acted as global mass diffusion and blew up above ~0.05.
        //
        // The drift does not have that problem. Where the solution is right,
        // the transported and EOS densities AGREE, so the drift is zero -- at
        // the physical contact included. Measured 2026-08-07 on this testbed:
        // median |drift| = 0.002 kg/m^3 over 3744 cells, with 5 cells (0.1%)
        // reaching ~137 in an odd-even (+,-,+ on adjacent cells) pattern at
        // y ~ 0.47 mm, which is NOT on the contact (y ~ 0.60). A 5e4 contrast
        // between physical and spurious, versus none for |grad rho|.
        //
        // This is what Kawai's Eq. (31) diffusion is FOR; his sensor is a 4th
        // derivative with a Gaussian filter, unavailable on an unstructured
        // mesh. The drift is a better sensor anyway: it measures the very
        // inconsistency the diffusion exists to remove, so the loop is
        // self-extinguishing -- smoothing the drift turns the sensor off.
        // Requires rhoTransport (without it rho IS thermo.rho() and the drift
        // is identically zero). ladDriftEps = drift fraction of rho at which
        // the sensor saturates.
        const Switch ladDriftSensor
        (
            pimple.dict().lookupOrDefault<Switch>("ladDriftSensor", false)
        );
        const scalar ladDriftEps
        (
            pimple.dict().lookupOrDefault<scalar>("ladDriftEps", scalar(0.01))
        );
        // ladDucros (2026-08-08): the CONTACT-discontinuity sensor of
        // Jain, Agrawal & Moin, Phys. Rev. Fluids 9 (2024) 024609
        // (arXiv:2307.03257), Eq. (17):
        //
        //     f_D = |grad rho|^2 / ( |grad rho|^2
        //                          + a (theta^2 + omega_i omega_i)(rho/|u|)^2
        //                          + eps )
        //
        // theta = div(u) is the dilatation and omega = curl(u). It fires on a
        // density jump and switches OFF where dilatation (shocks) or enstrophy
        // (vortical motion) is large -- which is exactly the discrimination
        // every sensor in this solver has lacked. That paper states the problem
        // in our own words: a plain density/internal-energy indicator "will not
        // only detect contact discontinuities, but will also detect shocks and
        // vortical motions ... hence unnecessarily dissipative".
        //
        // Two properties make it usable here where Kawai's is not: it needs NO
        // Gaussian filtering of the artificial properties (the step that blocks
        // unstructured meshes), and it is validated with a SECOND-ORDER central
        // scheme at C_D = 0.5, not only with high-order compact differencing.
        //
        // Implemented multiplied through by |u|^2 so nothing divides by a
        // velocity that vanishes at a no-slip wall:
        //     f_D = A|u|^2 / ( A|u|^2 + a B rho^2 + eps ),  A=|grad rho|^2,
        //                                                   B=theta^2+|omega|^2
        // NOTE this swaps ONLY the sensor; the Dr envelope
        // (LADrhoCoeff*V^(1/3)*|U|) is left as-is so the drift-vs-fD comparison
        // varies one thing. Jain's own envelope uses (|u|+c_s) and an r-th
        // derivative of rho instead.
        const Switch ladDucros
        (
            pimple.dict().lookupOrDefault<Switch>("ladDucros", false)
        );
        const scalar ladDucrosA
        (
            pimple.dict().lookupOrDefault<scalar>("ladDucrosA", scalar(2))
        );
        if (ladDucros)
        {
            const scalarField A(magSqr(fvc::grad(rho))().primitiveField());
            const scalarField B
            (
                (sqr(fvc::div(U)) + magSqr(fvc::curl(U)))().primitiveField()
            );
            const scalarField U2(magSqr(U)().primitiveField());
            const scalarField R2(sqr(rho.primitiveField()));
            sensor = (A*U2)/(A*U2 + ladDucrosA*B*R2 + SMALL);
            label nAct = 0;
            forAll(sensor, celli)
            {
                if (sensor[celli] > 0.1) nAct++;
            }
            reduce(nAct, sumOp<label>());
            Info<< "LAD-rho ducrosContact: max = " << gMax(sensor)
                << ", cells active(>0.1) = " << nAct << endl;
        }
        if (ladDriftSensor && rhoTransport)
        {
            const scalarField drift
            (
                mag(rho.primitiveField() - thermo.rho()().primitiveField())
            );
            sensor =
                min(drift/(ladDriftEps*rho.primitiveField()), scalar(1));
            label nAct = 0;
            forAll(sensor, celli)
            {
                if (sensor[celli] > 0.1) nAct++;
            }
            reduce(nAct, sumOp<label>());
            Info<< "LAD-rho driftSensor: max = " << gMax(sensor)
                << ", cells active(>0.1) = " << nAct << endl;
        }
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

    // Kawai Eq. (31)-(32) pair: the artificial mass flux A_rho = Dr*grad(rho)
    // enters continuity as +div(A_rho) AND momentum as +div(A_rho (x) u). The
    // second half is what makes a uniform velocity field EXACTLY invariant
    // under the smoothing (constant-velocity condition) -- without it the
    // artificial mass diffusion generates spurious momentum at a contact.
    // Publish the face flux Sf & A_rho here so momentumPredictor (which runs
    // BEFORE this function in the PIMPLE loop) can pick it up; the one-
    // corrector lag is the same convention the manifold/EOS refresh already
    // uses. Zero whenever LADrhoCoeff = 0, so the term vanishes by itself.
    if (!mesh.foundObject<surfaceScalarField>("phiA_amd"))
    {
        mesh.objectRegistry::store
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phiA_amd",
                    mesh.time().name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimMass/dimTime, 0)
            )
        );
    }
    mesh.lookupObjectRef<surfaceScalarField>("phiA_amd") =
        fvc::interpolate(Dr)*fvc::snGrad(rho)*mesh.magSf();

    fvScalarMatrix pDDtEqn
    (
        massConservativeP
      ? fvScalarMatrix
        (
            fvc::ddt(rho)
          + psi*correction(fvm::ddt(p))
          + fvc::div(phiHbyAm)
        )
      : rhoTransport
        // Kawai-faithful split: Eq. (33) carries NO A_rho term -- the
        // artificial mass diffusion lives in Eqs. (31)-(32) only. With
        // rhoTransport on, continuity (31) already applies div(Dr grad rho);
        // keeping the same term here DOUBLE-COUNTS it: rho smooths while the
        // pressure equation receives the diffusion again as a separate
        // volumetric source. Measured 2026-08-08: with the term in both
        // places, ANY Dr > 0 (C = 0.5, either sensor) pinned |U| at the 400
        // clamp from t = 1.8e-7 onward, monotonically worse with C.
      ? fvScalarMatrix
        (
            psis*fvm::ddt(p)
          + fvc::div(phiHbyAv)
        )
      : fvScalarMatrix
        (
            psis*fvm::ddt(p)
          + fvc::div(phiHbyAv)
          - fvc::laplacian(Dr, rho)/rho
        )
    );

    // --- pepRealCoeff: REAL-FLUID internal-energy source coefficient --------
    // Kawai Eq. (33) carries the (tau:S - div q) source with the coefficient
    // (1/rho)*alpha_p/(c_v*beta_t), where alpha_p = -(1/rho)(d rho/dT)_p is the
    // thermal expansivity and beta_t = (1/rho)(d rho/dp)_T = psi/rho the
    // isothermal compressibility. Our pepRHS/pepTauS used (gamma - 1), which is
    // exactly that group's IDEAL-GAS limit and is wrong by orders of magnitude
    // in dense cold LOx. alpha_p is not exposed by the thermo interface, but
    // the standard relation  c_p - c_v = T*alpha_p^2/(rho*beta_t)  inverts it
    // with only Cp, Cv, T and psi (note rho*beta_t = psi):
    //     alpha_p = sqrt(psi*(Cp - Cv)/T)
    //     X = alpha_p/(rho*c_v*beta_t) = sqrt((Cp - Cv)/(T*psi))/Cv
    // Ideal-gas check: Cp - Cv = R, T*psi = 1/R  ->  X = R/Cv = gamma - 1.
    // Default off so the A/B against the ideal-gas coefficient stays clean.
    const Switch pepRealCoeff
    (
        pimple.dict().lookupOrDefault<Switch>("pepRealCoeff", false)
    );
    tmp<volScalarField> tExpCoeff;
    if (pepRealCoeff)
    {
        const volScalarField& Cp = thermo.Cp();
        const volScalarField& Cv = thermo.Cv();
        // Cp - Cv >= 0 thermodynamically; a table-interpolated pair can dip
        // marginally negative, which would make the sqrt NaN. Floor at zero.
        tExpCoeff =
            sqrt(max(Cp - Cv, dimensionedScalar(Cp.dimensions(), 0))
                /(thermo.T()*thermo.psi()))/Cv;
        Info<< "pepRealCoeff: X min/max = "
            << gMin(tExpCoeff().primitiveField()) << " / "
            << gMax(tExpCoeff().primitiveField())
            << "  (ideal-gas (gamma-1) min/max = "
            << gMin(thermo.gamma()().primitiveField()) - 1 << " / "
            << gMax(thermo.gamma()().primitiveField()) - 1 << ")" << endl;
    }
    else
    {
        tExpCoeff = thermo.gamma() - dimensionedScalar(dimless, 1);
    }
    const volScalarField& expCoeff = tExpCoeff();

    // --- pepFull: pressure ADVECTION + energy-diffusion RHS ----------------
    // (1) psis*(u.grad p) in exact flux form: div(F p) - p div(F) = F.grad(p)
    //     with F = psisf * (volumetric flux). Uses the CURRENT mass flux
    //     phi/rhof (last corrector) as the advecting velocity.
    // (2) RHS (gamma-1)*psis*laplacian(Dh, h): the diffusive expansion source
    //     of the full pressure-evolution equation, discretised with the SAME
    //     Deff("h") as the h transport for consistency.
    if (!massConservativeP && pepAdvect)
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
    if (!massConservativeP && pepRHS)
    {
        if (fgmTable_.useEnthalpy())
        {
            const volScalarField& h = hPtr_();
            const volScalarField DhP("DhP", Deff("h"));
            pDDtEqn -= expCoeff*psis*fvc::laplacian(DhP, h);
        }
    }
    // --- pepTauS: VISCOUS DISSIPATION tau:S in the pressure-evolution RHS ---
    // Kawai, Terashima & Negishi JCP 300 (2015) Eq. (19) carries the internal-
    // energy source as (tau:S - div q); our pepFull (Terashima-Koshi JCP 231
    // (2012)) supplied only the -div q half (the Deff-grad-h laplacian above).
    // Kawai Sec. 2.3/6.3 identify exactly this tau:S inconsistency in Refs
    // [16,17] -- "non-negligible spurious impact on the thermodynamic
    // properties" -- and the original authors issued a corrigendum (JCP 283
    // (2015) 609). tau:S is the double inner product of the viscous stress and
    // the strain rate = viscous heating [W/m^3], >= 0, and scales as mu*S^2:
    // negligible in the bulk but 3-4 decades larger inside wall prism layers
    // (du/dy over a 3 um first cell) -- i.e. exactly where the injector's dt
    // bottleneck sits. Default off so the A/B is clean.
    if (!massConservativeP && pepTauS)
    {
        const volScalarField muEff
        (
            "muEff", thermo.mu() + rho*momentumTransport().nut()
        );
        const volSymmTensorField Sr(symm(fvc::grad(U)));
        const volScalarField tauS
        (
            "tauS",
            2.0*muEff*(Sr && Sr)
          - (2.0/3.0)*muEff*sqr(fvc::div(U))
        );
        pDDtEqn -= expCoeff*psis*tauS;
        Info<< "pepTauS: tau:S max = " << gMax(tauS.primitiveField())
            << " W/m^3" << endl;
    }

    // Snapshot p before the solve so massConservativeP can apply the stock
    // SIMPLErho density increment below (see the block after the clamp).
    autoPtr<volScalarField> pPre;
    if (massConservativeP)
    {
        pPre.set(new volScalarField("pPre", p));
    }

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            pDDtEqn
          - fvm::laplacian(massConservativeP ? rhorAUf : rAUf, p)
        );

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
            // (mass-conservative branch already solves in mass units).
            phi =
                massConservativeP
              ? (phiHbyAm + pEqn.flux())
              : rhof*(phiHbyAv + pEqn.flux());
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

    // massConservativeP: the stock pressure equation carries fvc::ddt(rho),
    // so its algorithm REQUIRES the matching rho += psi*dp increment after the
    // solve -- otherwise the density the next ddt(rho) differences is not the
    // density the pressure equation just solved for, and the mismatch is a
    // spurious source sitting on the contact. Running the stock pDDtEqn with
    // the PEP path's disabled increment is an inconsistent hybrid; this closes
    // it. Guarded by mcpCorrectRho (default true when massConservativeP is on)
    // because the increment is exactly what blew rho negative in the flame-zone
    // spike that got it disabled for PEP: psi_gas ~1e-5 x dp ~ -4e8 -> drho
    // ~ -4000. Floored at a small positive density for that reason.
    if (massConservativeP && pPre.valid())
    {
        const Switch mcpCorrectRho
        (
            pimple.dict().lookupOrDefault<Switch>("mcpCorrectRho", true)
        );
        if (mcpCorrectRho)
        {
            rho = max
            (
                rho + psi*(p - pPre()),
                dimensionedScalar(dimDensity, 1e-3)
            );
            rho.correctBoundaryConditions();
            Info<< "mcpCorrectRho: rho min/max = "
                << gMin(rho.primitiveField()) << " / "
                << gMax(rho.primitiveField()) << " kg/m^3" << endl;
        }
    }

    // Thermodynamic density update: the stock SIMPLErho increment
    // correctRho(psi*dp) is DISABLED in the PEP path. At a flame-zone
    // pressure spike the increment is huge (psi_gas ~1e-5 x dp ~ -4e8 ->
    // drho ~ -4000) and drove the THERMO density to -4869 kg/m3 in written
    // states -- and the per-corrector updateManifold()+thermo_.correct()+
    // rho_ re-sync at the top of this function already provides the full
    // EOS-consistent density update (one-corrector lag, nOuter >= 3).

    // --- rhoTransport: Kawai Eq. (31) --------------------------------------
    //     d(rho)/dt + div(rho u) = div(A_rho),     A_rho = Dr grad(rho)
    // solved with the mass flux phi the pressure correctors just produced.
    //
    // Why this is the piece that was missing. The two routes tried before were
    // treated as an either/or:
    //   - PEP pressure equation (33): well conditioned, but NO continuity
    //     equation is solved anywhere -> ~9% of throughput leaks (2026-08-06).
    //   - stock mass-conserving pressure equation: conserves mass, but its
    //     fvc::ddt(rho) source is O(8e9) and must cancel div(rho u) to 2e-5
    //     relative or the residual is amplified by 1/psi ~ 2.5e5 into bar-level
    //     pressure. Four attempts to tame it (face-scheme matching, the
    //     psi*dp increment, per-corrector thermo, resync blending) all failed.
    // Kawai has BOTH because his continuity equation is NOT his pressure
    // equation. A segregated solver only needs its pressure equation to
    // constrain div(u) -- (33) does that -- leaving (31) free to stand alone.
    //
    // The update is diagonal-in-time and explicit in the flux, i.e.
    //     rho^{n+1} = rho^n - dt*div(phi) + dt*div(Dr grad rho)
    // so summing over cells telescopes div(phi) to the boundary fluxes: total
    // mass conservation is an IDENTITY here, not an approximation. Written
    // directly rather than via solve() because rho's boundary condition is
    // `calculated` and would not evaluate in a matrix solve; boundary values
    // are taken from the EOS, where the BCs fix the state anyway.
    //
    // rho now differs from thermo.rho(). That drift is the quasi-conservative
    // inconsistency, relocated from mass (where it was 9%) to the
    // thermodynamic state (where it is measured below). No relaxation term is
    // applied: the PIMPLE outer loop is itself the fixed-point iteration, and
    // whether that map contracts is the experiment.
    if (rhoTransport)
    {
        const scalar dt = mesh.time().deltaTValue();
        tmp<volScalarField> tdiv(fvc::div(phi));
        if (LADrhoCoeff > 0)
        {
            tdiv.ref() -= fvc::laplacian(Dr, rho);
        }
        rho.primitiveFieldRef() =
            rho.oldTime().primitiveField() - dt*tdiv().primitiveField();
        rho.boundaryFieldRef() == thermo.rho()().boundaryField();
        // Diagnostic field: the transported-vs-EOS density discrepancy, i.e.
        // the quasi-conservative inconsistency this formulation relocates from
        // mass into the thermodynamic state. Registered AUTO_WRITE so the
        // SPATIAL distribution can be inspected -- a max value alone cannot
        // tell an interface-local artefact from a domain-wide bias.
        if (!mesh.foundObject<volScalarField>("rhoDrift"))
        {
            mesh.objectRegistry::store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "rhoDrift",
                        mesh.time().name(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    rho - thermo.rho()
                )
            );
        }
        else
        {
            mesh.lookupObjectRef<volScalarField>("rhoDrift") =
                rho - thermo.rho();
        }
        // Interval-consistent mass diagnostic, measured AT THE UPDATE POINT
        // with the SAME phi the update used -- immune to the two artifacts
        // found 2026-08-08: (i) the offline budget compares interval-averaged
        // dM/dt against WRITE-INSTANT fluxes, which the adjustableRunTime
        // landing step biases; (ii) the thermophysicalPredictor print pairs
        // ddt(rho) with a phi that has since been re-solved. Here
        // dM = sum((rho - rho.old)*V)/dt and F = boundary sum of phi are the
        // exact pair the update telescoped, so their mismatch is the TRUE
        // per-step conservation error of this formulation.
        {
            scalar dM =
                gSum
                (
                    (rho.primitiveField() - rho.oldTime().primitiveField())
                   *mesh.V()
                )/mesh.time().deltaTValue();
            // Physical patches only -- processor patches differ in number
            // per rank, and a collective (gSum) inside this loop deadlocks.
            // Local sums, then ONE reduce.
            scalar F = 0;
            forAll(phi.boundaryField(), patchi)
            {
                if (!mesh.boundary()[patchi].coupled())
                {
                    F += sum(phi.boundaryField()[patchi]);
                }
            }
            reduce(F, sumOp<scalar>());
            Info<< "rhoTransport mass: dM/dt = " << dM
                << "  bFlux = " << F
                << "  err = " << dM + F << " kg/s" << endl;
        }
        // rhoAnchorCoeff (2026-08-08): WEAK EOS anchor after the transport
        // update,  rho <- (1-a)*rho + a*thermo.rho(),  a << 1.
        //
        // Motivation: full transport (rhoTransportFull) closes the mass budget
        // (0.3% of throughput) but kills the physical 10 kHz vortex-shedding
        // tone at probe P1 (RMSu1 0.63 -> 0.107; the mode and its 20 kHz
        // harmonic vanish from the spectrum), while the accidental hybrid
        // (per-corrector full wipe) preserves the shedding but muddles the
        // budget. A weak anchor interpolates between them: the anchoring
        // timescale is dt/a, so with a ~ 0.05 the density is pulled to the
        // EOS over ~20 steps -- far shorter than the ~3000-step shedding
        // period (so the instability dynamics keep their EOS-consistent
        // density) yet far longer than one step (so the within-step continuity
        // increment survives and the budget stays closed to O(a)).
        // a = 0 keeps pure transport; a = 1 reproduces the full wipe.
        const scalar rhoAnchorCoeff
        (
            pimple.dict().lookupOrDefault<scalar>("rhoAnchorCoeff", scalar(0))
        );
        if (rhoAnchorCoeff > 0)
        {
            rho =
                (1 - rhoAnchorCoeff)*rho + rhoAnchorCoeff*thermo.rho();
        }
        Info<< "rhoTransport: rho min/max = "
            << gMin(rho.primitiveField()) << " / "
            << gMax(rho.primitiveField())
            << "  |rho - rhoEOS| max = "
            << gMax(mag(rho - thermo.rho())().primitiveField())
            << " kg/m^3" << endl;
    }

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
    //
    // rhoTransportFull (2026-08-08): this is the SECOND EOS re-sync in this
    // function (the comment above the first one even warns it used to run
    // twice). The rhoTransport update happens BEFORE this point, so with
    // SIMPLErho yes this line wiped the transported density at the end of
    // every corrector -- measured: with the AMD double-count fixed, runs with
    // Dr = 0 and Dr > 0 were BIT-IDENTICAL in p, U and written rho; only the
    // pre-wipe diagnostics differed. What survived (via rho.relax() 0.7 on
    // non-final correctors) was a within-step 30% blend -- an accidental
    // HYBRID: EOS-anchored each step, continuity-corrected within correctors.
    // That hybrid is what produced the s6/s9 improvements (mass 5x, pHi -31%).
    // rhoTransportFull skips the wipe entirely, making the transport genuine
    // and letting drift accumulate across steps -- the A/B against the hybrid
    // decides whether the EOS anchor was load-bearing.
    const Switch rhoTransportFull
    (
        pimple.dict().lookupOrDefault<Switch>("rhoTransportFull", false)
    );
    if (pimple.simpleRho() && !(rhoTransport && rhoTransportFull))
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
