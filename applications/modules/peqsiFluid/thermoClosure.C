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

\*---------------------------------------------------------------------------*/

#include "peqsiFluid.H"
#include "tabulatedRealGasMixture.H"
#include "thermodynamicConstants.H"
#include "fvcLaplacian.H"
#include "fvcGrad.H"
#include <algorithm>
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::peqsiFluid::updateCoefficients()
{
    // ------------------------------------------------------------------
    // Exact real-gas coefficients (WKK Eqs. 9, 11, 12) assembled from the
    // validated thermo outputs (cp, cv, psi = (drho/dp)_T, rho, T) via
    // thermodynamic identities -- no new EoS-level code:
    //
    //   (dp/dv)_T = -rho^2/psi
    //   (dp/dT)_v = rho sqrt( (cp - cv)/(T psi) )
    //               [ cp - cv = -T (dp/dT)_v^2 (dv/dp)_T; positive root ]
    //   xi        = (dh/dT)_v = cv + v (dp/dT)_v
    //   (dh/dv)_T = T (dp/dT)_v + v (dp/dv)_T
    //               [ (du/dv)_T = T (dp/dT)_v - p ]
    //
    //   alpha = (dp/dT)_v / (rho xi)                        (WKK Eq. 11)
    //   beta  = [ (dp/dv)_T - (dp/dT)_v (dh/dv)_T / xi ]/rho  (Eq. 12)
    //
    // Cell-wise consistency check: beta/(1-alpha) == -rho c^2 = -gamma/psi
    // (WKK App. D).  Pointwise evaluation; zero-gradient boundaries.
    // ------------------------------------------------------------------

    const auto tCp = thermo_.Cp();
    const auto tCv = thermo_.Cv();

    const scalarField& cp = tCp().primitiveField();
    cpPrev_ = cp;   // reused by the LAD stability cap next step
    const scalarField& cv = tCv().primitiveField();
    const scalarField& psi = thermo_.psi().primitiveField();
    const scalarField& T = thermo_.T().primitiveField();
    const scalarField& rho = rho_.primitiveField();

    scalarField& a = alpha_.primitiveFieldRef();
    scalarField& b = beta_.primitiveFieldRef();

    scalar maxDev = 0;

    // The sweep's v = 1/rho and every derivative behind it are
    // unphysical when the transported rho has drifted below a physical
    // floor (the lip episodes reach rho < 0): a sign-flipped alpha/beta
    // feeds the Helmholtz diagonal with the wrong sign, the operator
    // goes indefinite and one solve detonates.  Those cells get the
    // coefficients of the EOS state (p, T) instead.  The bad count is
    // reduced on every rank BEFORE thermo_.rho() so the fetch stays
    // rank-uniform.
    const scalar rhoFloor
    (
        pimple.dict().lookupOrDefault<scalar>("peqsiRhoCoeffFloor", 0.01)
    );
    // Upper wall, symmetric with the floor: SRK leaves its domain at
    // the covolume, v = b_mix, and past it p(v, T) flips sign (measured:
    // 3.5% over gives -4.7 GPa -- the -GPa dp signature of the lip
    // events).  The wall is temperature-independent and exact in
    // composition: rho_wall = W/bM, i.e. the violation test is
    // rho * bM / W > 1.  bM comes from the Tier-2 fields the closure
    // already fills per cell; with a 2% standoff so the sweep never
    // evaluates a state on the sign flip itself.  Cells past the wall
    // get the EOS-state coefficients, exactly like the floor cells.
    const Switch covolGuard =
        pimple.dict().lookupOrDefault<Switch>("peqsiCovolGuard", true);

    // Covolume wall from the TRANSPORTED composition.  The wall is
    // rho_wall = 1/sum_k(Y_k b_k) with per-species mass covolumes b_k
    // (SRK+Peneloux, bisected on the actual p sign flip; table session
    // 2026-08-26).  The table-bM route was measured to MISS the one
    // confirmed real violation: an off-manifold fuel-gallery cell
    // (3.1% H2) whose true wall (704) the manifold composition
    // overestimates by 25% -- rho is a transported quantity, so the
    // composition the wall is judged on must be transported too.
    static bool bkBuilt = false;
    static labelList bkIdx;
    static scalarList bkVal;
    if (covolGuard && !bkBuilt)
    {
        const HashTable<scalar> bk
        ({
            {"O2", 6.713045e-04}, {"CO", 9.433968e-04},
            {"NC10H22", 1.156827e-03}, {"CO2", 5.989275e-04},
            {"CH4", 1.818376e-03}, {"H2O", 9.802690e-04},
            {"C6H6", 9.470105e-04}, {"PHC3H7", 1.004178e-03},
            {"C2H2", 1.466095e-03}, {"C2H4", 1.514882e-03},
            {"CYC9H18", 1.052476e-03}, {"TOLUEN", 9.835832e-04},
            {"OH", 6.542857e-04}, {"C3H6", 1.369163e-03},
            {"H2", 9.783274e-03}, {"C4H6", 1.253522e-03},
            {"PC3H4", 1.299461e-03}, {"O", 6.955082e-04},
            {"STYREN", 9.651829e-04}, {"C2H6", 1.413318e-03},
            {"C3H3", 1.175729e-03}, {"CH2CO", 9.066221e-04},
            {"CH2O", 1.324750e-03}, {"AC7H14", 1.211563e-03},
            {"C4H8", 1.291707e-03}, {"AC3H5", 1.267570e-03},
            {"AC5H10", 1.257264e-03}, {"C4H2", 1.477722e-03},
            {"IC4H3", 1.448554e-03}, {"DC10H21", 1.165081e-03},
            // raw SRK (no Peneloux, +2.8% -> fires slightly early,
            // the safe direction) for the default specie
            {"N2", 9.560000e-04}
        });
        const PtrList<volScalarField>& Ysp =
            thermo_.Y();
        DynamicList<label> di;
        DynamicList<scalar> dv;
        forAll(Ysp, k)
        {
            if (bk.found(Ysp[k].name()))
            {
                di.append(k);
                dv.append(bk[Ysp[k].name()]);
            }
        }
        bkIdx.transfer(di);
        bkVal.transfer(dv);
        Info<< "PEQSI covolume guard: wall from transported Y over "
            << bkIdx.size() << " species (+ backstop rho > 1670)"
            << endl;
        bkBuilt = true;
    }

    scalarField sumYb;
    if (covolGuard && bkIdx.size())
    {
        sumYb.setSize(rho.size(), 0.0);
        const PtrList<volScalarField>& Ysp =
            thermo_.Y();
        forAll(bkIdx, j)
        {
            const scalarField& Yk = Ysp[bkIdx[j]].primitiveField();
            const scalar bkj = bkVal[j];
            forAll(sumYb, i) sumYb[i] += bkj*Yk[i];
        }
    }

    // Threshold 0.95, not 0.98: the pressure response is so steep near
    // the wall that margin 0.97 already means ~820 MPa (156x chamber)
    // in cold O2 -- while legitimate 60 K LOX sits at margin 0.883, so
    // the usable window is [0.89, 0.97] and 0.95 (413 MPa, before the
    // sign flip; 7% above normal LOX) is the working point.  This is a
    // LAST LINE, not a fix: a cell at 0.95 is already physically wrong,
    // the guard only prevents the FPE/e18 runaway and buys time for
    // the real cure (front-normalisation of the rho/h undershoot).
    const scalar covolMargin =
        pimple.dict().lookupOrDefault<scalar>("peqsiCovolMargin", 0.95);
    auto pastWall = [&](const label i) -> bool
    {
        if (!covolGuard) return false;
        // composition-blind backstop: no mixture's wall exceeds
        // 1669.7 kg/m3 (pure CO2)
        if (rho[i] > 1670.0) return true;
        return sumYb.size() && rho[i]*sumYb[i] > covolMargin;
    };

    // Consistency band (peqsiDriftGuardLo/Hi, default 0 = off): the
    // floor catches rho < 0.01 and the wall catches rho near 1/sum(Yb),
    // but the whole interval between them is unguarded for CONSISTENCY
    // -- the measured worst error (rho 38.6 vs rho_EOS 854.8, 22x) sits
    // exactly there and passes both.  When armed, a cell whose
    // rho/rho_EOS leaves [lo, hi] gets the EOS-state coefficients like
    // the other two.  Thresholds must come from the quiet-phase drift
    // histogram (front cells carry large chronic drift legitimately);
    // armed only when both knobs are set.
    const scalar driftLo =
        pimple.dict().lookupOrDefault<scalar>("peqsiDriftGuardLo", 0.0);
    const scalar driftHi =
        pimple.dict().lookupOrDefault<scalar>("peqsiDriftGuardHi", 0.0);
    const bool driftGuard = (driftLo > 0 && driftHi > driftLo);
    tmp<volScalarField> trhoEg;
    const scalarField* rhoEg = nullptr;
    if (driftGuard)
    {
        trhoEg = thermo_.rho();
        rhoEg = &trhoEg().primitiveField();
    }
    auto offBand = [&](const label i) -> bool
    {
        if (!driftGuard) return false;
        const scalar re = (*rhoEg)[i];
        if (re <= vSmall) return false;   // wall guard's territory
        const scalar r = rho[i]/re;
        return r < driftLo || r > driftHi;
    };

    // Rate guard (peqsiRateGuard, default 0 = off): the per-step
    // |drho|/rho separates the lethal transient (measured 30%/step at
    // the blowing cell) from chronic front drift (quiet-phase 500-step
    // accumulations top out at 61%, i.e. ~0.1%/step) by 2-3 orders of
    // magnitude -- the one signal where the two populations do not
    // overlap, unlike any level-based band.  A cell moving faster than
    // the threshold gets the EOS-state coefficients for that step.
    const scalar rateGuard =
        pimple.dict().lookupOrDefault<scalar>("peqsiRateGuard", 0.0);
    // Observe-only census (peqsiRateCensus): per-step |drho|/rho
    // quantiles with no intervention -- the measurement that sets the
    // guard threshold and, against the physical bound U dt/delta =
    // 0.1-3 %/step, judges whether the front scheme itself is inside
    // its legitimate range.
    const Switch rateCensus =
        pimple.dict().lookupOrDefault<Switch>("peqsiRateCensus", false);
    const bool haveRate =
        rateGuard > 0
     && rhoPrevGuard_.valid()
     && rhoPrevGuard_().size() == rho.size();
    if
    (
        rateCensus
     && rhoPrevGuard_.valid()
     && rhoPrevGuard_().size() == rho.size()
    )
    {
        const label every =
            pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
        static label nC = 0;
        if (every > 0 && (nC++ % every) == 0)
        {
            // quantiles via a fixed log histogram (64 bins, 1e-8..1)
            labelList bins(64, label(0));
            const scalarField& rp = rhoPrevGuard_();
            forAll(rho, i)
            {
                const scalar r =
                    mag(rho[i] - rp[i])/max(mag(rp[i]), scalar(1e-3));
                label b =
                    r <= 1e-8 ? 0
                  : min(label(64.0*(log10(r) + 8.0)/9.0), label(63));
                bins[b]++;
            }
            reduce(bins, sumOp<labelList>());
            label total = 0;
            forAll(bins, b) total += bins[b];
            auto quant = [&](const scalar q) -> scalar
            {
                label target = label(q*total), c = 0;
                forAll(bins, b)
                {
                    c += bins[b];
                    if (c >= target)
                    {
                        return pow(10.0, -8.0 + 9.0*(b + 1)/64.0);
                    }
                }
                return 1.0;
            };
            Info<< "PEQSI rate census: |drho|/rho per step p50 < "
                << quant(0.50) << ", p99 < " << quant(0.99)
                << ", p99.99 < " << quant(0.9999)
                << " (log-bin upper bounds)" << endl;
        }
    }
    if
    (
        (rateGuard > 0 || rateCensus)
     && (!rhoPrevGuard_.valid()
      || rhoPrevGuard_().size() != rho.size())
    )
    {
        rhoPrevGuard_.reset(new scalarField(rho));
    }
    auto tooFast = [&](const label i) -> bool
    {
        if (!haveRate) return false;
        const scalar rp = rhoPrevGuard_()[i];
        return mag(rho[i] - rp) > rateGuard*max(mag(rp), scalar(1e-3));
    };

    label nBad = 0;
    label nWallFired = 0;
    label nBandFired = 0;
    label nRateFired = 0;
    scalar wallSumYbMin = great, wallSumYbMax = -great;
    forAll(rho, i)
    {
        if (pastWall(i))
        {
            nBad++;
            nWallFired++;
            if (sumYb.size())
            {
                wallSumYbMin = min(wallSumYbMin, sumYb[i]);
                wallSumYbMax = max(wallSumYbMax, sumYb[i]);
            }
        }
        else if (rho[i] < rhoFloor) nBad++;
        else if (offBand(i)) { nBad++; nBandFired++; }
        else if (tooFast(i)) { nBad++; nRateFired++; }
    }
    if ((rateGuard > 0 || rateCensus) && rhoPrevGuard_.valid())
    {
        rhoPrevGuard_() = rho;
    }

    {
        // Fired-cell audit (table session: a corrupted bM makes the
        // guard fire LATE -- log the bM range of fired cells so a
        // late fire is identifiable after the fact)
        const label every =
            pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
        static label n = 0;
        label nW = nWallFired;
        reduce(nW, sumOp<label>());
        if (every > 0 && (n++ % every) == 0 && nW > 0)
        {
            reduce(wallSumYbMin, minOp<scalar>());
            reduce(wallSumYbMax, maxOp<scalar>());
            Info<< "PEQSI covolume guard: " << nW
                << " cells past the wall, fired-cell sumYb in ["
                << wallSumYbMin << ", " << wallSumYbMax << "]" << endl;
        }
        label nB = nBandFired;
        reduce(nB, sumOp<label>());
        if (every > 0 && (n % every) == 0 && nB > 0)
        {
            Info<< "PEQSI drift-band guard: " << nB
                << " cells outside [" << driftLo << ", " << driftHi
                << "] of rho_EOS" << endl;
        }
        label nR = nRateFired;
        reduce(nR, sumOp<label>());
        if (every > 0 && (n % every) == 0 && nR > 0)
        {
            Info<< "PEQSI rate guard: " << nR
                << " cells moved more than " << rateGuard
                << " of their rho in one step" << endl;
        }
    }

    // Crossing-speed diagnostic: the margin of the worst fired cell on
    // the PREVIOUS step tells how fast the undershoot travels from the
    // legitimate 0.883 toward the wall -- the target number for the
    // front-normalisation work.  If a cell crosses in a single step the
    // guard itself cannot catch it in time.
    if (covolGuard && sumYb.size())
    {
        label worstI = -1;
        scalar worstM = 0;
        forAll(rho, i)
        {
            const scalar m = rho[i]*sumYb[i];
            if (m > covolMargin && m > worstM)
            {
                worstM = m;
                worstI = i;
            }
        }
        if (worstI >= 0)
        {
            const scalar mPrev =
                (covolMarginPrev_.valid()
              && covolMarginPrev_().size() == rho.size())
              ? covolMarginPrev_()[worstI] : -1.0;
            Pout<< "PEQSI covolume crossing: margin "
                << mPrev << " -> " << worstM
                << " in one step at " << mesh.C()[worstI] << endl;
        }
        if
        (
            !covolMarginPrev_.valid()
         || covolMarginPrev_().size() != rho.size()
        )
        {
            covolMarginPrev_.reset(new scalarField(rho.size(), 0.0));
        }
        scalarField& mp = covolMarginPrev_();
        forAll(rho, i) mp[i] = rho[i]*sumYb[i];
    }
    tmp<volScalarField> trhoE;
    if (returnReduce(nBad, sumOp<label>()) > 0)
    {
        trhoE = thermo_.rho();
    }
    const scalarField* rhoEPtr =
        trhoE.valid() ? &trhoE().primitiveField() : nullptr;

    forAll(rho, i)
    {
        const scalar rhoC =
            (rho[i] < rhoFloor || pastWall(i) || offBand(i) || tooFast(i))
          ? max((*rhoEPtr)[i], rhoFloor)
          : rho[i];
        const scalar v = 1.0/rhoC;
        const scalar dpdv = -sqr(rhoC)/psi[i];
        const scalar dpdT = rhoC*sqrt(max((cp[i] - cv[i])/(T[i]*psi[i]), 0.0));
        const scalar xi = cv[i] + v*dpdT;
        const scalar dhdv = T[i]*dpdT + v*dpdv;

        a[i] = dpdT/(rhoC*xi);
        b[i] = (dpdv - dpdT*dhdv/xi)/rhoC;

        // App. D identity: beta/(1-alpha) == -rho c^2, c^2 = gamma/psi
        // (gamma = cp/cv from the already-evaluated fields -- the
        // thermo.gamma() call re-evaluates BOTH Cp and Cv sweeps)
        const scalar lhs = b[i]/(1.0 - a[i]);
        const scalar rhs = -rhoC*(cp[i]/cv[i])/psi[i];
        maxDev = max(maxDev, mag(lhs - rhs)/max(mag(rhs), small));
    }

    // ------------------------------------------------------------------
    // Manifold coefficients (peqsiUseTableCoeffs).  The runtime SRK
    // evaluation above stays as the fallback: it serves every cell the
    // manifold cannot -- cells whose specific volume has drifted off the
    // table state, and every cell at all when the switch is off.
    //
    // The motivation is cost, not accuracy.  Profiled on the
    // multicomponent path, the temperature Newton and the property
    // sweeps behind it take 84% of the step against 2.1% for a single
    // species, because thirty-species janaf and SRK departure functions
    // are evaluated per cell per iteration.  (An earlier claim that the
    // tabulated derivatives would also be more accurate was withdrawn:
    // the (cp - cv) identity route reproduces the analytic derivatives
    // to six figures.)
    // ------------------------------------------------------------------
    if
    (
        fgmActive_
     && alphaTab_.valid()
     && pimple.dict().lookupOrDefault<Switch>("peqsiUseTableCoeffs", false)
    )
    {
        const scalarField& aT = alphaTab_();
        const scalarField& bT = betaTab_();
        const boolList& useT = tabUsable_();

        label nTab = 0;
        forAll(a, i)
        {
            if (useT[i])
            {
                a[i] = aT[i];
                b[i] = bT[i];
                nTab++;
            }
        }
        {
            const label every =
                pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
            static label n = 0;
            if (every > 0 && (n++ % every) == 0)
            {
                reduce(nTab, sumOp<label>());
                label nAll = a.size();
                reduce(nAll, sumOp<label>());
                Info<< "PEQSI coefficients: from manifold on " << nTab
                    << "/" << nAll << " cells (runtime SRK elsewhere)"
                    << endl;
            }
        }

        alpha_.correctBoundaryConditions();
        beta_.correctBoundaryConditions();
    }

    // Interval-gated: the identity is structural (assembled from the
    // same fields it is checked against), so a violation is a CODE
    // regression, not a flow event -- catching it within diagInterval
    // steps loses nothing, and the reduction is a per-step sync point.
    {
        const label every =
            pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
        static label n = 0;
        if (every > 0 && (n++ % every) == 0)
        {
            reduce(maxDev, maxOp<scalar>());
            Info<< "PEQSI coefficients: "
                << "max |beta/(1-alpha) + rho c^2| / (rho c^2) = "
                << maxDev << endl;
        }
    }

    // What the density drift actually costs (peqsiCoeffDrift).
    //
    // The coefficients above mix two determinations of the same state:
    // rho is the TRANSPORTED density, while cp, cv and psi come from the
    // thermo, evaluated at the EOS state (p, T, Y).  Where the two
    // densities disagree -- the drift -- alpha and beta are assembled at
    // a state that is not a state.
    //
    // Recomputing both ways answers what that costs, and the answer is
    // nothing, for a reason worth recording.
    //
    //   measured: rho gap 0.257 max -> alpha 4.8e-16, beta 0.257
    //
    // alpha is invariant to machine precision and beta tracks the gap
    // exactly 1:1.  Both follow from the algebra: with
    // (dp/dv)_T = -rho^2/psi and (dp/dT)_v = rho S, S = sqrt((cp-cv)/(T psi)),
    //
    //   xi   = cv + v (dp/dT)_v = cv + S          (v rho = 1, so no rho)
    //   alpha = (dp/dT)_v/(rho xi) = S/xi         -> rho cancels
    //   beta  = rho [ -1/psi - S (T S - 1/psi)/xi ]  -> beta is LINEAR in rho
    //
    // (An earlier expectation that beta would amplify the gap by its
    // -7.6 to -14.8 v-exponent was wrong: that exponent describes beta
    // along the MANIFOLD, where composition and state both move, not
    // this linear prefactor.)
    //
    // Because beta is linear in rho, and every dynamical use of it is in
    // the combination rho/beta or beta/rho --
    //
    //   advectSubstep:305   sqrt(-beta/((1-alpha) rhoN))   sound speed
    //   advectSubstep:787   -beta/((1-alpha) rho)          c^2
    //   correctPressure:93  rhoN (1-alpha)/beta            Helmholtz
    //                                                      coefficient and
    //                                                      the Eq. 11 update
    //
    // -- the drift CANCELS out of the dynamics exactly.  The
    // overdetermined thermodynamic state is real, but PEQSI's leading
    // terms are structurally immune to it, and the drift is a
    // consistency diagnostic rather than an error that propagates.
    //
    // What is NOT covered: rho_EOS used directly (the boundary refresh),
    // and the fact that the transport properties are a consistent
    // function of (p, T, Y) while the dynamics run on the transported
    // rho -- those two states differ by the drift, and nothing here makes
    // them agree.
    if (pimple.dict().lookupOrDefault<Switch>("peqsiCoeffDrift", false))
    {
        const scalarField& rhoE = thermo_.rho()().primitiveField();
        scalar dAmax = 0, dBmax = 0, dAsum = 0, dBsum = 0, dRmax = 0;
        forAll(rho, i)
        {
            const scalar re = max(rhoE[i], small);
            const scalar ve = 1.0/re;
            // A diagnostic must not be able to kill a run: psi and xi
            // are positive for a stable fluid, but this is exactly the
            // path one enables when the state is already suspect.
            const scalar psiE = max(psi[i], vSmall);
            const scalar dpdvE = -sqr(re)/psiE;
            const scalar dpdTE =
                re*sqrt(max((cp[i] - cv[i])/max(T[i]*psiE, vSmall), 0.0));
            const scalar xiE = max(cv[i] + ve*dpdTE, small);
            const scalar dhdvE = T[i]*dpdTE + ve*dpdvE;
            const scalar aE = dpdTE/(re*xiE);
            const scalar bE = (dpdvE - dpdTE*dhdvE/xiE)/re;

            const scalar da = mag(aE - a[i])/max(mag(a[i]), small);
            const scalar db = mag(bE - b[i])/max(mag(b[i]), small);
            dAmax = max(dAmax, da); dBmax = max(dBmax, db);
            dAsum += da; dBsum += db;
            dRmax = max(dRmax, mag(re - rho[i])/max(rho[i], small));
        }
        label n = rho.size();
        reduce(n, sumOp<label>());
        reduce(dAmax, maxOp<scalar>()); reduce(dBmax, maxOp<scalar>());
        reduce(dRmax, maxOp<scalar>());
        reduce(dAsum, sumOp<scalar>()); reduce(dBsum, sumOp<scalar>());
        Info<< "PEQSI coeff drift: rho gap max = " << dRmax
            << " -> alpha " << dAmax << " max / " << dAsum/max(n, 1)
            << " mean, beta " << dBmax << " max / " << dBsum/max(n, 1)
            << " mean" << endl;
    }

    alpha_.correctBoundaryConditions();
    beta_.correctBoundaryConditions();
}


void Foam::solvers::peqsiFluid::invertTemperature()
{
    static cpuTime clk;
    static scalar tNewtLoop_ = 0, tNewtCorrect_ = 0;
    // The "loop" figure lumped the manifold closure in with the Newton
    // (the increment ran from function entry) -- it read as "the Newton
    // is slow" when most of it was the closure.  Split them.
    static scalar tClosure_ = 0, tSY_ = 0;
    clk.cpuTimeIncrement();

    // Stage 2a: composition (and the coefficient diagnostics) from the
    // manifold BEFORE the temperature inversion, so the (h, p) Newton
    // and every property behind it run on the looked-up mixture.
    if (fgmActive_)
    {
        const word zvarMode
        (
            pimple.dict().lookupOrDefault<word>("peqsiZvar", "algebraic")
        );
        if (zvarMode != "algebraic" && zvarMode != "transport")
        {
            FatalErrorInFunction
                << "peqsiZvar must be 'algebraic' or 'transport', got '"
                << zvarMode << "'.  A typo silently leaves the segregation "
                << "axis unclosed." << exit(FatalError);
        }
        if (zvarMode == "algebraic")
        {
            updateSegregation();
        }

        const bool wantSY =
            pimple.dict().lookupOrDefault<Switch>("peqsiCompSource", false);

        fgmClosure();
        tClosure_ += clk.cpuTimeIncrement();

        if (wantSY)
        {
            updateCompositionSource();
        }
        tSY_ += clk.cpuTimeIncrement();

        // Seeding the inversion from the manifold: available, OFF, and
        // measured to be a pessimisation on the cases tried.
        //
        // The reasoning that motivated it was that the temperature
        // Newton is 84% of the multicomponent step.  It is -- but not
        // because it iterates many times.  Starting from the previous
        // step's temperature, which is what the solver already does, an
        // advecting field converges in 2 iterations; the manifold guess
        // lands |T_tbl - T| ~ 1.2 K away from that and costs 5.  On the
        // gate case: T-Newton 12.23 s -> 20.60 s, 68% worse.
        //
        // The 84% is the cost of ONE iteration -- thirty-species janaf
        // and SRK departure functions over the whole field -- so it
        // yields to cheaper evaluation, not to fewer iterations.  Kept
        // switchable because a restart, or a case where the state moves
        // far per step, has no good previous temperature to start from,
        // and there the ordering may reverse.
        if
        (
            Tguess_.valid()
         && pimple.dict().lookupOrDefault<Switch>("peqsiTseed", false)
        )
        {
            volScalarField& Tw0 = const_cast<volScalarField&>(thermo_.T());
            scalarField& Tf0 = Tw0.primitiveFieldRef();
            const scalar Tmin =
                pimple.dict().lookupOrDefault<scalar>("peqsiTmin", 100);
            const scalarField& Tg = Tguess_();
            forAll(Tf0, i)
            {
                if (Tg[i] > 0)
                {
                    Tf0[i] = min(max(Tg[i], Tmin), scalar(4000));
                }
            }
            Tw0.correctBoundaryConditions();
        }
    }

    // ------------------------------------------------------------------
    // Temperature closure on the transported (h, p) state, reusing the
    // validated heRhoThermo machinery: seed he with the transported h,
    // let thermo.correct() invert T and refresh psi/mu/kappa, then
    // RESTORE the transported density (thermo.correct() would otherwise
    // overwrite it with the EOS density and destroy the scheme's
    // structural mass conservation).
    //
    // CLOSURE PAIR (peqsiConstantV, default OFF -- see the measured
    // failure at the end of this comment).
    //
    // WKK Fig. 3 inverts h(T, v) at fixed specific volume with slope xi.
    // This module used to invert at fixed p with a cp-slope Newton, on
    // the assumption that the difference was "of EOS-residual order".
    // Measured on the 1-D validation case A, it is not: the resulting
    // gap between the transported density and the EOS density ran to a
    // volume MEAN of 45%, with every cell above 1% and the worst cell
    // carrying rho = 317 against rhoEOS = 91 at the same (p, T).  The
    // fixed-p pair leaves p, T and h mutually consistent and lets rho --
    // the variable that carries mass and momentum -- float.
    //
    // The constant-v closure instead holds the TRANSPORTED rho and h and
    // solves for the pair (pEos, T) that the equation of state maps to
    // them, by a 2x2 Newton on
    //
    //     F1 = rhoEos(pEos, T) - rho_transported = 0
    //     F2 = hEos  (pEos, T) - h_transported   = 0
    //
    //     J = [ psi        drho/dT|p ]      drho/dT|p = rho^2 dpdT/dpdv
    //         [ dh/dp|T    cp        ]      dh/dp|T   = v + T dpdT/dpdv
    //
    // Converged, rhoEos == rho_transported identically, so the thermo
    // state -- and with it psi, cp, cv and therefore alpha and beta --
    // is evaluated on ONE consistent state instead of mixing the
    // transported density with properties taken at (p_transported, T).
    //
    // The inconsistency does not vanish; it MOVES to pressure, as
    // pEos - p_transported.  That is the residual PEQSI is built to
    // carry (p is the projected variable), but it is reported below so
    // the trade is measured rather than assumed.
    //
    // MEASURED, AND IT DOES NOT WORK IN THIS FORM (1-D case A):
    //   drift    0.446 -> 0.731 volume mean  (worse)
    //   |pEos-p|/p max 5.07, mean |pEos-p| 12.8 MPa against a 5 MPa field
    //   Newton pinned at 60 iterations, 998 non-convergence warnings
    // The Jacobian is right (checked against the ideal-gas limit, where
    // dh/dp|T collapses to zero as it must); what fails is linearising
    // the density residual in p.  The step it is asked to take is not
    // small: at the worst cell the transported rho is 317 where the EOS
    // gives 91 at the same T, and closing that by pressure alone needs
    // tens of MPa, far outside the range over which rho(p) is linear.
    //
    // Doing this properly needs h(T, v) evaluated directly -- which is
    // what WKK actually do -- rather than reached through h(p, T).  For
    // SRK that is analytic (ideal part plus departure at fixed v) but it
    // lives in src/SRKGas, not in the basicThermo interface, so it is a
    // deeper change than this one.  Left switchable and OFF until then;
    // the drift diagnostic below stays on either way, because a 45%
    // volume-mean gap between the mass-carrying density and the EOS is
    // worth watching whichever pair closes the state.
    //
    // The discarded EOS density is the scheme's "energy-side parking":
    // report the drift as the standing audit metric.
    // ------------------------------------------------------------------

    // The RGP-13 custom RhoFluidThermo::correct() evaluates properties AT
    // THE GIVEN T (fgmFluid contract: T comes from the manifold, not from
    // he inversion) -- seeding he and calling correct() therefore does NOT
    // invert temperature.  Do the Newton inversion here, field-wise, with
    // the pure evaluation he(p, T), then hand the converged T to the
    // thermo and let correct() refresh psi/mu/kappa at (p, T).
    {
        volScalarField& Tw = const_cast<volScalarField&>(thermo_.T());

        const scalar tol = 1e-6;    // relative, WKK Fig. 3 threshold
        const scalar dTmax = 25;    // per-pass step clamp [K]: keeps the
                                    // iterate inside the SRK validity range
                                    // across the interface cells
        // Search floor (peqsiTmin).  The default is 100 K, not the old
        // 50: at 52.5 bar the oxygen melting line sits at ~54.5 K, so a
        // 50 K result was already in the SOLID region -- SRK "solves"
        // there only because it knows nothing about the solid phase.
        // 100 K is reachable for every Z <= 0.95 mixing-line state
        // (peer-measured SRK minima: 50-95 K over that range), so a
        // cell ON this clamp now marks a genuine anomaly instead of
        // being silently parked outside physics.  Pure fuel (Z = 1)
        // has an SRK minimum of ~110 K, but its stream enters at
        // 283.4 K with a freezing point of ~245 K -- reaching 100 K
        // there is itself the pathology the clamp should flag.
        const scalar Tmin =
            pimple.dict().lookupOrDefault<scalar>("peqsiTmin", 100);
        const scalar Tmax = 4000;
        label iter = 0;
        scalar maxRel = great;
        label nSaturated_ = 0;
        label nStall_ = 0;
        scalar stallTLo_ = 0, stallTHi_ = 0;

        const Switch constantV
        (
            pimple.dict().lookupOrDefault<Switch>("peqsiConstantV", false)
        );

        scalarField& Tf = Tw.primitiveFieldRef();
        const scalarField& hf = h_.primitiveField();

        if (constantV)
        {
            // ----------------------------------------------------------
            // Constant-v closure -- WKK Fig. 3, literally: hold the
            // TRANSPORTED (rho, h) and invert temperature at fixed
            // specific volume.
            //
            // The pressure at fixed (v, T) is analytic for SRK, so the
            // 2x2 (p, T) Newton the first two attempts fought with
            // collapses to a 1-D Newton in T:
            //
            //     p(v, T) = RR T/(vs - b) - aAlpha(T)/(vs (vs + b)),
            //     vs = W v + c(T)                     [molar, Peneloux]
            //
            // and the enthalpy at that state is the LIBRARY's own
            // he(p(v,T), T) -- only the 15 lines of cubic algebra are
            // replicated, and they are cross-checked cell-by-cell
            // against the library state on first use (srkCheck below):
            // if the replica disagreed with SRKGas (root blending,
            // Peneloux c(T), unit conventions), the run refuses to
            // continue rather than silently closing on a different EOS.
            //
            // Slope: xi = (dh/dT)_v = cv + vs (dp/dT)_v with the
            // analytic (dp/dT)_v = RR/(vs-b) - aAlpha'(T)/(vs(vs+b)).
            //
            // Scope: single species (the guard below).  The FGM stage
            // replaces the replica with table-baked coefficients.
            // ----------------------------------------------------------
            if (!srkReplicaValid_)
            {
                FatalErrorInFunction
                    << "peqsiConstantV requires the single-species SRK "
                    << "replica (multi-species closes via the manifold "
                    << "coefficients in the FGM stage)" << exit(FatalError);
            }

            const volScalarField pSave("PEQSI:pSave", p_);
            const scalarField& rhoTf = rho_.primitiveField();
            scalarField& pf = p_.primitiveFieldRef();

            const scalar RR = constant::thermodynamic::RR;

            // Peneloux translation c(T) -- replica of SRKGas::cAt
            auto cAt = [this](const scalar T) -> scalar
            {
                if (srkCTlo_ <= 0) return srkC_;
                if (T <= srkCTlo_)
                    return srkCq0_ + srkCq1_*T + srkCq2_*T*T;
                if (T >= srkCThi_) return srkC_;
                const scalar cL =
                    srkCq0_ + srkCq1_*srkCTlo_ + srkCq2_*sqr(srkCTlo_);
                scalar s = (T - srkCTlo_)/(srkCThi_ - srkCTlo_);
                s = s*s*(3.0 - 2.0*s);
                return cL*(1.0 - s) + srkC_*s;
            };

            // p(v, T) and (dp/dT)_v, molar algebra; v in m^3/kg
            auto pOfVT =
                [this, &cAt, RR](const scalar v, const scalar T,
                                 scalar& dpdT) -> scalar
            {
                const scalar vs = srkW_*v + cAt(T);      // m^3/kmol
                const scalar sqrtT = sqrt(max(T, small));
                const scalar aAlpha =
                    srkCoef1_ - srkCoef2_*sqrtT + srkCoef3_*T;
                const scalar daAlpha =
                    -srkCoef2_/(2*sqrtT) + srkCoef3_;
                const scalar den1 = vs - srkB_;
                const scalar den2 = vs*(vs + srkB_);
                dpdT = RR/den1 - daAlpha/den2;
                return RR*T/den1 - aAlpha/den2;
            };

            // ---- one-time replica cross-check against the library ----
            // The incoming state is EOS-consistent by construction:
            // thermo_.rho() at the current (p, T) is what SRKGas says.
            // p(1/rhoEos, T) from the replica must reproduce p.
            if (!srkChecked_)
            {
                const scalarField& rhoEf0 =
                    thermo_.rho()().primitiveField();
                const scalarField& Tf0 = Tw.primitiveField();
                scalar worst = 0;
                forAll(Tf0, i)
                {
                    scalar dpdTdum;
                    const scalar pRep =
                        pOfVT(1.0/max(rhoEf0[i], small), Tf0[i], dpdTdum);
                    worst =
                        max(worst, mag(pRep - pf[i])/max(mag(pf[i]), small));
                }
                reduce(worst, maxOp<scalar>());
                Info<< "PEQSI closure: SRK p(v,T) replica check, "
                    << "max |p_rep - p|/p = " << worst << endl;
                if (worst > 1e-6)
                {
                    FatalErrorInFunction
                        << "SRK replica disagrees with the library EOS ("
                        << worst << "); refusing to close on a different "
                        << "equation of state" << exit(FatalError);
                }
                srkChecked_ = true;
            }

            // ---- Newton in T at fixed v, field-synchronous ----------
            // he(p, T) is a field call, so the iteration is field-wise:
            // every pass writes the trial p(v, T) into p_ (which IS the
            // thermo's p) and evaluates the library enthalpy there.
            const auto tCv = thermo_.Cv();
            const scalarField& cvf = tCv().primitiveField();

            for (; iter < 60 && maxRel > tol; ++iter)
            {
                maxRel = 0;

                forAll(Tf, i)
                {
                    scalar dpdT;
                    pf[i] =
                        max
                        (
                            pOfVT(1.0/max(rhoTf[i], small), Tf[i], dpdT),
                            1e4
                        );
                }
                p_.correctBoundaryConditions();

                const volScalarField hk(thermo_.he(p_, Tw));
                const scalarField& hkf = hk.primitiveField();

                forAll(Tf, i)
                {
                    scalar dpdT;
                    const scalar v = 1.0/max(rhoTf[i], small);
                    pOfVT(v, Tf[i], dpdT);
                    const scalar xi =
                        max(cvf[i] + (srkW_*v)*dpdT/srkW_, small);

                    scalar dT = (hf[i] - hkf[i])/xi;
                    dT = min(max(dT, -dTmax), dTmax);
                    const scalar Tnew = Tf[i] + dT;

                    if ((Tnew <= Tmin && dT < 0) || (Tnew >= Tmax && dT > 0))
                    {
                        Tf[i] = Tnew <= Tmin ? Tmin : Tmax;
                        nSaturated_++;
                        continue;
                    }
                    Tf[i] = min(max(Tnew, Tmin), Tmax);
                    maxRel = max(maxRel, mag(dT)/max(Tf[i], small));
                }

                reduce(maxRel, maxOp<scalar>());
                Tw.correctBoundaryConditions();
            }

            // ---- reachability fallback --------------------------------
            // At transcritical interface cells the transported (v, h)
            // pair lies OUTSIDE the image of the equation of state: no
            // temperature reproduces that enthalpy at that volume, and
            // p(v, T) runs into the spinodal (p <= 0, floored).
            // Measured on case A: forcing the closure there collapsed
            // pEos to the 1e4 Pa floor on 58% of the cells and the
            // garbage properties fed back through alpha/beta into the
            // dynamics.  Parking the inconsistency in p is only
            // possible where the state is reachable; elsewhere the cell
            // falls back to the fixed-p (h, p) inversion, and the count
            // is reported -- that count IS the measure of how much of
            // the field the constant-v closure can actually serve.
            label nFallback = 0;
            {
                const scalar dpFracMax = 0.5;   // |pEos-p|/p beyond which
                                                // the parking is deemed
                                                // unphysical

                DynamicList<label> bad(Tf.size()/8);
                forAll(Tf, i)
                {
                    const bool clamped =
                        (Tf[i] <= Tmin + small) || (Tf[i] >= Tmax - small);
                    const bool floored = (pf[i] <= 1e4 + small);
                    const bool wild =
                        mag(pf[i] - pSave[i])
                      > dpFracMax*max(mag(pSave[i]), small);

                    if (clamped || floored || wild)
                    {
                        bad.append(i);
                        pf[i] = pSave[i];
                    }
                }
                nFallback = bad.size();

                // fixed-p Newton on the fallback cells only
                if (returnReduce(nFallback, sumOp<label>()) > 0)
                {
                    const auto tCp = thermo_.Cp();
                    const scalarField& cpf = tCp().primitiveField();

                    for (label k = 0; k < 60; ++k)
                    {
                        scalarField Tsub(bad.size());
                        forAll(bad, m) Tsub[m] = Tf[bad[m]];
                        const tmp<scalarField> thk(thermo_.he(Tsub, bad));
                        const scalarField& hkf2 = thk();

                        scalar mrel = 0;
                        forAll(bad, m)
                        {
                            const label i = bad[m];
                            scalar dT =
                                (hf[i] - hkf2[m])/max(cpf[i], small);
                            dT = min(max(dT, -dTmax), dTmax);
                            Tf[i] = min(max(Tf[i] + dT, Tmin), Tmax);
                            mrel = max(mrel, mag(dT)/max(Tf[i], small));
                        }
                        if (mrel < tol) break;
                    }
                }

                reduce(nFallback, sumOp<label>());
            }

            // p_ now holds pEos(v, T*) on the reachable cells and the
            // transported pressure on the fallback cells; refresh every
            // property THERE, then hand the transported pressure back
            // to the solver.
            thermo_.correct();

            scalar maxDpRel = 0, sumAbs = 0;
            forAll(pf, i)
            {
                const scalar d = pf[i] - pSave[i];
                maxDpRel = max(maxDpRel, mag(d)/max(mag(pSave[i]), small));
                sumAbs += mag(d);
            }
            reduce(maxDpRel, maxOp<scalar>());
            reduce(sumAbs, sumOp<scalar>());
            label nc = pf.size();
            reduce(nc, sumOp<label>());

            p_ = pSave;
            p_.correctBoundaryConditions();

            Info<< "PEQSI closure (constant-v): |pEos - p|/p max = "
                << maxDpRel << ", mean |pEos - p| = " << sumAbs/max(nc, 1)
                << " Pa, T-Newton iters = " << iter
                << ", maxRel = " << maxRel
                << ", fallback cells = " << nFallback << endl;
        }
        else
        {
            // cp slope evaluated ONCE per step (adequate: cp varies
            // slowly over the per-step temperature increments; profiling
            // showed the per-iteration Cp() re-evaluation doubled the
            // Newton cost)
            const auto tCp = thermo_.Cp();
            const scalarField& cpf = tCp().primitiveField();

            // First pass over the whole field, then only over the cells
            // still moving.  The multicomponent enthalpy is the entire
            // cost of this closure -- 95% of the step on the 2-D case --
            // and it is paid per cell per iteration, so after the first
            // pass the field-wide he() is mostly re-deriving values that
            // have already converged.  he(T, cells) evaluates the subset
            // (the constant-v branch's fallback loop already uses it).
            DynamicList<label> active(Tf.size()/4);

            // Secant memory.  The fixed cp slope stalls exactly where
            // the closure matters most: measured 2026-08-23, 11 cells at
            // T = 165-179 K -- near-critical O2, where cp peaks through
            // pseudo-boiling -- rode the warning for 913 of 1000 steps
            // (S_Y off: zero).  cp there is both huge and fast-varying,
            // so a slope evaluated once per step points the step wrong
            // every iteration.  The secant slope costs nothing: it is
            // built from the h evaluations the loop already pays for.
            scalarField Tsec(Tf.size(), great);   // great = no history yet
            scalarField hsec(Tf.size(), 0.0);

            // Safeguarded Newton (Numerical Recipes Sec. 9.4,
            // "rtsafe": Newton-Raphson combined with bisection inside a
            // shrinking bracket of known-opposite-sign residual).  Near
            // a pseudo-boiling Widom crossing cp both peaks and varies
            // rapidly with T, so a slope sampled across one dTmax=25 K
            // span can badly misestimate the local curvature and the
            // secant overshoots back and forth without converging --
            // measured on the rd0110 3M restart (peqsiRestartH manifold
            // landed ~250k O2-branch cells in [50, 175] K, T_c = 154.6 K
            // at this pressure): flat at maxRel 0.2-0.5 for 10+ steps.
            // A first attempt (bisect on a dT SIGN flip alone, no
            // verified bracket) measured ZERO improvement -- a sign
            // flip proves the SLOPE ESTIMATE reversed, not that the
            // root lies between the last two T's, so that heuristic
            // does not actually bound the search.  This version tracks
            // the true bracket from the residual VALUES: the moment any
            // two evaluated (T, h) pairs straddle the target (which a
            // genuinely oscillating cell produces on its own), every
            // later Newton/secant step is accepted only if it lands
            // inside that bracket, otherwise the step is bisection --
            // guaranteed to converge geometrically by monotonicity
            // (cp > 0), using the SAME evaluations the loop already
            // pays for.  Cells that are converging monotonically never
            // form a bracket and are untouched by any of this.
            scalarField TloF, ThiF;
            boolList bracketed;

            // Take the temperature straight from the manifold wherever
            // the cell's specific volume still matches the table's.
            // Enthalpy is an AXIS of this table, so the tabulated T is
            // the exact inversion at the table's own specific volume --
            // the Newton on those cells is re-deriving, through 106
            // janaf and SRK departure evaluations per cell per
            // iteration, a number the table already holds.  Only the
            // volume offset needs work, and that is one linear term.
            // Cells off the manifold (useF false) keep the full Newton.
            const bool tableT =
                fgmActive_ && Tguess_.valid() && tabUsable_.valid()
             && pimple.dict().lookupOrDefault<Switch>("peqsiTableT", false);

            if (tableT)
            {
                const scalarField& Tg = Tguess_();
                const boolList& useT = tabUsable_();
                label nT = 0;
                forAll(Tf, i)
                {
                    if (useT[i] && Tg[i] > Tmin && Tg[i] < Tmax)
                    {
                        Tf[i] = Tg[i];
                        nT++;
                    }
                    else
                    {
                        active.append(i);
                    }
                }
                {
                    const label every =
                        pimple.dict().lookupOrDefault<label>
                        (
                            "peqsiDiagInterval", 10
                        );
                    static label n = 0;
                    if (every > 0 && (n++ % every) == 0)
                    {
                        reduce(nT, sumOp<label>());
                        label nA = Tf.size();
                        reduce(nA, sumOp<label>());
                        Info<< "PEQSI closure: T from manifold on " << nT
                            << "/" << nA << " cells (Newton elsewhere)"
                            << endl;
                    }
                }
                if (returnReduce(active.empty(), andOp<bool>()))
                {
                    maxRel = 0;
                }
            }

            for (; iter < 60 && maxRel > tol; ++iter)
            {
                maxRel = 0;

                if (iter == 0 && !tableT)
                {
                const volScalarField hk(thermo_.he(p_, Tw));
                const scalarField& hkf = hk.primitiveField();

                forAll(Tf, i)
                {
                    Tsec[i] = Tf[i];
                    hsec[i] = hkf[i];

                    scalar dT = (hf[i] - hkf[i])/max(cpf[i], small);
                    dT = min(max(dT, -dTmax), dTmax);
                    const scalar Tnew = Tf[i] + dT;

                    // Saturated at a clamp: the transported h lies
                    // outside the EOS window (front undershoot) -- pin
                    // and treat as converged instead of oscillating for
                    // 60 iterations (measured 934 warnings before the
                    // 2-D vortex blow-up)
                    if ((Tnew <= Tmin && dT < 0) || (Tnew >= Tmax && dT > 0))
                    {
                        Tf[i] = Tnew <= Tmin ? Tmin : Tmax;
                        nSaturated_++;
                        continue;
                    }

                    Tf[i] = min(max(Tnew, Tmin), Tmax);
                    const scalar rel = mag(dT)/max(Tf[i], small);
                    if (rel > tol) active.append(i);
                    maxRel = max(maxRel, rel);
                }
                }
                else
                {
                    if (TloF.empty())
                    {
                        TloF.setSize(Tf.size(), 0.0);
                        ThiF.setSize(Tf.size(), 0.0);
                        bracketed.setSize(Tf.size(), false);
                    }

                    scalarField Tsub(active.size());
                    forAll(active, m) Tsub[m] = Tf[active[m]];
                    const tmp<scalarField> thk(thermo_.he(Tsub, active));
                    const scalarField& hkf = thk();

                    DynamicList<label> next(active.size());
                    forAll(active, m)
                    {
                        const label i = active[m];

                        // Residual at the CURRENT Tf[i] (hkf[m] was just
                        // evaluated there).  Fold it into the bracket
                        // before deciding this iteration's step: either
                        // it completes a not-yet-bracketed cell's first
                        // sign-changing pair, or it shrinks an existing
                        // bracket -- both using an evaluation the loop
                        // needed anyway.
                        const scalar rCur = hkf[m] - hf[i];
                        if (!bracketed[i])
                        {
                            if
                            (
                                Tsec[i] < great
                             && rCur*(hsec[i] - hf[i]) < 0.0
                            )
                            {
                                if (rCur < 0.0)
                                {
                                    TloF[i] = Tf[i]; ThiF[i] = Tsec[i];
                                }
                                else
                                {
                                    TloF[i] = Tsec[i]; ThiF[i] = Tf[i];
                                }
                                bracketed[i] = true;
                            }
                        }
                        else if (rCur < 0.0)
                        {
                            TloF[i] = Tf[i];
                        }
                        else
                        {
                            ThiF[i] = Tf[i];
                        }

                        scalar slope = max(cpf[i], small);
                        if (Tsec[i] < great && mag(Tf[i] - Tsec[i]) > small)
                        {
                            const scalar sec =
                                (hkf[m] - hsec[i])/(Tf[i] - Tsec[i]);
                            // h(T) rises with T; a non-positive secant is
                            // table/EOS noise, keep the cp fallback.
                            if (sec > small) slope = sec;
                        }
                        Tsec[i] = Tf[i];
                        hsec[i] = hkf[m];

                        scalar dT = (hf[i] - hkf[m])/slope;
                        dT = min(max(dT, -dTmax), dTmax);
                        scalar Tnew = Tf[i] + dT;

                        // rtsafe safeguard: once bracketed, a step that
                        // would leave (Tlo, Thi) is untrustworthy --
                        // bisect instead.  Monotonicity (cp > 0)
                        // guarantees the bisection point is a strictly
                        // better bracket, so this cannot diverge.
                        if
                        (
                            bracketed[i]
                         && (Tnew <= TloF[i] || Tnew >= ThiF[i])
                        )
                        {
                            Tnew = 0.5*(TloF[i] + ThiF[i]);
                            dT = Tnew - Tf[i];
                        }

                        if ((Tnew <= Tmin && dT < 0) || (Tnew >= Tmax && dT > 0))
                        {
                            Tf[i] = Tnew <= Tmin ? Tmin : Tmax;
                            nSaturated_++;
                            continue;
                        }

                        Tf[i] = min(max(Tnew, Tmin), Tmax);
                        const scalar rel = mag(dT)/max(Tf[i], small);
                        if (rel > tol) next.append(i);
                        maxRel = max(maxRel, rel);
                    }
                    active.transfer(next);
                }

                reduce(maxRel, maxOp<scalar>());
            }

            // Stall census for the warning below (active goes out of
            // scope here): how many cells are still moving and in what
            // temperature band.
            nStall_ = active.size();
            stallTLo_ = great; stallTHi_ = -great;
            forAll(active, m)
            {
                stallTLo_ = min(stallTLo_, Tf[active[m]]);
                stallTHi_ = max(stallTHi_, Tf[active[m]]);
            }
        }

        tNewtLoop_ += clk.cpuTimeIncrement();

        {
            const label every =
                pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
            static label n = 0;
            if (every > 0 && (n++ % every) == 0)
            {
                label it2 = iter;
                reduce(it2, maxOp<label>());
                Info<< "PEQSI thermo closure: T-Newton iterations = " << it2
                    << ", closure " << tClosure_ << " s, S_Y " << tSY_
                    << " s, newton " << tNewtLoop_ << " s, correct() "
                    << tNewtCorrect_ << " s" << endl;
            }
        }

        // The reduction must be UNCONDITIONAL: gating a collective on a
        // per-rank count deadlocks the ranks that skip it -- exactly what
        // hung the first 48-rank rd0110 bring-up (all ranks spinning in
        // MPI progress, one saturated rank waiting in this reduce).
        reduce(nSaturated_, sumOp<label>());
        if (nSaturated_ > 0)
        {
            Info<< "PEQSI thermo closure: " << nSaturated_
                << " clamp-saturated Newton cell-iterations" << endl;
        }

        if (maxRel > tol)
        {
            // The bare maxRel cannot distinguish one pathological cell
            // from a field-wide failure, and the failing runs only ever
            // showed the number.  Say who is stalled: how many cells and
            // in what temperature band.  (S_Y off: zero warnings in
            // 1000 steps; S_Y on: a warning nearly every step -- the
            // stall is driven by the dp kick moving h between closures.)
            reduce(nStall_, sumOp<label>());
            reduce(stallTLo_, minOp<scalar>());
            reduce(stallTHi_, maxOp<scalar>());

            WarningInFunction
                << "T Newton inversion not converged: maxRel = "
                << maxRel << " after " << iter << " iterations on "
                << nStall_ << " cells, T in [" << stallTLo_ << ", "
                << stallTHi_ << "] K" << endl;
        }
        Tw.correctBoundaryConditions();
    }

    // The solver's continuity density rho_ is SEPARATE storage from the
    // thermo's EOS density: correct() refreshes the thermo state (psi,
    // mu, kappa and its own rho) at (p, T) and never touches rho_.  The
    // gap between the two densities is the family's energy-side parking
    // (cf. the transported-vs-EOS drift of the rhoTransport campaign).
    // With the constant-v closure this refresh has already happened, at
    // the CONSISTENT state (pEos, T); repeating it here would re-evaluate
    // everything at the transported pressure and throw that away.
    if (!pimple.dict().lookupOrDefault<Switch>("peqsiConstantV", false))
    {
        clk.cpuTimeIncrement();
        thermo_.correct();
        tNewtCorrect_ += clk.cpuTimeIncrement();
    }

    // Boundary values of the transported density: its 'calculated'
    // patches are updated by NO equation (the base postSolve sync that
    // used to refresh them was removed to protect the interior mass
    // ledger), so without this they stay frozen at the construction-time
    // state.  Refresh the BOUNDARY only from the EOS state so boundary
    // fluxes (inflow mass through div(phiv, rho), the Fdp convective
    // coefficient) see the current (p, T) surface state; the interior
    // stays transported.  The 1-D periodic validation could never expose
    // this; the 2-D jet's Dirichlet-state inlet does.
    {
        const volScalarField::Boundary& rhoEosB =
            thermo_.rho()().boundaryField();
        forAll(rho_.boundaryFieldRef(), patchi)
        {
            // PHYSICAL patches only: processor/cyclic patches must keep
            // their COUPLED exchange values (the neighbour cell's
            // transported density) -- blanket-assigning the one-sided
            // EOS value there broke every boundary operator at
            // processor faces (measured: parallel-vs-serial differences
            // pinned to processor-boundary cell rows, seeding the 2-D
            // shear-layer blow-up at a processor boundary)
            if (!rho_.boundaryField()[patchi].coupled())
            {
                rho_.boundaryFieldRef()[patchi] == rhoEosB[patchi];
            }
        }
        rho_.correctBoundaryConditions();
    }

    // Drift between the transported density (which carries the mass
    // ledger) and the EOS density at (p, T).  The max alone cannot say
    // whether a large value is one pathological cell or a whole region,
    // and that distinction decides whether the gap is worth closing by
    // iteration: reported here with its location and a volume-weighted
    // distribution.
    scalar maxDrift = 0;
    label maxCell = -1;
    scalar volAbove10 = 0, volAbove1 = 0, volTot = 0, volWeighted = 0;
    label nAbove10 = 0, nAbove1 = 0;
    // Interval-gated like the conservation audit, and for the same two
    // reasons: seven global reductions per step are synchronisation
    // points at high rank counts, and the per-step Pout below is exactly
    // the unbounded many-rank pattern the dpExtreme budget was added to
    // stop (both cirius MPI_ERR_TRUNCATE deaths followed Pout bursts).
    const label driftDiagN =
        pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
    static label driftCount = 0;
    const bool driftReport =
        driftDiagN > 0 && (driftCount++ % driftDiagN) == 0;
    if (driftReport)
    {
        driftHist_.setSize(64, 0);
        forAll(driftHist_, b) driftHist_[b] = 0;

        const scalarField& re = thermo_.rho()().primitiveField();
        const scalarField& rt = rho_.primitiveField();
        const scalarField& vol = mesh.V();

        // Log-binned distribution alongside max and mean.  Those two
        // bracket the population from opposite ends and neither sees
        // the part that matters: the maximum is one cell and carries
        // the boot transient, the volume mean saturates (measured
        // t^0.105, so a decade of time buys 1.26x).  Every failure this
        // campaign chased lived between them -- the cold cells grew
        // 42x in count while their mean deficit moved 31%, and the
        // reachability violations grew 56% at the worst while the count
        // moved 9%.  Sixty-four bins over 1e-6..10.
        forAll(rt, i)
        {
            const scalar d = mag(re[i] - rt[i])/max(rt[i], small);
            if (d > maxDrift)
            {
                maxDrift = d;
                maxCell = i;
            }
            volTot += vol[i];
            volWeighted += d*vol[i];
            if (d > 0.10) { volAbove10 += vol[i]; nAbove10++; }
            if (d > 0.01) { volAbove1 += vol[i]; nAbove1++; }

            const label b =
                d <= 1e-6 ? 0
              : min(label(64.0*(log10(d) + 6.0)/7.0), label(63));
            driftHist_[b]++;
        }

        scalar maxDriftLocal = maxDrift;
        reduce(maxDrift, maxOp<scalar>());
        reduce(volAbove10, sumOp<scalar>());
        reduce(volAbove1, sumOp<scalar>());
        reduce(volTot, sumOp<scalar>());
        reduce(volWeighted, sumOp<scalar>());
        reduce(nAbove10, sumOp<label>());
        reduce(nAbove1, sumOp<label>());
        Pstream::listCombineGather(driftHist_, plusEqOp<label>());
        Pstream::listCombineScatter(driftHist_);

        // Only the rank actually holding the maximum reports its state
        if (maxCell >= 0 && mag(maxDriftLocal - maxDrift) < small)
        {
            Pout<< "PEQSI drift peak: " << maxDrift
                << " at " << mesh.C()[maxCell]
                << " rho=" << rt[maxCell]
                << " rhoEOS=" << re[maxCell]
                << " p=" << p_[maxCell]
                << " T=" << thermo_.T()[maxCell]
                << " h=" << h_[maxCell] << endl;
        }
    }

    // ------------------------------------------------------------------
    // Closed-domain uniform pressure mode (peqsiClosedDomain).
    //
    // In a domain with no pressure-setting boundary -- fully periodic, or
    // all-wall -- the SPATIALLY UNIFORM part of the pressure is in the
    // null space of the Helmholtz operator, so the acoustic substep
    // cannot set it.  Measured directly on the 1-D reacting gate: the
    // mixture burns from 2659 K to 3688 K and dp stays identically 0 for
    // the whole run, while the EOS and transported densities part by 27%.
    // The gas has nowhere to put its expansion.
    //
    // The low-Mach literature has the constraint this needs.  Splitting
    // p = p0(t) + p1(x,t), the uniform p0 is fixed by a GLOBAL
    // thermodynamic constraint, not by the Poisson/Helmholtz solve:
    // constant for an open system (Cook & Riley, Phys. Fluids 8:2868,
    // 1996) and, for a closed one, by requiring the total mass to hold
    // (Nicoud, JCP 158:71, 2000; Nicoud & Poinsot 1999).  Kotturshettar,
    // Alcobia, Boldini, Pecnik & Costa (arXiv:2607.29224) write the
    // real-gas version as an inverse problem, Eq. (11),
    //     M0 = int rho(p0(t), T) dV = const,
    // solved by Newton with the isothermal compressibility
    // chi = (1/rho) drho/dp|_T, which they note may come from tabulated
    // data as well as from an analytic EOS.
    //
    // PEQSI needs that constraint with a DIFFERENT right-hand side.  In
    // the low-Mach formulation rho is evaluated from the EOS, so pinning
    // the mass pins p0.  Here rho is TRANSPORTED and its mass is already
    // exact by construction -- that is the whole point of the scheme --
    // so mass carries no information about p0.  What the uniform mode
    // must do instead is make the EOS hold IN THE MEAN:
    //
    //     p0 <- p0 + (int rho dV - int rhoEOS(p,T,Y) dV) / int rho*chi dV
    //
    // Same Newton step, same chi, target changed from the initial mass
    // to the transported mass.  It moves only the uniform mode, so the
    // mass ledger is untouched (a uniform shift adds no flux), and it is
    // the one degree of freedom that can absorb a mean EOS mismatch.
    //
    // ------------------------------------------------------------------
    // SUPERSEDED by peqsiCompSource, and NOT thermodynamically innocent.
    //
    // This shifts p uniformly to satisfy the EOS in the mean, and does
    // not touch h.  Since e = h - p/rho, that changes the internal
    // energy by -dp0/rho with nothing on the other side of the ledger.
    // Measured on the closed reacting gate, energy audit de/e:
    //
    //     S_Y only         +5.8e-09      drift 0.027   p 7.337 MPa
    //     this only        -1.57e-01     drift 0       p 7.192 MPa
    //     both             -1.64e-02     drift 1.5e-16
    //
    // It buys an exactly consistent EOS by breaking a conservation law.
    // S_Y buys the same pressure rise with energy exact to machine
    // precision, because it drives p and h together (the Lp/Lrh pair).
    //
    // This was written before S_Y existed, when the uniform mode was the
    // only way to get any pressure response at all out of a closed
    // domain, and it was the right patch then.  With the composition
    // source in place it is strictly worse.  Kept switchable for the A/B
    // above and for a non-reacting closed case, where S_Y is identically
    // zero and this is the only mechanism there is -- but never with
    // peqsiCompSource, and the warning below says so.
    //
    // Off by default; with any pressure-setting boundary the level is
    // already determined and it must not run either way.
    // ------------------------------------------------------------------
    if (pimple.dict().lookupOrDefault<Switch>("peqsiClosedDomain", false))
    {
        if (pimple.dict().lookupOrDefault<Switch>("peqsiCompSource", false))
        {
            static bool warned = false;
            if (!warned)
            {
                WarningInFunction
                    << "peqsiClosedDomain with peqsiCompSource: the "
                    << "composition source already supplies the pressure "
                    << "response, and this shifts p without a matching "
                    << "enthalpy change, so it removes energy (-1.6e-2 "
                    << "measured together, -1.6e-1 alone, against 5.8e-9 "
                    << "for the source by itself).  Turn one off." << endl;
                warned = true;
            }
        }

        const scalarField& V = mesh.V();
        const scalarField& rt = rho_.primitiveField();
        const tmp<volScalarField> trhoE(thermo_.rho());
        const scalarField& re = trhoE().primitiveField();
        const scalarField& psiF = thermo_.psi().primitiveField();

        scalar mT = 0, mE = 0, den = 0;
        forAll(rt, i)
        {
            mT += rt[i]*V[i];
            mE += re[i]*V[i];
            den += psiF[i]*V[i];        // rho*chi == psi == drho/dp|_T
        }
        reduce(mT, sumOp<scalar>());
        reduce(mE, sumOp<scalar>());
        reduce(den, sumOp<scalar>());

        if (mag(den) > vSmall)
        {
            const scalar dp0 = (mT - mE)/den;
            p_ += dimensionedScalar(p_.dimensions(), dp0);
            p_.correctBoundaryConditions();

            const label every =
                pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
            static label n = 0;
            if (every > 0 && (n++ % every) == 0)
            {
                Info<< "PEQSI closed domain: mean EOS residual "
                    << (mT - mE)/max(mT, vSmall)
                    << " -> dp0 = " << dp0 << " Pa" << endl;
            }
        }
    }

    // Internal-energy audit (peqsiEnergyAudit).
    //
    // rho e = rho h - p.  In a CLOSED adiabatic domain -- periodic, or
    // all-wall/slip -- with no work crossing the boundary, int(rho e) dV
    // must hold exactly.  Mass and rho*h are already reported; neither
    // is the conserved quantity once the pressure moves, and rho*h in
    // particular MUST grow in a constant-volume burn.
    //
    // The global figure alone hides a leak that cancels between cells,
    // so the largest single-cell change is reported with it: a
    // conservative scheme moves energy between cells and a broken one
    // creates it somewhere.  Both are normalised by the initial global
    // mean so they are directly comparable.
    if (pimple.dict().lookupOrDefault<Switch>("peqsiEnergyAudit", false))
    {
        const scalarField& V = mesh.V();
        const scalarField& rf = rho_.primitiveField();
        const scalarField& hf = h_.primitiveField();
        const scalarField& pf = p_.primitiveField();

        scalarField re(rf.size());
        forAll(re, i) re[i] = rf[i]*hf[i] - pf[i];

        scalar E = 0, Vt = 0;
        forAll(re, i) { E += re[i]*V[i]; Vt += V[i]; }
        reduce(E, sumOp<scalar>());
        reduce(Vt, sumOp<scalar>());

        if (!reInit_.valid())
        {
            reInit_.reset(new scalarField(re));
            E0_ = E;
            Info<< "PEQSI energy audit: reference int(rho e) = " << E
                << " J (mean " << E/max(Vt, vSmall) << " J/m3)" << endl;
        }
        else if (reInit_().size() == re.size())
        {
            const scalar eMean = mag(E0_)/max(Vt, vSmall);
            scalar dMax = 0;
            label iMax = -1;
            // Relative twin of the absolute metric.  |rho e| spans two
            // orders across the domain (dense fuel 1.45e9 J/m3 against
            // 8.3e6 at low Z, measured 2026-08-23 on this table), so the
            // mean-normalised maximum always lands on the fuel side no
            // matter where the defect is -- a 0.5% error there outranks
            // 20% at Z_st.  The per-cell relative maximum is the one
            // that can actually point somewhere.
            scalar rMax = 0;
            label iRel = -1;
            forAll(re, i)
            {
                const scalar d = mag(re[i] - reInit_()[i])/max(eMean, vSmall);
                if (d > dMax) { dMax = d; iMax = i; }

                const scalar r =
                    mag(re[i] - reInit_()[i])
                   /max(mag(reInit_()[i]), small*eMean);
                if (r > rMax) { rMax = r; iRel = i; }
            }
            const scalar dMaxL = dMax;
            const scalar rMaxL = rMax;
            reduce(dMax, maxOp<scalar>());
            reduce(rMax, maxOp<scalar>());

            const label every =
                pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
            static label n = 0;
            if (every > 0 && (n++ % every) == 0)
            {
                Info<< "PEQSI energy audit: global d(int rho e)/|int rho e| = "
                    << (E - E0_)/max(mag(E0_), vSmall)
                    << ", worst cell |d(rho e)|/mean = " << dMax << endl;
                if (iMax >= 0 && dMaxL == dMax)
                {
                    Pout<< "    worst at " << mesh.C()[iMax]
                        << "  rho e " << reInit_()[iMax] << " -> " << re[iMax]
                        << endl;
                }
                if (iRel >= 0 && rMaxL == rMax)
                {
                    Pout<< "    worst relative " << rMax << " at "
                        << mesh.C()[iRel]
                        << "  rho e " << reInit_()[iRel] << " -> " << re[iRel]
                        << endl;
                }
            }
        }
    }

    if (driftReport)
    {
        // Who owns the coldest cells?  The manifold cannot return below
        // its own T floor (80 K in the dp9lo2/dd2 lineage), so a cell at
        // the solver's 50 K clamp must have come out of the Newton
        // branch -- report which branch holds the minimum, and how many
        // cells sit ON the clamp, so the question is answered from the
        // log instead of by inference.
        {
            const scalar TminC =
                pimple.dict().lookupOrDefault<scalar>("peqsiTmin", 100);
            const scalarField& Tfd = thermo_.T().primitiveField();
            const boolList* useT =
                tabUsable_.valid() ? &tabUsable_() : nullptr;
            label nAtFloor = 0, nFloorTab = 0;
            scalar tMinTab = great, tMinNewt = great;
            forAll(Tfd, i)
            {
                const bool fromTab = useT ? (*useT)[i] : false;
                if (fromTab) tMinTab = min(tMinTab, Tfd[i]);
                else tMinNewt = min(tMinNewt, Tfd[i]);
                if (Tfd[i] < TminC + 0.5)
                {
                    nAtFloor++;
                    if (fromTab) nFloorTab++;
                }
            }
            // Clamp-boundary chatter census: count cells that changed
            // clamp state since the previous step.  A persistent-state
            // mismatch keeps the same cells clamped (crossings ~ 0); a
            // chatter instability flips cells every step (crossings
            // large and oscillating).
            {
                boolList cur(Tfd.size(), false);
                forAll(Tfd, i) cur[i] = (Tfd[i] < TminC + 0.5);
                label nCross = 0;
                if (clampPrev_.valid()
                 && clampPrev_().size() == cur.size())
                {
                    const boolList& prev = clampPrev_();
                    forAll(cur, i) if (cur[i] != prev[i]) nCross++;
                }
                reduce(nCross, sumOp<label>());
                Info<< "PEQSI clamp chatter: " << nCross
                    << " cells crossed the floor boundary since the "
                    << "previous census" << endl;
                clampPrev_.reset(new boolList(cur));
            }

            reduce(nAtFloor, sumOp<label>());
            reduce(nFloorTab, sumOp<label>());
            reduce(tMinTab, minOp<scalar>());
            reduce(tMinNewt, minOp<scalar>());
            Info<< "PEQSI T floor census: on clamp (<" << TminC + 0.5
                << " K) " << nAtFloor
                << " cells, of which manifold-path " << nFloorTab
                << "; Tmin manifold-path = " << tMinTab
                << ", Newton-path = " << tMinNewt << " K" << endl;

            // Attribution: a clamped cell holds a (p, h) pair the EOS
            // has no root for.  h enters from two places -- the
            // advective substep and the acoustic Eq. 24 update -- so
            // report, for the WORST clamped cell, its state and how
            // much each stage moved h this step.  hN_ is the step's
            // starting enthalpy and hStarSaved_ the value the substep
            // handed to the acoustic stage, so the split is exact.
            if (nAtFloor > 0)
            {
                scalar worstH = great;
                label worstCell = -1;
                forAll(Tfd, i)
                {
                    if (Tfd[i] < TminC + 0.5 && h_[i] < worstH)
                    {
                        worstH = h_[i];
                        worstCell = i;
                    }
                }
                scalar rep[8] = {great, 0, 0, 0, 0, 0, 0, 0};
                if (worstCell >= 0)
                {
                    const scalar hAdv =
                        hStarSaved_.valid()
                      ? hStarSaved_()[worstCell] - hN_()[worstCell]
                      : 0.0;
                    const scalar hAco =
                        hStarSaved_.valid()
                      ? h_[worstCell] - hStarSaved_()[worstCell]
                      : 0.0;
                    rep[0] = worstH;
                    rep[1] = p_[worstCell];
                    rep[2] = rho_[worstCell];
                    rep[3] = fgmActive_ ? Z_()[worstCell] : 0.0;
                    rep[4] = fgmActive_ ? Yc_()[worstCell] : 0.0;
                    rep[5] = hAdv;
                    rep[6] = hAco;
                    rep[7] = hN_()[worstCell];
                }
                // rank with the coldest h reports; ties are harmless
                scalar key = rep[0];
                reduce(key, minOp<scalar>());
                if (rep[0] == key && worstCell >= 0)
                {
                    // Where does that state sit on the table's own 4th
                    // axis?  dhF_ is the enthalpy defect the closure
                    // computed for this cell; comparing it with the
                    // axis ends says whether the cell is simply OFF the
                    // table (h too low for its own Z) rather than the
                    // victim of a bad step.
                    scalar dhCell = 0, dhLo = 0, dhHi = 0;
                    if (dhF_.valid() && fgmTable_.valid())
                    {
                        dhCell = dhF_()[worstCell];
                        const List<scalar>& ax = fgmTable_().chiAxis();
                        if (ax.size())
                        {
                            dhLo = ax[0];
                            dhHi = ax[ax.size() - 1];
                        }
                    }
                    Pout<< "PEQSI T floor worst cell: h = " << rep[0]
                        << " J/kg (h^n = " << rep[7]
                        << "), p = " << rep[1]
                        << ", rho = " << rep[2]
                        << ", Z = " << rep[3]
                        << ", Yc = " << rep[4]
                        << " | dh advective = " << rep[5]
                        << ", dh acoustic = " << rep[6]
                        << " | dh coord = " << dhCell
                        << " in axis [" << dhLo << ", " << dhHi << "]"
                        << " at " << mesh.C()[worstCell]
                        << endl;
                }
            }
        }

        // Physically-impossible-state census (table session request,
        // 2026-08-25).  Pure O2 at 52.5 bar cannot reach dh below
        // -163 kJ/kg for ANY temperature down to the triple point, and
        // isentropic expansion of the liquid streams cools by ~2 K, not
        // tens -- so a near-pure-oxidiser cell below that defect, or a
        // near-pure-fuel cell below 240 K (n-decane freezes at 243.5),
        // is a numerical undershoot by construction, not physics the
        // table should chase.
        if (fgmActive_ && dhF_.valid())
        {
            const scalarField& Zf2 = Z_().primitiveField();
            const scalarField& Tf2 = thermo_.T().primitiveField();
            const boolList* useT2 =
                tabUsable_.valid() ? &tabUsable_() : nullptr;
            label nOxBad = 0, nOxBadTab = 0, nFuelCold = 0;
            scalar worstDh = 0;
            label worstC = -1;
            forAll(Zf2, i)
            {
                if (Zf2[i] < 0.05 && dhF_()[i] < -163e3)
                {
                    nOxBad++;
                    if (useT2 && (*useT2)[i]) nOxBadTab++;
                    if (dhF_()[i] < worstDh)
                    {
                        worstDh = dhF_()[i];
                        worstC = i;
                    }
                }
                if (Zf2[i] > 0.95 && Tf2[i] < 240.0) nFuelCold++;
            }
            reduce(nOxBad, sumOp<label>());
            reduce(nOxBadTab, sumOp<label>());
            reduce(nFuelCold, sumOp<label>());
            Info<< "PEQSI impossible-state census: ox (Z<0.05, dh<-163k) "
                << nOxBad << " cells (manifold-path " << nOxBadTab
                << "), fuel (Z>0.95, T<240) " << nFuelCold
                << " cells" << endl;
            scalar key = worstDh;
            reduce(key, minOp<scalar>());
            if (worstC >= 0 && worstDh == key && key < 0)
            {
                // Both closure paths' temperatures: if the two-path
                // density mismatch hypothesis is right (manifold T vs
                // thermo-clamped T differing by 20+ K = ~10% rho on the
                // same cell, regenerated every iteration), their gap
                // widens just before an event.
                // c^2 from the ACTUAL acoustic coefficients the
                // dynamics used (App. D identity), plus the cached
                // compressibility: the A''-mechanism signature is c^2
                // sagging and psi inflating just before an event.
                const scalar c2coef =
                    -beta_[worstC]
                    /max((1.0 - alpha_[worstC])
                        *max(rho_[worstC], scalar(0.01)), vSmall);
                Pout<< "PEQSI impossible-state worst: dh = " << worstDh
                    << " J/kg at " << mesh.C()[worstC]
                    << " T(thermo) = " << Tf2[worstC]
                    << " T(manifold) = "
                    << (Tguess_.valid() ? Tguess_()[worstC] : -1.0)
                    << " p = " << p_[worstC]
                    << " c2(coef) = " << c2coef
                    << " alpha = " << alpha_[worstC]
                    << " beta = " << beta_[worstC]
                    << " rho(transported) = " << rho_[worstC]
                    << " psi = " << thermo_.psi()[worstC];
                // Covolume-wall margin: SRK dies at v = b_mix, i.e.
                // rho_wall ~ W/b (Peneloux shift cM printed alongside
                // for the exact wall).  margin > ~1 means the
                // transported rho has left the EOS domain -- the
                // sign-flip trigger the table session identified
                // (p(v > wall) < 0, the observed -GPa signature).
                if (RGfields_.size() > 4 && RGfields_[0].size())
                {
                    const scalar bMc = RGfields_[0][worstC];
                    const scalar cMc = RGfields_[4][worstC];
                    const scalar Wc = thermo_.W()()[worstC];
                    Pout<< " bM = " << bMc
                        << " cM = " << cMc
                        << " W = " << Wc
                        << " wallMargin(rho*bM/W) = "
                        << (Wc > vSmall
                              ? rho_[worstC]*bMc/Wc : -1.0);
                }
                Pout<< endl;
            }
        }

        // Compression bound alongside the temperature range.  In a
        // non-reacting run the only ways a cell can exceed the
        // isentropic temperature for the pressure it actually sits at
        // are viscous heating (bounded at ~1 K here, from the measured
        // muSgs and shear) and error.  Reporting p_max next to T_max
        // makes the comparison a subtraction instead of a separate
        // analysis, and it caught the rd0110 degradation at t = 23 us
        // where the threshold diagnostics stayed silent until 300 us.
        Info<< "PEQSI thermo closure: T = ["
            << gMin(thermo_.T().primitiveField()) << ", "
            << gMax(thermo_.T().primitiveField())
            << "] K, p = ["
            << gMin(p_.primitiveField())/1e6 << ", "
            << gMax(p_.primitiveField())/1e6
            << "] MPa, rho = ["
            << gMin(rho_.primitiveField()) << ", "
            << gMax(rho_.primitiveField())
            << "], rho drift (EOS vs transported) max = "
            << maxDrift
            << ", vol-mean = " << volWeighted/max(volTot, vSmall)
            << ", vol frac >1% = " << volAbove1/max(volTot, vSmall)
            << ", >10% = " << volAbove10/max(volTot, vSmall) << endl;

        // Quantiles of the same drift, from the log histogram.  Report
        // the upper edge of the bin the quantile falls in, as the rate
        // census does.
        {
            label tot = 0;
            forAll(driftHist_, b) tot += driftHist_[b];
            if (tot > 0)
            {
                auto q = [&](const scalar f)
                {
                    const label target = label(f*tot);
                    label acc = 0;
                    forAll(driftHist_, b)
                    {
                        acc += driftHist_[b];
                        if (acc >= target)
                        {
                            return pow(10.0, -6.0 + 7.0*(b + 1)/64.0);
                        }
                    }
                    return scalar(10);
                };
                // p99 is reported but carries no discrimination here:
                // it lands in the same bin as the median, because the
                // bulk of the field sits inside a single 1.29x bin.
                // The A/B signal is entirely in the two deeper
                // quantiles, and two of them are needed rather than
                // one -- with a single tail point there is no way to
                // separate a real divergence from bin noise.  At 3.2M
                // cells p99.9 averages 3220 cells and p99.99 averages
                // 322, both large enough to be stable.
                Info<< "PEQSI drift quantiles: p50 < " << q(0.50)
                    << ", p99 < " << q(0.99)
                    << ", p99.9 < " << q(0.999)
                    << ", p99.99 < " << q(0.9999)
                    << " (log-bin upper bounds)" << endl;
            }
        }
    }

    // Cell-COUNT fractions as well as volume fractions.  On a strongly
    // graded mesh the two say very different things: the drift lives in
    // the transcritical interface, which sits in the SMALLEST cells,
    // while the volume is dominated by the large quiescent cells far
    // downstream -- so a volume weighting flatters the number by roughly
    // the cell-size ratio.  Reporting both stops that reading.
    if (driftReport)
    {
        label nTot = rho_.primitiveField().size();
        reduce(nTot, sumOp<label>());
        Info<< "PEQSI drift by cell count: frac >1% = "
            << scalar(nAbove1)/max(nTot, 1)
            << ", >10% = " << scalar(nAbove10)/max(nTot, 1)
            << "  (" << nAbove1 << " / " << nAbove10
            << " of " << nTot << " cells)" << endl;
    }
}



void Foam::solvers::peqsiFluid::fgmClosure()
{
    // ------------------------------------------------------------------
    // Manifold lookup (stage 2a).  Per cell:
    //   c   = Yc / Cnorm(Z)
    //   dh  = h - [(1-Z) hOx + Z hFuel] - dhRef(Z, gZ, c, dh=0)
    //   Y_k, W, sourceYc, and the PEQSI coefficient diagnostics from
    //   the shared 4-axis stencil.
    //
    // What is CONSUMED at this stage: composition (written into the
    // thermo, which then serves he/psi/Cp/Cv to the closure) and the
    // PV source.  The coefficient blocks are COMPARED against the
    // runtime evaluation and reported -- switching the solver onto them
    // (peqsiUseTableCoeffs) is the follow-up step once the wiring
    // numbers are seen.  The dimensionless/raw phase rule from the
    // dp/p study is applied in the comparison so the diagnostic tests
    // the form that would actually be used.
    // ------------------------------------------------------------------
    FGMTable& tbl = fgmTable_();

    const scalarField& Zf = Z_().primitiveField();
    const scalarField& gZf = Zvar_().primitiveField();
    const scalarField& Ycf = Yc_().primitiveField();
    const scalarField& hfld = h_.primitiveField();
    const scalarField& pfld = p_.primitiveField();
    const scalarField& rhof = rho_.primitiveField();
    const scalarField& Tf = thermo_.T().primitiveField();

    scalarField& src = sourceYc_().primitiveFieldRef();

    // Re-arm on a mesh change before anything reads the cached fields.
    // Everything below is sized to the CURRENT cell count; the lookups
    // armed at the bottom of this function cache pointers into it.
    const scalar* const yData =
        thermo_.Y().size() ? thermo_.Y()[0].primitiveField().begin() : nullptr;
    const bool meshMoved =
        armedNCells_ >= 0
     && (mesh.nCells() != armedNCells_ || yData != armedYData_);

    if (meshMoved)
    {
        alphaTab_.clear();
        betaTab_.clear();
        tabUsable_.clear();
        cNormF_.clear();
        dhF_.clear();
        Tguess_.clear();
        RGfields_.clear();
        baseBlendArmed_ = false;
        tier2Armed_ = false;
        Info<< "PEQSI manifold: mesh changed (" << armedNCells_ << " -> "
            << mesh.nCells() << " cells) -- re-arming the lookups" << endl;
    }

    if (!alphaTab_.valid())
    {
        alphaTab_.reset(new scalarField(Zf.size(), 0.0));
        betaTab_.reset(new scalarField(Zf.size(), 0.0));
        tabUsable_.reset(new boolList(Zf.size(), false));
        cNormF_.reset(new scalarField(Zf.size(), 0.0));
        dhF_.reset(new scalarField(Zf.size(), 0.0));
        cnormZ_.reset(new scalarField(Zf.size(), 1.0));
        if (pimple.dict().lookupOrDefault<Switch>("peqsiCompSource", false))
        {
            rhoTabF_.reset(new scalarField(Zf.size(), 0.0));
        }

        if
        (
            tbl.hasRealGasCoeffs()
         && pimple.dict().lookupOrDefault<Switch>("peqsiTier2", false)
        )
        {
            RGfields_.setSize(tabulatedRealGasMixture::nCoeffs_);
            forAll(RGfields_, k)
            {
                RGfields_[k].setSize(Zf.size(), 0.0);
            }
        }
        Tguess_.reset(new scalarField(Zf.size(), 0.0));
    }
    scalarField& aTabF = alphaTab_();
    scalarField& bTabF = betaTab_();
    boolList& useF = tabUsable_();
    scalarField& TgF = Tguess_();

    const bool haveDTdv = tbl.hasOptTable("dTdv_h");
    const List<scalar>* dTdvT = haveDTdv ? &tbl.optTable("dTdv_h") : nullptr;

    const wordList& spn = tbl.speciesNames();
    PtrList<volScalarField>& Yall = thermo_.Y();
    const wordList& thSp = thermo_.species();
    List<scalarField*> Yp(spn.size());
    // Table pointers hoisted too: Ytable(word) is a HashTable lookup, and
    // inside the cell loop it was hashing the same 30 species names once
    // per cell -- 60k string hashes per step on the 2-D case.
    List<const List<scalar>*> Ysrc(spn.size());
    forAll(spn, k)
    {
        const label ti = findIndex(thSp, spn[k]);
        Yp[k] = &Yall[ti].primitiveFieldRef();
        Ysrc[k] = &tbl.Ytable(spn[k]);
    }

    const scalar hOx = tbl.hOx();
    const scalar hFuel = tbl.hFuel();
    const scalar pTbl = tbl.pRef();
    const scalar RR = constant::thermodynamic::RR;

    // Coefficient diagnostic.  A single max hides where the table is
    // weak, and "where" is the actionable half: the tabulation session
    // reports self-reproduction over the WHOLE manifold, this reports
    // it over the states a run actually visits, so binning by the
    // manifold coordinates turns the solver into an independent probe
    // of table quality.  Bins: c (progress) and T, which is where the
    // known weaknesses live (hot core, high c).
    scalar dAlpha = 0, dBeta = 0, dT = 0;
    label nGas = 0, nDense = 0;

    label dbCell = -1;
    FGMTable::FGMStencil stDb;
    scalar dbTab = 0, dbRun = 0, dbAlTab = 0, dbAlRun = 0;

    const label nBin = 4;
    FixedList<scalar, nBin> dTbinC(0.0), dTbinT(0.0);
    FixedList<label, nBin> nBinC(0), nBinT(0);

    const bool haveCoeffs = tbl.hasOptTable("PEQSI_alpha");
    const List<scalar>* aT =
        haveCoeffs ? &tbl.optTable("PEQSI_alpha") : nullptr;
    const List<scalar>* bnT =
        haveCoeffs ? &tbl.optTable("PEQSI_betan") : nullptr;
    const List<scalar>* bT =
        haveCoeffs ? &tbl.optTable("PEQSI_beta") : nullptr;
    const List<scalar>* xiT =
        haveCoeffs ? &tbl.optTable("PEQSI_xi") : nullptr;
    const List<scalar>* cvT =
        haveCoeffs ? &tbl.optTable("PEQSI_cv") : nullptr;
    const List<scalar>* dpdTT =
        haveCoeffs ? &tbl.optTable("PEQSI_dpdT_v") : nullptr;
    const List<scalar>* rhoTabOpt =
        tbl.hasOptTable("PEQSI_rho") ? &tbl.optTable("PEQSI_rho") : nullptr;
    const List<scalar>* dpdvnT =
        haveCoeffs && tbl.hasOptTable("PEQSI_dpdvn")
      ? &tbl.optTable("PEQSI_dpdvn") : nullptr;
    const List<scalar>* dpdvT_ =
        haveCoeffs && tbl.hasOptTable("PEQSI_dpdv_T")
      ? &tbl.optTable("PEQSI_dpdv_T") : nullptr;
    const List<scalar>* WT =
        tbl.hasOptTable("W") ? &tbl.optTable("W") : nullptr;
    const List<scalar>& Ttbl = tbl.Ttable();

    // ---- native-chi 4th axis (UFPV tables) ----
    // The 4th stencil coordinate is chi_st instead of the enthalpy
    // defect, VERBATIM from the validated fgmFluid path (fgmFluid.C
    // updateManifold): chi_tilde = 2 (D_eff/rho) |grad Z|^2 with
    // D_eff = mu/Le (+ rho nut/Sct under LES), mapped to the
    // stoichiometric value with the Pitsch-Steiner polynomial shape
    // chi_st = chi_tilde Zst(1-Zst)/(Z(1-Z)) -- the Z(1-Z) form, not
    // the erfc profile, because that is what the reference
    // implementation runs and validated.  Everything downstream --
    // dhF_ (which the S_Y perturbation loop and the Opt-1/Tier-2
    // mixture arming read as "the 4th coordinate"), the shared
    // stencil, the packed lookup -- is coordinate-agnostic.
    const bool useChi = tbl.hasChi();
    scalar shapeZst = 0, chiClampMin = 0, chiClampMax = great;
    tmp<volScalarField> tG2;
    const scalarField* g2p = nullptr;
    tmp<volScalarField> tMuLe;
    const scalarField* DeffZp = nullptr;
    if (useChi)
    {
        const scalar Zst = tbl.lookupOrDefault<scalar>("Z_st", 0.0625);
        shapeZst = max(Zst*(1.0 - Zst), small);
        chiClampMin = tbl.lookupOrDefault<scalar>("chiClampMin", 0.0);
        chiClampMax = tbl.lookupOrDefault<scalar>("chiClampMax", great);

        tG2 = magSqr(fvc::grad(Z_()));
        g2p = &tG2().primitiveField();

        // D_eff for Z: mu/Le, unity Le unless the table names one.
        // (The chi tables to date are unity-Lewis families.)
        scalar LeZ = 1.0;
        if (tbl.found("Le") && tbl.isDict("Le"))
        {
            LeZ = tbl.subDict("Le").lookupOrDefault<scalar>("Z", 1.0);
        }
        tMuLe = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject("PEQSI:DeffZ", mesh.time().name(), mesh),
                mesh,
                dimensionedScalar(dimMass/dimLength/dimTime, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        tMuLe.ref().primitiveFieldRef() = thermo_.mu()().primitiveField()/LeZ;
        if (sgsActive_ > 0)
        {
            tMuLe.ref().primitiveFieldRef() +=
                rho_.primitiveField()
               *momentumTransport->nut()().primitiveField()
               /pimple.dict().lookupOrDefault<scalar>("peqsiSct", 0.7);
        }
        DeffZp = &tMuLe().primitiveField();
    }

    // ---- packed lookup arming (once) ----
    // One interleaved gather per cell instead of ~50 scattered table
    // reads at the same stencil.  Bit-identical (see interpolatePacked);
    // costs a duplicate of the packed tables in memory, so switchable.
    const bool usePk =
        pimple.dict().lookupOrDefault<Switch>("peqsiPackedLookup", true);

    if (usePk && !packedArmed_)
    {
        DynamicList<word> pkKeys;
        auto slot = [&](const word& k) -> label
        {
            pkKeys.append(k);
            return pkKeys.size() - 1;
        };

        pkSrc_ = slot("sourcePV");
        pkT_ = slot("T");
        pkY_.setSize(spn.size());
        forAll(spn, k) pkY_[k] = slot("Y_" + spn[k]);
        pkRG_.setSize(RGfields_.size());
        forAll(RGfields_, k)
        {
            pkRG_[k] = slot("RG_" + tabulatedRealGasMixture::coeffNames()[k]);
        }
        if (rhoTabOpt) pkRho_ = slot("PEQSI_rho");
        if (WT) pkW_ = slot("W");
        if (aT) pkA_ = slot("PEQSI_alpha");
        if (bT) pkB_ = slot("PEQSI_beta");
        if (bnT) pkBn_ = slot("PEQSI_betan");
        if (xiT) pkXi_ = slot("PEQSI_xi");
        if (cvT) pkCv_ = slot("PEQSI_cv");
        if (dpdTT) pkDpdT_ = slot("PEQSI_dpdT_v");
        if (dpdvT_) pkDpdvT_ = slot("PEQSI_dpdv_T");
        if (dpdvnT) pkDpdvn_ = slot("PEQSI_dpdvn");
        if (dTdvT) pkDTdv_ = slot("dTdv_h");

        tbl.armPacked(pkKeys);
        pkVals_.setSize(pkKeys.size());
        packedArmed_ = true;
    }
    const bool pkOn = usePk && packedArmed_;


    // Discriminating test for the density-drift hypothesis
    // (peqsiFreezeY).  The claim under test is that the drift is the
    // composition's contribution to Drho/Dt, which the pressure equation
    // does not carry: alpha and beta are derived at FIXED composition,
    // so a manifold that moves Y produces a density change the flow
    // never sees.  Freezing Y after the first pass zeroes dY/dt while
    // leaving EVERYTHING else identical -- same manifold, same
    // coefficients, same temperature inversion -- so the difference
    // against the live run isolates exactly that one term.
    //
    // This is a diagnostic, not a model: with Y frozen the manifold no
    // longer closes the thermodynamics.
    static bool frozenOnce = false;
    const bool freezeY =
        pimple.dict().lookupOrDefault<Switch>("peqsiFreezeY", false)
     && frozenOnce;
    frozenOnce = true;

    const scalar vGateTol =
        pimple.dict().lookupOrDefault<scalar>("peqsiManifoldVtol", 0.02);

    const scalar sinkScale =
        pimple.dict().lookupOrDefault<scalar>("peqsiSinkScale", 1.0);

    // Density rescale of the tabulated source (peqsiSrcRhoExp, default
    // 0 = off, bit-identical).  The table is baked at one pressure, but
    // a reaction rate is a function of concentrations [X] = rho Y/W, so
    // omega ~ rho^(m-1): bimolecular m = 2 gives exponent 1, and the
    // table session measured n = 1.14 [1.10, 1.22] over the band
    // carrying 98% of the source (2026-08-25).  Scale omega by
    // (rho_cell/rho_tbl)^n.  Gated to non-negligible sources: below
    // 0.3% of the table maximum the fitted exponent diverges (12%
    // sign reversals) while the band carries 0.1% of the total -- an
    // equilibrium residual not worth amplifying with a wrong power.
    const scalar srcRhoExp =
        pimple.dict().lookupOrDefault<scalar>("peqsiSrcRhoExp", 0.0);
    scalar srcGate = 0.0;
    if (srcRhoExp != 0)
    {
        static scalar srcMaxCache = -1;
        static const void* srcMaxKey = nullptr;
        if (srcMaxKey != &tbl || srcMaxCache < 0)
        {
            const List<scalar>& sv = tbl.sourcePVTable();
            scalar m = 0;
            forAll(sv, i) m = max(m, mag(sv[i]));
            srcMaxCache = m;
            srcMaxKey = &tbl;
        }
        srcGate = 3e-3*srcMaxCache;
    }

    const scalar pvMin =
        pimple.dict().lookupOrDefault<scalar>("peqsiSourcePVMin", 0.0);
    const scalar pvMax =
        pimple.dict().lookupOrDefault<scalar>("peqsiSourcePVMax", 0.0);
    const bool pvClipOn = (pvMin < 0.0 || pvMax > 0.0);
    label nPvClip = 0;
    label nTfloor = 0;
    label nStiff1 = 0, nStiff2 = 0, nStiff3 = 0;
    scalar subMax = 0;
    const scalar dtChem = mesh.time().deltaTValue();

    const scalar Zdense =
        pimple.dict().lookupOrDefault<scalar>("peqsiZdense", 0.5);

    // Progress-variable clamp census.  The table builder was found on
    // 2026-08-23 to have wiped omega_c to exactly zero over c >~ 0.94 (a
    // window MAX that admitted the c=1 zero row as a candidate), deleting
    // the near-equilibrium SINK that pulls c back down.  Without it c has
    // nothing to stop it piling up against the c=1 clamp, which is an
    // ACCUMULATING failure -- invisible early, violent late.  That matches
    // this solver's divergence shape (step 1 dp +-0.4 Pa, step 176 +-6.5
    // kPa, then runaway) better than any instantaneous stiffness does.
    //
    // Outcome alone cannot confirm the mechanism: a new table that simply
    // runs is consistent with several stories.  The pile-up is the part
    // that is specific to this one, so measure it directly -- the count
    // above the wiped band and the count actually hitting the clamp.
    label nHiC = 0, nClampC = 0;
    label nNegPv = 0, nHiBand = 0, nHiNegPv = 0;
    // Where the highest c sits.  The restored sink's strength varies
    // with Z, so the plateau height should too -- a global max alone
    // cannot show that, and the Z-dependence is what makes the sink
    // the explanation rather than a coincidence.
    scalar zAtCmax = 0, cMaxLoZ = 0, cMaxHiZ = 0;
    // Bin by c, not by a single c >= 0.94 lump.  The restored sink lives
    // in [0.94, 1) -- at c = 1 the source is zero by construction, that
    // row being equilibrium.  A cell whose raw c has passed 1 is clamped
    // onto that row and therefore CANNOT be pulled back: it no longer
    // looks up the sink that would decelerate it.  So the lump count
    // conflates two opposite situations -- a cell sitting in the sink
    // and a cell that escaped past it -- and only the split says whether
    // the fix is reaching cells at all or merely reaching them too late.
    label nBand1 = 0, nBand1Neg = 0, nBand2 = 0, nOver = 0;
    scalar cRawMax = 0;

    forAll(Zf, celli)
    {
        const scalar Zcl = min(max(Zf[celli], 0.0), 1.0);
        const scalar gz = max(gZf[celli], 0.0);

        const scalar Cn = tbl.hasCnorm() ? tbl.interpolateCnorm(Zcl) : 1.0;
        const scalar cRaw = Ycf[celli]/max(Cn, small);
        const scalar Ccl = min(max(cRaw, 0.0), 1.0);

        if (cRaw >= 0.94) nHiC++;
        if (cRaw > 1.0) nClampC++;
        if (cRaw > cRawMax) { cRawMax = cRaw; zAtCmax = Zcl; }
        if (cRaw >= 0.94 && cRaw < 0.99) nBand1++;
        else if (cRaw >= 0.99 && cRaw <= 1.0) nBand2++;
        else if (cRaw > 1.0) nOver++;
        if (Zcl < 0.31) { if (cRaw > cMaxLoZ) cMaxLoZ = cRaw; }
        else            { if (cRaw > cMaxHiZ) cMaxHiZ = cRaw; }

        // 4th coordinate: enthalpy defect (dh tables, with the dh = 0
        // slice reference -- see FGMTable.H) or chi_st (native-chi
        // tables, mapping above).
        scalar coord4;
        if (useChi)
        {
            const scalar Dl =
                (*DeffZp)[celli]/max(rhof[celli], small);
            const scalar chiT = 2.0*Dl*(*g2p)[celli];
            const scalar shape = max(Zcl*(1.0 - Zcl), small);
            coord4 =
                max(chiClampMin,
                    min(chiClampMax, chiT*shapeZst/shape));
        }
        else
        {
            const scalar hMix = (1.0 - Zcl)*hOx + Zcl*hFuel;
            coord4 =
                hfld[celli] - hMix
              - (tbl.hasDhRef()
                  ? tbl.interpolateDhRef(Zcl, gz, Ccl) : 0.0);
        }

        cNormF_()[celli] = Ccl;
        cnormZ_()[celli] = Cn;
        dhF_()[celli] = coord4;

        FGMTable::FGMStencil st;
        tbl.makeStencil(Zcl, gz, Ccl, coord4, st);

        if (pkOn)
        {
            tbl.interpolatePacked(st, pkVals_.begin());
        }
        auto pk = [&](const label sl, const List<scalar>& t) -> scalar
        {
            return pkOn ? pkVals_[sl] : tbl.interpolate(t, st);
        };

        forAll(RGfields_, k)
        {
            RGfields_[k][celli] = pk(pkRG_[k], tbl.RGtable(k));
        }

        // Base-state table density on the stencil this loop already
        // built.  The composition-source block used to rebuild the same
        // stencil to fetch it.
        if (rhoTabF_.valid())
        {
            rhoTabF_()[celli] =
                rhoTabOpt ? pk(pkRho_, *rhoTabOpt) : 0.0;
        }

        if (!freezeY)
        {
            forAll(spn, k)
            {
                (*Yp[k])[celli] = pk(pkY_[k], *Ysrc[k]);
            }
        }
        // Diagnostic bound on the tabulated progress-variable source
        // (peqsiSourcePVMin/Max, both 0 = off).  NOT a correctness guard
        // -- binding it suppresses real physics.
        //
        // The table stores omega_c = omega_C/Cnorm(Z), because the axis
        // is c = C/Cnorm(Z).  Cnorm goes to zero with Z, so the SAME
        // physical reaction rate is amplified by 1/Cnorm, which reaches
        // 130-630 at the low-Z edge.  Verified directly on both tables --
        // undoing the normalisation lands every peak back in the physical
        // range:
        //
        //     table          Z      1/Cnorm   max|src|   max|src|*Cnorm
        //     dp2b      0.2354        1.6     9.40e+04      5.92e+04
        //     dp6       0.0379       58.6     3.98e+05      6.80e+03
        //     dp6       0.0735        7.1     2.87e+05      4.05e+04
        //     dp6       0.2903        1.7     1.05e+05      6.33e+04
        //
        // So dp6's large low-Z source is physical.  At Z = 0.0379 dp2b
        // carries 88.7 against dp6's 6803 unnormalised -- dp2b is
        // under-representing low-Z reaction by 77x, not dp6 overshooting.
        //
        // The stiffness is therefore intrinsic to the normalised progress
        // variable and cannot be tabulated away.  A CONSTANT bound is the
        // wrong shape for it: it cuts hardest exactly where 1/Cnorm is
        // largest and the physics is real.  A bound that had to exist
        // would scale with 1/Cnorm(Z).
        //
        // What the switch is good for is the count below -- how many
        // cells sit in the amplified low-Z band -- and for the A/B that
        // established the above.  Off by default.
        scalar pv = pk(pkSrc_, tbl.sourcePVTable());

        // DIAGNOSTIC ONLY -- never for production.  The table builder picks
        // the combustion branch by taking the MAX over a window, which is
        // right for a positive source but selects the WEAKEST sink once the
        // whole window is negative.  Measured 2026-08-23: the near-
        // equilibrium sink is tabulated 40-390x weaker than the flamelet
        // median over Z = 0.031-0.066.  Rebuilding the table to fix that
        // costs a bake, so scale the negative branch here first and see
        // whether a sink of the right ORDER actually stops c -- if it does
        // not, the table is exonerated and the bake would have been wasted.
        // Default 1 leaves the table untouched.
        if (pv < 0 && sinkScale != 1.0) pv *= sinkScale;

        if (srcRhoExp != 0 && rhoTabF_.valid() && mag(pv) > srcGate)
        {
            const scalar rT = rhoTabF_()[celli];
            if (rT > small)
            {
                pv *= pow(rhof[celli]/rT, srcRhoExp);
            }
        }

        if (pvClipOn)
        {
            const scalar pvRaw = pv;
            pv = min(max(pv, pvMin), pvMax);
            if (pv != pvRaw) nPvClip++;
        }
        src[celli] = rhof[celli]*pv;

        // Chemical-stiffness census.  Dc/Dt = sourcePV/Cnorm is a
        // POINTWISE ODE -- no spatial coupling -- so how many cells are
        // under-resolved, and by how much, decides whether the fix has to
        // be global (a dt limit everyone pays for) or local (sub-cycling
        // only the stiff cells).
        {
            const scalar DcDt = pv/max(Cn, small);
            const scalar frac = mag(DcDt)*dtChem/max(Ccl, 1e-3);
            if (frac > 0.1) nStiff1++;
            if (frac > 1.0) nStiff2++;
            if (frac > 10.0) nStiff3++;
            subMax = max(subMax, frac);

            // Near-equilibrium SINK census.  The pile-up at the c=1 clamp
            // is a missing sink, so counting cells that actually receive
            // one separates "c cannot come back" from "c is not being
            // pushed".  Split at c >= 0.94 because that is the band the
            // builder's window MAX wiped to exactly zero.
            if (pv < 0) nNegPv++;
            if (pv < 0 && Ccl >= 0.94 && Ccl < 0.99) nBand1Neg++;
            if (Ccl >= 0.94)
            {
                nHiBand++;
                if (pv < 0) nHiNegPv++;
            }
        }

        if (haveCoeffs && WT)
        {
            const scalar Wc = pk(pkW_, *WT);
            const scalar v = 1.0/max(rhof[celli], small);
            const scalar Zcomp =
                pfld[celli]*v*Wc/(RR*max(Tf[celli], small));

            // phase rule from the reactive dp/p study: gas cells rescale
            // the pressure-sensitive coefficients with the CELL pressure;
            // dense cells use the raw bake
            scalar betaTab;
            const bool gasCell = (Zcomp > Zdense);
            if (gasCell)
            {
                betaTab = pk(pkBn_, *bnT)*pfld[celli];
                nGas++;
            }
            else
            {
                betaTab = pk(pkB_, *bT);
                nDense++;
            }
            const scalar alphaTab = pk(pkA_, *aT);


            dAlpha =
                max(dAlpha,
                    mag(alphaTab - alpha_[celli])
                   /max(mag(alpha_[celli]), small));

            // Compare beta at MATCHED specific volume.  Raw beta is
            // steeply v-sensitive in the dense branch (measured power
            // -7.6 to -14.8 by the tabulation session), so a raw ratio
            // reports the interpolation error of v amplified by that
            // exponent, not a disagreement between table and runtime --
            // verified: the runtime assembly applied to the table's own
            // raw materials reproduces the table beta to 1.0000 at every
            // node, and the raw materials themselves agree with an
            // independent SRK evaluation to 0.6%.  beta*v is flat in v
            // (exponent -1.005 in the gas branch), so scaling by the
            // table's own v removes the amplification and leaves the
            // genuine difference.
            const scalar vTab =
                (pk(pkXi_, *xiT) - pk(pkCv_, *cvT))
               /max(pk(pkDpdT_, *dpdTT), small);
            const scalar vCell = 1.0/max(rhof[celli], small);
            const scalar rb =
                mag(betaTab*vTab - beta_[celli]*vCell)
               /max(mag(beta_[celli]*vCell), small);
            // Publish for updateCoefficients.  Usable only where the
            // cell's specific volume matches the table's: beta carries a
            // v exponent of -7.6 to -14.8 in the dense branch, so a cell
            // whose state has drifted off the manifold would receive a
            // coefficient amplified by that power.  The runtime SRK
            // evaluation serves those cells instead.
            aTabF[celli] = alphaTab;
            bTabF[celli] = betaTab;   // phase rule already applied above
            // Reachability gate.  2% in v was chosen to protect beta's
            // steep v-exponent; measured, it also lets the RAW table
            // temperature be 131 K off in the hot branch before a cell
            // is handed to the Newton, so the linear (dT/dv)_h
            // correction runs out well inside it.  Exposed because
            // tightening trades speed for accuracy directly: every cell
            // the gate rejects costs a Newton and gains the exact
            // inversion.
            // The table is baked at ONE pressure, so vTab is the specific
            // volume at pTbl.  A cell at a different pressure has a
            // genuinely different volume, and comparing the two raw was
            // measuring the pressure excursion rather than any departure
            // from the manifold.  It dominates: with the gate at 2% in v,
            //
            //     case                   |dp|/p   cells kept
            //     2-D shear, no reaction   0.9%      95.2%
            //     2-D counterflow, burning 5.6%      11.6%
            //     1-D constant-p, burning 23.5%       0.0%
            //
            // and refining dt fourfold changes none of it (35.9, 35.9,
            // 35.8%), so it is not an accumulating error -- it is where
            // the pressure went.
            //
            // Corrected with the phase rule beta already uses.  In the
            // gas branch v scales as 1/p at fixed T; in the dense branch
            // the fluid is nearly incompressible and normalising by p is
            // the wrong physics.  The tabulator measured exactly that
            // split: over +-19% in pressure the normalised form holds to
            // 0.30% in the gas phase where the raw value moves 19-43%,
            // and in the liquid fuel below ~400 K it reverses -- the
            // normalised form is 15-22% off where the raw is 2-3%.
            // v(p) by the tabulated compressibility rather than the
            // ideal-gas limit.  PEQSI_dpdvn = (dp/dv)_T v/P is
            // dimensionless and identically -1 for an ideal gas, so
            // holding it constant and integrating
            //     d(ln v)/d(ln p) = 1/dpdvn
            // gives  v = vTab (p/pTbl)^(1/dpdvn),  which collapses to
            // vTab*pTbl/p at dpdvn = -1.  Treating the exponent as
            // constant is exactly the assumption the normalisation was
            // built on -- the tabulator measured it holding to 0.30%
            // over +-19% in pressure.
            scalar vTabP = vTab;
            if (gasCell)
            {
                const scalar dpdvn =
                    dpdvnT ? pk(pkDpdvn_, *dpdvnT) : -1.0;
                // dpdvn must be negative for a stable fluid; anything
                // else is a table artefact, so fall back to ideal.
                const scalar e = (dpdvn < -small) ? 1.0/dpdvn : -1.0;
                vTabP = vTab*pow(max(pfld[celli], vSmall)/pTbl, e);
            }

            useF[celli] =
                mag(vTabP - vCell) < vGateTol*max(vCell, small);
            if (rb > dBeta)
            {
                dBeta = rb;
                stDb = st;
                dbCell = celli;
                dbTab = betaTab;
                dbRun = beta_[celli];
                dbAlTab = alphaTab;
                dbAlRun = alpha_[celli];
            }
        }

        // Newton seed.  The manifold temperature is the answer at the
        // TABLE's specific volume; the cell sits at its own, so the
        // table's (dT/dv)_h carries the guess the rest of the way.  The
        // Newton still runs and still decides -- this only changes where
        // it starts, which is worth doing because the temperature
        // inversion and the property sweeps behind it are 84% of the
        // step on the multicomponent path (2.1% for a single species):
        // thirty-species janaf and SRK departure functions, per cell,
        // per iteration.
        const scalar Ttab = pk(pkT_, Ttbl);
        if (!haveDTdv && haveCoeffs && dpdvT_)
        {
            // (dT/dv)_h = -(dh/dv)_T / (dh/dT)_v, and
            // (dh/dv)_T = T (dp/dT)_v + v (dp/dv)_T  [from du = T ds - p dv].
            // Every factor is already tabulated, so the derivative needs
            // no extra bake: a dedicated dTdv_h block would only be this
            // same product, evaluated at the same nodes.
            const scalar xiT_ = pk(pkXi_, *xiT);
            const scalar vT =
                (xiT_ - pk(pkCv_, *cvT))
               /max(pk(pkDpdT_, *dpdTT), small);
            const scalar vC = 1.0/max(rhof[celli], small);
            const scalar dhdv =
                Ttab*pk(pkDpdT_, *dpdTT)
              + vT*pk(pkDpdvT_, *dpdvT_);
            TgF[celli] = Ttab - (dhdv/max(xiT_, small))*(vC - vT);
        }
        else if (haveDTdv)
        {
            const scalar vT =
                (pk(pkXi_, *xiT) - pk(pkCv_, *cvT))
               /max(pk(pkDpdT_, *dpdTT), small);
            const scalar vC = 1.0/max(rhof[celli], small);
            TgF[celli] =
                Ttab + pk(pkDTdv_, *dTdvT)*(vC - vT);
        }
        else
        {
            TgF[celli] = Ttab;
        }

        // |T_tbl - T| is the gap between the manifold's raw temperature
        // and the one the solver is using.  It reads like an error
        // metric and is not one: where the cell is ON the manifold the
        // solver takes T from the table and the gap is ~0 by
        // construction, so what the maximum actually reports is how far
        // the RAW table value would have been wrong on the cells the
        // tabUsable_ gate rejected -- that is, the benefit of the gate.
        //
        // Measured on the 1-D reacting case, the correlation is exact:
        //     20/20 cells tabulated -> 0.3 K
        //     18-19/20             -> 131 K
        // The 131 K belongs entirely to the one or two cells whose
        // specific volume has left the manifold by more than the gate's
        // 2%, and those cells are on the Newton, not the table.  It is
        // steep -- 2% in v costing 131 K in the hot branch says the
        // linear (dT/dv)_h correction runs out well before the gate
        // does -- but the fallback catches it, and the cold 2-D case
        // sits at 4 K.
        // T_FLOOR contamination counter.  add_enthalpy_axis.py fills
        // unreachable dh nodes with 80 K on the assumption they are "a
        // corner the CFD never visits".  The tabulation session measured
        // that assumption failing: where dhRef is large an ordinary
        // mixing state looks them up, and interpolating a physical node
        // against an 80 K one puts T out by 139 K (measured at Z = 0.8).
        // The defect is in every table built before 2026-08-23.
        //
        // A run cannot tell from its own output whether it entered that
        // corner, so count it.  Anything the manifold returns below 100 K
        // is either the floor itself or an interpolation against it --
        // no state in these cases is legitimately that cold.
        if (Ttab < 100.0)
        {
            nTfloor++;
        }

        const scalar dTc = mag(Ttab - Tf[celli]);
        dT = max(dT, dTc);

        // c bins [0,.25,.5,.75,1], T bins [<600, <1200, <2000, >=2000]
        const label ic = min(label(Ccl*nBin), nBin - 1);
        dTbinC[ic] = max(dTbinC[ic], dTc); nBinC[ic]++;

        const scalar Tc2 = Tf[celli];
        const label it =
            Tc2 < 600 ? 0 : (Tc2 < 1200 ? 1 : (Tc2 < 2000 ? 2 : 3));
        dTbinT[it] = max(dTbinT[it], dTc); nBinT[it]++;
    }

    // Manifold-direction composition derivatives for S_Y: evaluate the
    // EOS at Y(c) and at Y(c + dc), same cell and same instant, so the
    // difference is a pure composition derivative with no time and
    // therefore no advection in it.  Costs two extra field evaluations
    // and two Y writes per step; affordable since Opt-1/Tier-2.
    if (pimple.dict().lookupOrDefault<Switch>("peqsiCompSource", false))
    {
        // A POSITIVE value here restores the old fixed-step difference and
        // is kept only to re-run the convergence diagnostic that condemned
        // it (see FGMTable::cInterval).  The default 0 means "difference
        // across the c node interval", which is the interpolant's exact
        // local slope and cannot depend on a step the user picked.
        const scalar dc =
            pimple.dict().lookupOrDefault<scalar>("peqsiCompSourceDc", 0.0);

        if (!dYdcRho_.valid())
        {
            dYdcRho_.reset(new scalarField(Zf.size(), 0.0));
            dYdcH_.reset(new scalarField(Zf.size(), 0.0));
        }

        // Perturb the STENCIL COORDINATE, not the Y fields.  With Opt-1
        // armed the mixture reads its composition from the stencil
        // (CfieldPtr_ == cNormF_) and never looks at Y, so a perturbation
        // written into Y is invisible to the thermo -- that mistake made
        // the derivative come out identically zero.
        const scalarField cSave(cNormF_());
        const List<scalar>* rhoTab =
            tbl.hasOptTable("PEQSI_rho") ? &tbl.optTable("PEQSI_rho") : nullptr;

        // Without the table's density there is no composition response to
        // difference, so S_Y comes out identically zero -- and a source
        // that is silently inert looks exactly like a source that is
        // correctly small.  Say so once.
        if (!rhoTab)
        {
            static bool warnedRho = false;
            if (!warnedRho)
            {
                WarningInFunction
                    << "peqsiCompSource is on but the table carries no "
                    << "PEQSI_rho block, so the composition source is "
                    << "identically zero.  Re-bake the table or turn the "
                    << "switch off." << endl;
                warnedRho = true;
            }
        }

        // The enthalpy pair is only consumed by the B term, and the
        // default drops B (peqsiCompSourceParts pressure).  Evaluating
        // he() over the whole field twice for a value nothing reads cost
        // 42% of the step.
        const word parts
        (
            pimple.dict().lookupOrDefault<word>("peqsiCompSourceParts",
                                                "pressure")
        );
        if (parts != "pressure" && parts != "both")
        {
            FatalErrorInFunction
                << "peqsiCompSourceParts must be 'pressure' or 'both', got '"
                << parts << "'.  A typo here silently selects the default "
                << "and changes the physics." << exit(FatalError);
        }
        const bool needH = (parts == "both");

        scalarField h0, h1;
        if (needH) h0 = thermo_.he(p_, thermo_.T())().primitiveField();

        forAll(Zf, celli)
        {
            if (dc > 0)
            {
                cNormF_()[celli] = min(max(cSave[celli] + dc, 0.0), 1.0);
                continue;
            }

            // Land on the far end of this cell's own c node interval.  Both
            // ends and the cell's c lie in one linear piece, so the slope
            // measured from c to that end equals the slope across the whole
            // interval -- no need to re-evaluate the low end.  When c sits
            // exactly on the upper node the step would vanish, so step down
            // instead; the interpolant is continuous, only its slope jumps.
            scalar cLo, cHi;
            tbl.cInterval(cSave[celli], cLo, cHi);
            const scalar cPert =
                (cHi - cSave[celli] > small) ? cHi : cLo;
            cNormF_()[celli] = min(max(cPert, 0.0), 1.0);
        }

        if (needH) h1 = thermo_.he(p_, thermo_.T())().primitiveField();

        forAll(Zf, celli)
        {
            const scalar step = cNormF_()[celli] - cSave[celli];
            const scalar inv = mag(step) > small ? 1.0/step : 0.0;
            dYdcH_()[celli] = needH ? (h1[celli] - h0[celli])*inv : 0.0;

            // Density response from the table.  thermo::rho() is cached
            // and would need a correct() to answer, so the manifold's own
            // density is differenced instead and scaled to the cell by
            // the ratio -- what enters S_Y is a relative response.
            scalar dR = 0;
            if (rhoTab && mag(step) > small)
            {
                const scalar Zcl = min(max(Zf[celli], 0.0), 1.0);
                const scalar gz = max(gZf[celli], 0.0);
                FGMTable::FGMStencil s1;
                tbl.makeStencil(Zcl, gz, cNormF_()[celli], dhF_()[celli], s1);
                const scalar r0 = rhoTabF_.valid() ? rhoTabF_()[celli] : 0.0;
                const scalar r1 = tbl.interpolate(*rhoTab, s1);
                if (mag(r0) > vSmall)
                {
                    dR = (r1 - r0)*inv*(rhof[celli]/r0);
                }
            }
            dYdcRho_()[celli] = dR;
        }

        cNormF_() = cSave;
    }

    forAll(spn, k)
    {
        Yall[findIndex(thSp, spn[k])].correctBoundaryConditions();
    }
    sourceYc_().correctBoundaryConditions();

    if (pvClipOn)
    {
        reduce(nPvClip, sumOp<label>());
        label nAllPv = Zf.size();
        reduce(nAllPv, sumOp<label>());
        const label every =
            pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
        static label nPv = 0;
        if (every > 0 && (nPv++ % every) == 0)
        {
            Info<< "PEQSI sourcePV clip: " << nPvClip << "/" << nAllPv
                << " cells bound to [" << pvMin << ", " << pvMax << "]"
                << endl;
        }
    }

    if (pimple.dict().lookupOrDefault<Switch>("peqsiStiffCensus", false))
    {
        reduce(nStiff1, sumOp<label>());
        reduce(nStiff2, sumOp<label>());
        reduce(nStiff3, sumOp<label>());
        reduce(subMax, maxOp<scalar>());
        label nAllS = Zf.size();
        reduce(nAllS, sumOp<label>());
        const label every =
            pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
        static label nS = 0;
        if (every > 0 && (nS++ % every) == 0)
        {
            reduce(nHiC, sumOp<label>());
            reduce(nClampC, sumOp<label>());
            reduce(cRawMax, maxOp<scalar>());
            reduce(cMaxLoZ, maxOp<scalar>());
            reduce(cMaxHiZ, maxOp<scalar>());
            reduce(nNegPv, sumOp<label>());
            reduce(nHiBand, sumOp<label>());
            reduce(nHiNegPv, sumOp<label>());
            reduce(nBand1, sumOp<label>());
            reduce(nBand1Neg, sumOp<label>());
            reduce(nBand2, sumOp<label>());
            reduce(nOver, sumOp<label>());
            Info<< "PEQSI c bins: [0.94,0.99) " << nBand1
                << " (sink on " << nBand1Neg << "), [0.99,1] " << nBand2
                << ", >1 " << nOver << endl;
            Info<< "PEQSI sink census: omega_c < 0 on " << nNegPv
                << "/" << nAllS << " cells; within c >= 0.94 band "
                << nHiNegPv << "/" << nHiBand << endl;
            Info<< "PEQSI c census: c >= 0.94 on " << nHiC
                << ", clamped (c > 1) on " << nClampC << "/" << nAllS
                << ", max raw c = " << cRawMax
                << " at Z = " << zAtCmax
                << " | max c: Z<0.31 " << cMaxLoZ
                << ", Z>=0.31 " << cMaxHiZ << endl;
            Info<< "PEQSI stiffness: |Dc/Dt| dt / c  > 0.1 on " << nStiff1
                << ", > 1 on " << nStiff2 << ", > 10 on " << nStiff3
                << " of " << nAllS << " cells; worst " << subMax << endl;
        }
    }

    reduce(nTfloor, sumOp<label>());
    if (nTfloor > 0)
    {
        label nAllT = Zf.size();
        reduce(nAllT, sumOp<label>());
        WarningInFunction
            << nTfloor << "/" << nAllT << " cells read a manifold "
            << "temperature below 100 K.  Tables built before 2026-08-23 "
            << "fill unreachable enthalpy nodes with an 80 K floor, and "
            << "interpolating against one is worth up to 139 K.  Results "
            << "from those cells are not trustworthy." << endl;
    }

    reduce(dAlpha, maxOp<scalar>());
    reduce(dBeta, maxOp<scalar>());
    reduce(dT, maxOp<scalar>());
    reduce(nGas, sumOp<label>());
    reduce(nDense, sumOp<label>());

    for (label b = 0; b < nBin; b++)
    {
        reduce(dTbinC[b], maxOp<scalar>());
        reduce(dTbinT[b], maxOp<scalar>());
        reduce(nBinC[b], sumOp<label>());
        reduce(nBinT[b], sumOp<label>());
    }

    // Worst-beta cell laid out in full: the ratio alone cannot say
    // whether the table or the runtime assembly is the odd one out, and
    // the two independently baked tables agreeing to four decimals on
    // this number says the difference is not in the table.
    // Interval-gated: the probe's raw-materials comparison evaluates
    // TWO FULL-FIELD property sweeps (Cv, Cp) for one cell's printout,
    // and the Pout is per rank per step -- measured 0.9 s/step of
    // "outer" time on the 48-rank rd0110 case before gating.
    static label dbTick = 0;
    const label dbEvery =
        pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
    const bool dbNow = dbEvery > 0 && (dbTick++ % dbEvery) == 0;

    if (dbCell >= 0 && dbNow)
    {
        const scalar v = 1.0/max(rhof[dbCell], small);
        const scalar cTab =
            sqrt(max(-dbTab/((1.0 - dbAlTab)*rhof[dbCell]), 0.0));
        const scalar cRun =
            sqrt(max(-dbRun/((1.0 - dbAlRun)*rhof[dbCell]), 0.0));
        // Raw materials on both sides.  beta is assembled from cv, psi
        // and (dp/dT)_v; if the assembled value disagrees while the
        // state matches, the disagreement lives in one of those three.
        const scalar cvRun = thermo_.Cv()().primitiveField()[dbCell];
        const scalar cpRun = thermo_.Cp()().primitiveField()[dbCell];
        const scalar psiRun = thermo_.psi().primitiveField()[dbCell];
        const scalar rhoC = rhof[dbCell];
        const scalar dpdvRun = -sqr(rhoC)/max(psiRun, small);
        const scalar dpdTRun =
            rhoC*sqrt(max((cpRun - cvRun)/(Tf[dbCell]*max(psiRun, small)), 0.0));
        const scalar xiRun = cvRun + (1.0/rhoC)*dpdTRun;

        Pout<< "PEQSI FGM raw: cv_run = " << cvRun
            << ", cv_tab = " << tbl.interpolate(tbl.optTable("PEQSI_cv"), stDb)
            << " | xi_run = " << xiRun
            << ", xi_tab = " << tbl.interpolate(tbl.optTable("PEQSI_xi"), stDb)
            << " | dpdT_run = " << dpdTRun
            << ", dpdT_tab = "
            << tbl.interpolate(tbl.optTable("PEQSI_dpdT_v"), stDb)
            << " | dpdv_run = " << dpdvRun
            << ", dpdv_tab = "
            << tbl.interpolate(tbl.optTable("PEQSI_dpdv_T"), stDb)
            << " | rho = " << rhoC
            << ", rho_tab = "
            << (tbl.interpolate(tbl.optTable("PEQSI_xi"), stDb)
              - tbl.interpolate(tbl.optTable("PEQSI_cv"), stDb))
              /max(tbl.interpolate(tbl.optTable("PEQSI_dpdT_v"), stDb), small)
            << endl;

        Pout<< "PEQSI FGM beta probe: beta_tab = " << dbTab
            << ", beta_run = " << dbRun
            << ", ratio = " << dbTab/dbRun
            << " | alpha_tab = " << dbAlTab
            << ", alpha_run = " << dbAlRun
            << " | c_tab = " << cTab << ", c_run = " << cRun
            << " | T = " << Tf[dbCell] << ", p = " << pfld[dbCell]
            << ", rho = " << rhof[dbCell] << ", v = " << v << endl;
    }

    if (dbNow)
    {
    Info<< "PEQSI FGM: coeff diag max rel dAlpha = " << dAlpha
        << ", dBeta(v-matched) = " << dBeta
        << ", |T_tbl - T| max = " << dT   // gate benefit, not error
        << " K, phase gas/dense = " << nGas << "/" << nDense << nl
        << "PEQSI FGM: |dT| by c  ";
    for (label b = 0; b < nBin; b++)
    {
        if (nBinC[b])
        {
            Info<< "[" << 0.25*b << "," << 0.25*(b + 1) << ")="
                << dTbinC[b] << "K/" << nBinC[b] << " ";
        }
    }
    Info<< nl << "PEQSI FGM: |dT| by T  ";
    {
        static const char* lab[4] =
            {"<600", "<1200", "<2000", ">=2000"};
        for (label b = 0; b < nBin; b++)
        {
            if (nBinT[b])
            {
                Info<< lab[b] << "=" << dTbinT[b] << "K/"
                    << nBinT[b] << " ";
            }
        }
    }
    Info<< endl;
    }

    // Opt-1.  thermo_.correct() rebuilds a 106-species base-thermo blend per
    // cell; under FPV the cell composition IS the table's node composition
    // interpolated on this same 16-corner stencil, so by linearity
    //   sum_corner w * nodeMix  ==  sum_species Y * thermo
    // and the blend can be read off pre-blended manifold nodes instead
    // (floating-point summation order apart).  This is the mixture-side
    // machinery fgmFluid already runs in production; peqsiFluid never armed
    // it.  Deferred to here, after the first pass has filled the stencil
    // coordinate fields, so the very first lookup sees valid coordinates.
    //
    // Armed independently of Tier-2: that path needs RG_* coefficient blocks
    // which this table does not carry, but Opt-1 needs only the internal-cell
    // recovery reference.
    if
    (
        !baseBlendArmed_
     && pimple.dict().lookupOrDefault<Switch>("peqsiBaseBlend", false)
    )
    {
        const tabulatedRealGasMixture* hook =
            dynamic_cast<const tabulatedRealGasMixture*>(&thermo_);

        if (!hook)
        {
            Info<< "PEQSI Opt-1: mixture does not implement "
                << "tabulatedRealGasMixture -> live species blend retained"
                << endl;
            baseBlendArmed_ = true;   // do not retry every step
        }
        else
        {
            // The manifold tabulates a subset of the mixture's species (30
            // of 106 here) and fgmClosure writes ONLY those; the rest are
            // never touched, so their node value is whatever they were
            // initialised to.  Zero is the physically right answer -- the
            // manifold defines the composition and a species off it is
            // absent -- but that is a claim about the case, not the table,
            // so it is checked rather than assumed.
            const wordList& sp = thermo_.species();
            const PtrList<volScalarField>& Yall = thermo_.Y();
            scalar Yoff = 0;
            label nOff = 0;
            forAll(sp, s)
            {
                if (!tbl.hasY(sp[s]))
                {
                    nOff++;
                    Yoff = max(Yoff, max(mag(Yall[s].primitiveField())));
                }
            }
            reduce(Yoff, maxOp<scalar>());

            if (Yoff > 1e-10)
            {
                Info<< "PEQSI Opt-1: " << nOff << " mixture species are not "
                    << "on the manifold and carry max |Y| = " << Yoff
                    << " -> live species blend retained" << endl;
            }
            else
            {
                if (nOff > 0)
                {
                    Info<< "PEQSI Opt-1: " << nOff << " of " << sp.size()
                        << " mixture species are off the manifold and"
                        << " identically zero (max |Y| = " << Yoff
                        << "); their node blend is exact" << endl;
                }
                List<List<scalar>> nodeY(sp.size());
                forAll(sp, s)
                {
                    nodeY[s] =
                        tbl.hasY(sp[s])
                      ? tbl.Ynodes(sp[s])
                      : List<scalar>(tbl.nTot(), scalar(0));
                }
                hook->armInternalRef(thermo_.Y()[0].primitiveField());
                hook->enableBaseBlendTabulation
                (
                    nodeY, tbl, Zf, gZf, cNormF_(), dhF_()
                );
            }

            // Tier-2.  Opt-1 removed the base-thermo species sum; the
            // O(n^2) calculateRealGas pair sum behind it is what is left,
            // and at 106 species that is 11236 pairs per cell.  The
            // mixture reads the 13 coefficients per cell instead, from
            // the RG_* blocks the tabulator bakes by running
            // calculateRealGas once per manifold node -- so this is the
            // same quantity looked up, not an approximation of it.
            if
            (
                tbl.hasRealGasCoeffs()
             && pimple.dict().lookupOrDefault<Switch>("peqsiTier2", false)
            )
            {
                const label nC = tabulatedRealGasMixture::nCoeffs_;
                List<const scalarField*> ptrs(nC);
                forAll(RGfields_, k)
                {
                    ptrs[k] = &RGfields_[k];
                }
                hook->armInternalRef(thermo_.Y()[0].primitiveField());
                hook->enableCoeffTabulation
                (
                    thermo_.Y()[0].primitiveField(), ptrs
                );
                tier2Armed_ = true;
                Info<< "PEQSI Tier-2: " << nC << " real-gas coefficients "
                    << "looked up per cell; the O(n^2) calculateRealGas "
                    << "pair sum is skipped on internal cells" << endl;
            }
            else if
            (
                !tbl.hasRealGasCoeffs()
             && pimple.dict().lookupOrDefault<Switch>("peqsiTier2", false)
            )
            {
                Info<< "PEQSI Tier-2: table carries no RG_* coefficients"
                    << " -> live pair sum retained" << endl;
            }
            baseBlendArmed_ = true;
        }
    }

    armedNCells_ = mesh.nCells();
    armedYData_ = yData;
}


Foam::tmp<Foam::volScalarField>
Foam::solvers::peqsiFluid::boundingArtDiffusivity
(
    const volScalarField& Y,
    const volScalarField& c,
    const scalar CY
) const
{
    // Terashima & Koshi, JCP 231 (2012) 6907, Eq. (27): artificial
    // diffusivity that switches on only where a [0,1] scalar leaves its
    // bounds,
    //   D* = C_Y c bar[(Y-1)H(Y-1) - Y(1-H(Y))] Delta_Y
    // with C_Y = 100 their standard value.  The bracket is the amount by
    // which the bound is violated (positive on both sides), so D* is
    // identically zero on a bounded field: resolved fronts are untouched
    // and only over/undershoots are diffused away.
    ensureDirGeometry();

    const PtrList<surfaceScalarField>& wDir = ladWDir_();
    const PtrList<volScalarField>& DeltaDir = ladDeltaDir_();

    // Bound violation, Eq. (27) bracket
    tmp<volScalarField> tviol
    (
        new volScalarField
        (
            IOobject("PEQSI:boundViol", mesh.time().name(), mesh),
            mesh,
            dimensionedScalar(dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& viol = tviol.ref();
    {
        const scalarField& yf = Y.primitiveField();
        scalarField& vf = viol.primitiveFieldRef();
        forAll(vf, i)
        {
            vf[i] = (yf[i] > 1 ? yf[i] - 1 : 0) + (yf[i] < 0 ? -yf[i] : 0);
        }
        viol.correctBoundaryConditions();
    }

    // Truncated-Gaussian overbar (same operator as the LAD terms)
    volScalarField violBar(viol);
    for (direction cmpt = 0; cmpt < 3; cmpt++)
    {
        violBar += sqr(DeltaDir[cmpt])/24.0*fvc::laplacian(wDir[cmpt], viol);
    }
    violBar = max(violBar, dimensionedScalar(dimless, 0));

    // Delta_Y = (dx dy dz)^(1/3): the cube root of the cell volume
    volScalarField DeltaY
    (
        IOobject("PEQSI:DeltaY", mesh.time().name(), mesh),
        mesh,
        dimensionedScalar(dimLength, 0),
        zeroGradientFvPatchScalarField::typeName
    );
    {
        const scalarField& V = mesh.V();
        scalarField& d = DeltaY.primitiveFieldRef();
        forAll(d, i) d[i] = cbrt(V[i]);
        DeltaY.correctBoundaryConditions();
    }

    return CY*c*violBar*DeltaY;
}


// ************************************************************************* //


void Foam::solvers::peqsiFluid::updateSegregation()
{
    // The segregation factor gZ is the table's second axis.  peqsiFluid
    // was ADVECTING it as a passive scalar with a null source -- so it
    // simply carried its initial value around and the axis never
    // responded to the flow.  That is not one of the two published
    // closures; it is a transport equation missing both of its terms.
    //
    // This is the algebraic one, ported from fgmFluid (which runs it in
    // production) so a single table and a single Cv serve both solvers:
    //
    //     Z''^2 = Cv L^2 |grad Z|^2,      gZ = Z''^2 / (Z (1 - Z))
    //
    // Cook & Riley (Phys. Fluids 6:2868, 1994) for the scale-similarity
    // form; Pierce & Moin (JFM 504:73, 2004) for L = Delta at LES
    // resolution.  Algebraic rather than a transported variance because
    // it needs no source inside the RK stages -- the fractional step
    // stays as PEQSI specifies it -- and because the transported form is
    // the rung above, warranted only once non-equilibrium mixing is
    // shown to matter.
    //
    // The advected Zvar_ is overwritten here rather than removed from
    // the RK sweep: advectSubstep.C is being edited concurrently for the
    // LAD IMEX work, and one redundant scalar advection is cheaper than
    // a merge conflict.
    if (!Zvar_.valid() || !Z_.valid())
    {
        return;
    }

    if (!varModelReported_)
    {
        const word simType(momentumTransport->lookup("simulationType"));
        if (simType == "LES")
        {
            varModel_ = varModel::les;
        }
        else if (simType == "laminar")
        {
            varModel_ = varModel::laminar;
        }
        else
        {
            // RAS needs L = k^(3/2)/epsilon and the clip fgmFluid
            // applies to it.  Not ported: every peqsiFluid case is LES
            // or laminar, and a silent laminar substitution would report
            // gZ = 0 for a case that has unresolved mixing.
            FatalErrorInFunction
                << "peqsiZvar algebraic closure supports simulationType "
                << "LES and laminar; got '" << simType << "'.  Port the "
                << "RAS mixing length from fgmFluid::varianceLengthSqr() "
                << "or set peqsiZvar transport." << exit(FatalError);
        }

        Cv_ = fgmTable_->lookupOrDefault<scalar>("Cv", 0.1);
        Cv_ = pimple.dict().lookupOrDefault<scalar>("peqsiCv", Cv_);

        Info<< "PEQSI variance: algebraic closure, simulationType "
            << simType << " -> L = "
            << (varModel_ == varModel::les ? "Delta = V^(1/3)" : "0 (laminar)")
            << ", Cv = " << Cv_ << endl;
        varModelReported_ = true;
    }

    scalarField& gZ = Zvar_().primitiveFieldRef();

    if (varModel_ == varModel::laminar)
    {
        gZ = 0.0;
        Zvar_().correctBoundaryConditions();
        return;
    }

    const volScalarField magSqrGradZ(magSqr(fvc::grad(Z_())));
    const scalarField& g2 = magSqrGradZ.primitiveField();
    const scalarField& V = mesh.V();
    const scalarField& Zf = Z_().primitiveField();

    forAll(gZ, celli)
    {
        // L^2 = (V^(1/3))^2 = V^(2/3)
        const scalar Lsqr = sqr(cbrt(V[celli]));
        const scalar Zvar = Cv_*Lsqr*g2[celli];
        const scalar Zcl = min(max(Zf[celli], 0.0), 1.0);
        gZ[celli] =
            min(max(Zvar/max(Zcl*(1.0 - Zcl), small), 0.0), 1.0);
    }

    Zvar_().correctBoundaryConditions();

    // gZ is a closure OUTPUT now, not a transported input, so its range
    // is the thing to watch: pinned at 1 means the manifold is being
    // read at maximum segregation everywhere (Cv or the filter width is
    // too large), pinned at 0 means the axis is doing nothing.
    static label segCount = 0;
    const label every =
        pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
    if (every > 0 && (segCount++ % every) == 0)
    {
        Info<< "PEQSI variance: gZ = [" << gMin(gZ) << ", " << gMax(gZ)
            << "], mean " << gAverage(gZ)
            << ", cells at the gZ = 1 clip = "
            << returnReduce
               (
                   label(std::count_if(gZ.begin(), gZ.end(),
                       [](const scalar v){ return v >= 1.0 - small; })),
                   sumOp<label>()
               )
            << endl;
    }
}


void Foam::solvers::peqsiFluid::updateCompositionSource()
{
    // S_Y, the term the pressure substep is missing under FPV.
    //
    // WKK eliminate DT/Dt between h(T,v) and p(T,v) to reach alpha and
    // beta.  Under FPV the composition is a THIRD independent variable,
    // so the same elimination leaves one more source:
    //
    //   Dp/Dt (1-alpha) = rho alpha L_h + beta div(u) + S_Y
    //   S_Y = sum_k [ (dp/dY_k)_{T,v} - rho alpha (dh/dY_k)_{T,v} ] DY_k/Dt
    //
    // ---------------------------------------------------------------
    // The rate is a MATERIAL derivative, and getting that wrong is the
    // trap.  A first version formed DY/Dt from the difference between
    // consecutive manifold lookups in the same cell.  That is the
    // EULERIAN rate dY/dt|_x, which differs from DY/Dt by u.grad(Y) --
    // and in an advecting flow that difference is the whole signal.  The
    // 1-D reacting gate did not catch it because u == 0 there, so the
    // two coincide; on the 2-D shear case the spurious advective part
    // drove the density drift from 0.388 to 8.03, a factor of 20.
    //
    // So the rate is taken from the TRANSPORT EQUATIONS instead, where
    // the material derivatives are what is actually solved for.  Under
    // FPV, Y = Y(Z, Zvar, c, dh), so
    //
    //   DY_k/Dt = (dY_k/dc) Dc/Dt + (dY_k/dZ) DZ/Dt + ...
    //
    // and for the reacting problem Dc/Dt = sourcePV/Cnorm is the
    // dominant term -- it is the only one carrying a SOURCE.  The Z and
    // gZ directions move by advection and diffusion alone: pure
    // advection contributes nothing to a material derivative, and with
    // speciesDiffusion off there is no diffusive source either, so those
    // directions are dropped.  That is also what makes the non-reacting
    // case come out at exactly zero, as it must.
    //
    // The composition derivative is then a MANIFOLD-DIRECTION difference
    // -- Y at c versus Y at c + dc, same cell, same instant -- so no
    // time and no advection enter it at all.
    //
    // Precedent: the composition contribution to the volumetric
    // constraint is standard in low-Mach reacting formulations.  Day &
    // Bell (Combust. Theory Modelling 4:535, 2000) Eq. (7) carries it as
    // the (W/W_m) omega_m term for an ideal gas -- note it is driven by
    // omega_m, the reaction rate, exactly as here -- after Majda &
    // Sethian (1985) and Najm, Wyckoff & Knio (JCP 143:381, 1998).  The
    // real-gas generalisation replaces the ideal-gas molecular-weight
    // algebra with EOS derivatives (cf. Kotturshettar, Alcobia, Boldini,
    // Pecnik & Costa, arXiv:2607.29224, who build the same constraint
    // from beta and chi but for a single-component fluid, explicitly
    // without composition terms).
    if (!sourceP_.valid())
    {
        sourceP_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "PEQSI:sourceP",
                    mesh.time().name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(p_.dimensions()/dimTime, 0)
            )
        );
    }

    scalarField& S = sourceP_().primitiveFieldRef();
    S = 0.0;

    if (!dYdcRho_.valid() || !sourceYc_.valid() || !cnormZ_.valid())
    {
        sourceP_().correctBoundaryConditions();
        return;
    }

    const scalarField& dRdc = dYdcRho_();     // (drho/dc)_{p,T} on the manifold
    const scalarField& dHdc = dYdcH_();       // (dh/dc)_{p,T}
    const scalarField& src = sourceYc_().primitiveField();   // rho * dYc/dt
    const scalarField& cn = cnormZ_();     // Cnorm(Z), not the normalised c

    const scalarField& aF = alpha_.primitiveField();
    const scalarField& rf = rho_.primitiveField();
    const scalarField& Tf = thermo_.T().primitiveField();
    const scalarField& psiF = thermo_.psi().primitiveField();
    const tmp<volScalarField> tCv(thermo_.Cv());
    const scalarField& cvF = tCv().primitiveField();

    // Only the DIRECT part of the composition response belongs here.
    //
    //   S_Y = sum_k [ (dp/dY_k)_{T,v} - rho alpha (dh/dY_k)_{T,v} ] DY_k/Dt
    //           \________ A ________/   \__________ B __________/
    //
    // The derivation puts both terms in, and both were carried at first.
    // Measurement says otherwise: A alone is stable and lands on the
    // constant-pressure endpoint, while A - rho alpha B overshoots by
    // about a factor of two and blows the run up.
    //
    // The reason is that B is already accounted for.  It is the
    // TEMPERATURE-MEDIATED half of the response -- composition changes,
    // the mixture enthalpy at fixed T changes, and T must move to keep
    // the transported h.  That is precisely what this closure's own
    // temperature inversion does: it solves h(T, p, Y) for T with the
    // manifold's Y every step, so the temperature response to a
    // composition change is captured there and alpha/beta are then
    // evaluated at the corrected state.  Adding B to the pressure
    // equation counts it a second time.
    //
    // A has no such counterpart.  It is the response at FIXED T and v,
    // and nothing else in the scheme supplies it, because p here is
    // transported rather than derived from the EOS.  That is the gap
    // the term exists to fill.
    //
    // "both" is kept for the A/B that established this.
    //
    // ---------------------------------------------------------------
    // Adversarial check on the above, because the first two arguments
    // for dropping B were weaker than they looked.
    //
    // The energy audit proves LESS than it appears to.  iL*S_Y is added
    // to Lp and to Lrh with the SAME value, so rho e = rho h - p gains
    // X - X = 0 whatever X is.  de/e ~ 6e-09 confirms the wiring is
    // consistent; it says nothing about the magnitude.  Any S_Y passes.
    //
    // The closed-box endpoint is not conservation-determined either --
    // that was the next guess and it is also wrong.  Doubling the source
    // moves it:
    //     parts=pressure   T 3866.7 K,  p 7.337 MPa
    //     parts=both       T 4000.0 K,  p 9.544 MPa
    // both conserving energy to ~1e-08.  The state is overdetermined, so
    // a given e admits different (h, p) splits and the inversion lands
    // on different T.
    //
    // What DOES pin the magnitude is the EOS, and it needs no table
    // temperature.  With e, v and Y fixed the true state is unique, so
    // the right source drives the transported state ONTO the EOS:
    //
    //     S_Y off          rho drift 0.2695
    //     parts=pressure   rho drift 0.0269      10x closer
    //     parts=both       rho drift 0.2227      overshot
    //
    // Read as a line through those three, the zero crossing sits just
    // below A, so the derived term is right to within about 10% and the
    // 2.7% residual is consistent with the rest of the discretisation
    // rather than with a missing factor.  This is independent of the
    // table lookup that the constant-pressure endpoint check relies on.
    // ---------------------------------------------------------------
    const word syParts =
        pimple.dict().lookupOrDefault<word>("peqsiCompSourceParts", "pressure");
    if (syParts != "pressure" && syParts != "both")
    {
        FatalErrorInFunction
            << "peqsiCompSourceParts must be 'pressure' or 'both', got '"
            << syParts << "'.  A typo would silently drop the enthalpy part."
            << exit(FatalError);
    }
    const bool pressureOnly = syParts != "both";

    scalar sMax = 0;
    scalar aPk = 0, bPk = 0, hPk = 0, dpk = 0, rPk = 0;
    label nSy1 = 0, nSy10 = 0, nSy100 = 0;
    scalar sFracMax = 0;
    const scalar dtCs = mesh.time().deltaTValue();
    const scalarField& pF = p_.primitiveField();
    // Cnorm cancellation probe.  c = Yc/Cnorm(Z), so dRdc = Cnorm dRho/dYc
    // and DcDt = (dYc/dt)/Cnorm: the product that enters A is Cnorm-free by
    // construction.  If the peak of |S_Y| nonetheless sits at low Z -- where
    // Cnorm is smallest -- then the cancellation is NOT happening and the
    // normalisation is amplifying the source rather than dividing out.
    // Reporting Cnorm and dRdc/Cnorm at the peak settles which it is: the
    // ratio is the physical dRho/dYc and should be O(1) across Z.
    const scalarField& Zpk = Z_().primitiveField();
    scalar zPk = 0, cnPk = 0, dcPk = 0;
    scalar rhoPk = 0, TPk = 0, psiPk = 0;
    forAll(S, i)
    {
        const scalar rho = max(rf[i], small);
        const scalar v = 1.0/rho;

        // Dc/Dt.  sourceYc_ is rho * dYc/dt and c = Yc/Cnorm(Z), so
        // Dc/Dt = sourceYc/(rho Cnorm).  Zero wherever the manifold's
        // progress-variable source is zero, which is the whole
        // non-reacting field.
        const scalar DcDt = src[i]/(rho*max(cn[i], small));

        // Composition rates at fixed (p, T), from the manifold direction
        const scalar dRho = dRdc[i]*DcDt;
        const scalar dh = dHdc[i]*DcDt;

        // (dp/dv)_T = -rho^2/psi, and
        //   sum_k (dp/dY_k)_{T,v} dY_k = -(dp/dv)_T dv|_{T,p}
        //                              =  (dp/dv)_T drho / rho^2
        const scalar dpdv = -sqr(rho)/max(psiF[i], vSmall);
        const scalar A = dpdv*dRho/sqr(rho);

        // (dh/dp)_T = v + T (dp/dT)_v/(dp/dv)_T, with xi = Cv/(1-alpha)
        // from xi = Cv + v (dp/dT)_v and (dp/dT)_v = rho alpha xi.
        const scalar oneMa = max(1.0 - aF[i], small);
        const scalar dhdp = v - Tf[i]*aF[i]*cvF[i]*psiF[i]/(rho*oneMa);

        const scalar B = dh + dhdp*A;

        // A/B discriminator (peqsiCompSourceParts).  The split shows the
        // pressure part A alone matching the physical estimate (6.6e11
        // against 5.0e11) and the enthalpy part roughly doubling it.
        // "pressure" drops the enthalpy contribution so the two can be
        // told apart: if A alone is stable, the defect is localised in
        // the (dh/dY)_{T,v} evaluation; if it is not, the problem is that
        // a source this stiff does not belong in an explicit RK stage.
        S[i] = pressureOnly ? A : (A - rho*aF[i]*B);

        // How much of the LOCAL pressure one step of this source asks
        // for.  The peak alone cannot say whether the scheme is being
        // asked something unreasonable in one cell or everywhere, and
        // that is the difference between a local pathology and a step
        // size that does not fit the physics.
        {
            const scalar frac = mag(S[i])*dtCs/max(pF[i], vSmall);
            if (frac > 0.01) nSy1++;
            if (frac > 0.10) nSy10++;
            if (frac > 1.00) nSy100++;
            sFracMax = max(sFracMax, frac);
        }

        if (mag(S[i]) > sMax)
        {
            sMax = mag(S[i]);
            aPk = A; bPk = -rho*aF[i]*B; hPk = dh; dpk = dhdp;
            rPk = dRdc[i];
            zPk = Zpk[i]; cnPk = cn[i]; dcPk = DcDt;
            rhoPk = rho; TPk = Tf[i]; psiPk = psiF[i];
        }
    }

    sourceP_().correctBoundaryConditions();

    reduce(sMax, maxOp<scalar>());
    reduce(nSy1, sumOp<label>());
    reduce(nSy10, sumOp<label>());
    reduce(nSy100, sumOp<label>());
    reduce(sFracMax, maxOp<scalar>());
    const label every =
        pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
    static label n = 0;
    if (every > 0 && (n++ % every) == 0)
    {
        Info<< "PEQSI S_Y vs p: dt|S_Y|/p > 1% on " << nSy1
            << ", > 10% on " << nSy10
            << ", > 100% on " << nSy100
            << " of " << returnReduce(S.size(), sumOp<label>())
            << " cells; worst " << sFracMax << nl
            << "PEQSI composition source: max |S_Y| = " << sMax
            << " Pa/s (dt*S_Y = " << sMax*mesh.time().deltaTValue()
            << " Pa)" << nl
            << "PEQSI S_Y split: A(pressure) = " << aPk
            << ", -rho*alpha*B(enthalpy) = " << bPk
            << " | dh/dt = " << hPk << ", dh/dp|T = " << dpk
            << ", drho/dt = " << rPk << nl
            << "PEQSI S_Y peak cell: Z = " << zPk
            << ", Cnorm = " << cnPk
            << ", Dc/Dt = " << dcPk
            << ", drho/dc = " << rPk
            << " | drho/dYc = drho/dc/Cnorm = "
            << rPk/max(cnPk, vSmall) << nl
            << "PEQSI S_Y peak EOS: rho = " << rhoPk
            << ", T = " << TPk
            << ", psi = " << psiPk
            << ", c_sound = " << sqrt(1.0/max(rhoPk*psiPk, vSmall))
            << " m/s" << endl;
    }
}
