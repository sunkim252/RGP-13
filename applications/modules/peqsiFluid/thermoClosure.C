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
#include "thermodynamicConstants.H"
#include "fvcLaplacian.H"
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

    forAll(rho, i)
    {
        const scalar v = 1.0/rho[i];
        const scalar dpdv = -sqr(rho[i])/psi[i];
        const scalar dpdT = rho[i]*sqrt(max((cp[i] - cv[i])/(T[i]*psi[i]), 0.0));
        const scalar xi = cv[i] + v*dpdT;
        const scalar dhdv = T[i]*dpdT + v*dpdv;

        a[i] = dpdT/(rho[i]*xi);
        b[i] = (dpdv - dpdT*dhdv/xi)/rho[i];

        // App. D identity: beta/(1-alpha) == -rho c^2, c^2 = gamma/psi
        // (gamma = cp/cv from the already-evaluated fields -- the
        // thermo.gamma() call re-evaluates BOTH Cp and Cv sweeps)
        const scalar lhs = b[i]/(1.0 - a[i]);
        const scalar rhs = -rho[i]*(cp[i]/cv[i])/psi[i];
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
        reduce(nTab, sumOp<label>());
        label nAll = a.size();
        reduce(nAll, sumOp<label>());

        Info<< "PEQSI coefficients: from manifold on " << nTab << "/"
            << nAll << " cells (runtime SRK elsewhere)" << endl;

        alpha_.correctBoundaryConditions();
        beta_.correctBoundaryConditions();
    }

    reduce(maxDev, maxOp<scalar>());
    Info<< "PEQSI coefficients: max |beta/(1-alpha) + rho c^2| / (rho c^2) = "
        << maxDev << endl;

    alpha_.correctBoundaryConditions();
    beta_.correctBoundaryConditions();
}


void Foam::solvers::peqsiFluid::invertTemperature()
{
    static cpuTime clk;
    static scalar tNewtLoop_ = 0, tNewtCorrect_ = 0;
    clk.cpuTimeIncrement();

    // Stage 2a: composition (and the coefficient diagnostics) from the
    // manifold BEFORE the temperature inversion, so the (h, p) Newton
    // and every property behind it run on the looked-up mixture.
    if (fgmActive_)
    {
        fgmClosure();

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
            const scalarField& Tg = Tguess_();
            forAll(Tf0, i)
            {
                if (Tg[i] > 0)
                {
                    Tf0[i] = min(max(Tg[i], scalar(50)), scalar(4000));
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
        const scalar Tmin = 50, Tmax = 4000;
        label iter = 0;
        scalar maxRel = great;
        label nSaturated_ = 0;

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
                reduce(nT, sumOp<label>());
                label nA = Tf.size();
                reduce(nA, sumOp<label>());
                Info<< "PEQSI closure: T from manifold on " << nT << "/"
                    << nA << " cells (Newton elsewhere)" << endl;
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
                    scalarField Tsub(active.size());
                    forAll(active, m) Tsub[m] = Tf[active[m]];
                    const tmp<scalarField> thk(thermo_.he(Tsub, active));
                    const scalarField& hkf = thk();

                    DynamicList<label> next(active.size());
                    forAll(active, m)
                    {
                        const label i = active[m];
                        scalar dT = (hf[i] - hkf[m])/max(cpf[i], small);
                        dT = min(max(dT, -dTmax), dTmax);
                        const scalar Tnew = Tf[i] + dT;

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
        }

        tNewtLoop_ += clk.cpuTimeIncrement();

        {
            label it2 = iter;
            reduce(it2, maxOp<label>());
            Info<< "PEQSI thermo closure: T-Newton iterations = " << it2
                << ", loop " << tNewtLoop_ << " s, correct() "
                << tNewtCorrect_ << " s" << endl;
        }

        if (nSaturated_ > 0)
        {
            reduce(nSaturated_, sumOp<label>());
            Info<< "PEQSI thermo closure: " << nSaturated_
                << " clamp-saturated Newton cell-iterations" << endl;
        }

        if (maxRel > tol)
        {
            WarningInFunction
                << "T Newton inversion not converged: maxRel = "
                << maxRel << " after " << iter << " iterations" << endl;
        }
        Tw.correctBoundaryConditions();

        if (nSaturated_ > 0)
        {
            reduce(nSaturated_, sumOp<label>());
            Info<< "PEQSI thermo closure: " << nSaturated_
                << " clamp-saturated Newton cell-iterations" << endl;
        }

        if (maxRel > tol)
        {
            WarningInFunction
                << "T Newton inversion not converged: maxRel = "
                << maxRel << " after " << iter << " iterations" << endl;
        }
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
    {
        const scalarField& re = thermo_.rho()().primitiveField();
        const scalarField& rt = rho_.primitiveField();
        const scalarField& vol = mesh.V();

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
        }

        scalar maxDriftLocal = maxDrift;
        reduce(maxDrift, maxOp<scalar>());
        reduce(volAbove10, sumOp<scalar>());
        reduce(volAbove1, sumOp<scalar>());
        reduce(volTot, sumOp<scalar>());
        reduce(volWeighted, sumOp<scalar>());
        reduce(nAbove10, sumOp<label>());
        reduce(nAbove1, sumOp<label>());

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

    Info<< "PEQSI thermo closure: T = ["
        << gMin(thermo_.T().primitiveField()) << ", "
        << gMax(thermo_.T().primitiveField())
        << "] K, rho drift (EOS vs transported) max = "
        << maxDrift
        << ", vol-mean = " << volWeighted/max(volTot, vSmall)
        << ", vol frac >1% = " << volAbove1/max(volTot, vSmall)
        << ", >10% = " << volAbove10/max(volTot, vSmall) << endl;

    // Cell-COUNT fractions as well as volume fractions.  On a strongly
    // graded mesh the two say very different things: the drift lives in
    // the transcritical interface, which sits in the SMALLEST cells,
    // while the volume is dominated by the large quiescent cells far
    // downstream -- so a volume weighting flatters the number by roughly
    // the cell-size ratio.  Reporting both stops that reading.
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

    if (!alphaTab_.valid())
    {
        alphaTab_.reset(new scalarField(Zf.size(), 0.0));
        betaTab_.reset(new scalarField(Zf.size(), 0.0));
        tabUsable_.reset(new boolList(Zf.size(), false));
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
    forAll(spn, k)
    {
        const label ti = findIndex(thSp, spn[k]);
        Yp[k] = &Yall[ti].primitiveFieldRef();
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
    const List<scalar>* dpdvT_ =
        haveCoeffs && tbl.hasOptTable("PEQSI_dpdv_T")
      ? &tbl.optTable("PEQSI_dpdv_T") : nullptr;
    const List<scalar>* WT =
        tbl.hasOptTable("W") ? &tbl.optTable("W") : nullptr;
    const List<scalar>& Ttbl = tbl.Ttable();

    const scalar Zdense =
        pimple.dict().lookupOrDefault<scalar>("peqsiZdense", 0.5);

    forAll(Zf, celli)
    {
        const scalar Zcl = min(max(Zf[celli], 0.0), 1.0);
        const scalar gz = max(gZf[celli], 0.0);

        const scalar Cn = tbl.hasCnorm() ? tbl.interpolateCnorm(Zcl) : 1.0;
        const scalar Ccl =
            min(max(Ycf[celli]/max(Cn, small), 0.0), 1.0);

        // dh with the dh = 0 slice reference (see FGMTable.H)
        const scalar hMix = (1.0 - Zcl)*hOx + Zcl*hFuel;
        const scalar dh =
            hfld[celli] - hMix - (tbl.hasDhRef() ? tbl.interpolateDhRef(Zcl, gz, Ccl) : 0.0);

        FGMTable::FGMStencil st;
        tbl.makeStencil(Zcl, gz, Ccl, dh, st);

        forAll(spn, k)
        {
            (*Yp[k])[celli] = tbl.interpolate(tbl.Ytable(spn[k]), st);
        }
        src[celli] = rhof[celli]*tbl.interpolate(tbl.sourcePVTable(), st);

        if (haveCoeffs && WT)
        {
            const scalar Wc = tbl.interpolate(*WT, st);
            const scalar v = 1.0/max(rhof[celli], small);
            const scalar Zcomp =
                pfld[celli]*v*Wc/(RR*max(Tf[celli], small));

            // phase rule from the reactive dp/p study: gas cells rescale
            // the pressure-sensitive coefficients with the CELL pressure;
            // dense cells use the raw bake
            scalar betaTab;
            if (Zcomp > Zdense)
            {
                betaTab = tbl.interpolate(*bnT, st)*pfld[celli];
                nGas++;
            }
            else
            {
                betaTab = tbl.interpolate(*bT, st);
                nDense++;
            }
            const scalar alphaTab = tbl.interpolate(*aT, st);


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
                (tbl.interpolate(*xiT, st) - tbl.interpolate(*cvT, st))
               /max(tbl.interpolate(*dpdTT, st), small);
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
            useF[celli] =
                mag(vTab - vCell) < 0.02*max(vCell, small);
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
        const scalar Ttab = tbl.interpolate(Ttbl, st);
        if (!haveDTdv && haveCoeffs && dpdvT_)
        {
            // (dT/dv)_h = -(dh/dv)_T / (dh/dT)_v, and
            // (dh/dv)_T = T (dp/dT)_v + v (dp/dv)_T  [from du = T ds - p dv].
            // Every factor is already tabulated, so the derivative needs
            // no extra bake: a dedicated dTdv_h block would only be this
            // same product, evaluated at the same nodes.
            const scalar xiT_ = tbl.interpolate(*xiT, st);
            const scalar vT =
                (xiT_ - tbl.interpolate(*cvT, st))
               /max(tbl.interpolate(*dpdTT, st), small);
            const scalar vC = 1.0/max(rhof[celli], small);
            const scalar dhdv =
                Ttab*tbl.interpolate(*dpdTT, st)
              + vT*tbl.interpolate(*dpdvT_, st);
            TgF[celli] = Ttab - (dhdv/max(xiT_, small))*(vC - vT);
        }
        else if (haveDTdv)
        {
            const scalar vT =
                (tbl.interpolate(*xiT, st) - tbl.interpolate(*cvT, st))
               /max(tbl.interpolate(*dpdTT, st), small);
            const scalar vC = 1.0/max(rhof[celli], small);
            TgF[celli] =
                Ttab + tbl.interpolate(*dTdvT, st)*(vC - vT);
        }
        else
        {
            TgF[celli] = Ttab;
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

    forAll(spn, k)
    {
        Yall[findIndex(thSp, spn[k])].correctBoundaryConditions();
    }
    sourceYc_().correctBoundaryConditions();

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
    if (dbCell >= 0)
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

    Info<< "PEQSI FGM: coeff diag max rel dAlpha = " << dAlpha
        << ", dBeta(v-matched) = " << dBeta
        << ", |T_tbl - T| max = " << dT
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
