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
#include "fvmLaplacian.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcFlux.H"
#include "fvcLaplacian.H"
#include "fvcDdt.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::peqsiFluid::pressureCorrector()
{
    // ------------------------------------------------------------------
    // Acoustic substep: modified Helmholtz equation for dp = p^{n+1}-p*
    // in the pressure-equilibrium-consistent form (PEQSI Eq. 19),
    // followed by the coupled updates (WKK Eqs. 22-24) and the
    // thermodynamic closure.  Runs once per time step (the fractional
    // step has no corrector iteration); set nCorrectors 1 in the case.
    // ------------------------------------------------------------------

    if (acousticTimeIndex_ == runTime.timeIndex())
    {
        return;
    }
    acousticTimeIndex_ = runTime.timeIndex();

    const bool timers
    (
        pimple.dict().lookupOrDefault<Switch>("peqsiTimers", false)
    );
    scalar tMark = timers ? runTime.elapsedCpuTime() : 0;
    auto mark = [&](scalar& acc)
    {
        if (timers)
        {
            const scalar now = runTime.elapsedCpuTime();
            acc += now - tMark;
            tMark = now;
        }
    };

    if (!rhoN_.valid())
    {
        FatalErrorInFunction
            << "acoustic substep called without a preceding advective "
            << "substep (momentumPredictor)" << exit(FatalError);
    }

    const dimensionedScalar& dt = runTime.deltaT();

    // Starred state (published by the advective substep)
    const volScalarField rhoStar("PEQSI:rhoStar", rho_);
    const volVectorField UStar("PEQSI:UStar", U_);
    const volScalarField pStar("PEQSI:pStar", p_);
    const volScalarField hStar("PEQSI:hStar", h_);

    // n-state coefficient combination:
    //   coef = rho^n (1-alpha)/beta   ( == -1/c^2 < 0 by WKK App. D )
    volScalarField coef
    (
        "PEQSI:coef",
        rhoN_()*(1.0 - alpha_)/beta_
    );
    // coef is -1/c^2 < 0 for any physical state (App. D identity); a
    // cell whose TRANSPORTED rho^n has drifted negative flips its sign
    // and the Helmholtz diagonal with it -- the operator goes indefinite
    // and one solve detonates (rd0110 lip episode: dp reached 1e18 Pa in
    // a single step).  Ceiling just below zero decouples such a cell
    // from the acoustic update (its dp interpolates from the neighbours,
    // the Eq. 11-13 increments vanish there) until advection restores
    // the transported rho.
    coef.min(dimensionedScalar(coef.dimensions(), -vSmall));

    // dp boundary types derived from the p field: a fixed-pressure
    // boundary (the 2-D jet outlet: p fixed to chamber pressure, PEQSI
    // Sec. III B) must hold dp = 0 there; everything else (walls, inlets,
    // constraint patches) is zero-gradient.
    wordList dpBcTypes
    (
        p_.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );
    forAll(p_.boundaryField(), patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(p_.boundaryField()[patchi]))
        {
            dpBcTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
        // A waveTransmissive outlet declares a far-field PRESSURE the
        // boundary relaxes toward -- for the acoustic correction that is
        // a pressure-fixed far field, which is exactly the reference's
        // own outlet treatment (PEQSI Sec. III B: outlet p fixed to the
        // chamber value).  The old zeroGradient mapping made it a hard
        // acoustic REFLECTOR instead: on the rd0110 restart the ringing
        // could never leave, and the non-telescoping h term rectified
        // the sustained oscillation into a systematic heating (88% of
        // cells +100 K in 300 ns, global rho e drifting monotonically)
        // until a boundary EOS evaluation died at step 352.
        else if (p_.boundaryField()[patchi].type() == "waveTransmissive")
        {
            dpBcTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField dp
    (
        IOobject("PEQSI:dp", runTime.name(), mesh),
        mesh,
        dimensionedScalar(dimPressure, 0),
        dpBcTypes
    );

    // Implicit convective-coefficient face flux:
    //   d/dx_i ( 2 rho^n u^n_i (1-alpha)/beta dp / dt )
    const surfaceScalarField Fdp
    (
        "PEQSI:Fdp",
        fvc::flux(2.0*coef*UN_())/dt
    );

    // Helmholtz RHS -- the PEQSI Eq. (18) consistency substitution
    // evaluated with the SUBSTEP'S OWN quadrature (default), or the
    // paper's trapezoidal form (A/B switch):
    //
    // default (substep-consistent): source = (4/dt) sComp, where sComp is
    //   the advective substep's accumulated compression bookkeeping
    //   (1/6 rho^n + 1/6 rho^(1) + 2/3 rho^(2)) div(phiv).  Integrating
    //   the Helmholtz equation over a periodic domain, the laplacian and
    //   convective terms telescope, leaving int(coef dp) = dt int(sComp)
    //   -- exactly the substep's mass change (its conservative-flux part
    //   telescopes too), so Eq. (22) restores global mass conservation
    //   identically, for any interpolation scheme.  Locally sComp is of
    //   the rho div(u) form (zero at a passive contact), so the
    //   pressure-equilibrium property is preserved.  This is Eq. (18)'s
    //   own principle -- "substitute what the substep actually solved" --
    //   applied to our SSP-RK3 discretisation instead of the reference's
    //   trapezoid.
    //
    // trapezoid form (PEQSI Eq. 19 literal): (2/dt)(rho^n div u^n +
    //   rho^* div u^*)/... exact only if the substep satisfies the
    //   trapezoidal advective continuity; with our substep the assumption
    //   error was measured at +5.6% mass per period (case A).  Kept for
    //   the A/B against the reference.
    //
    // (The un-substituted WKK Eq. 28 residual form was also tried: it is
    // mass-exact by the same telescoping argument, but its RHS carries
    // the interface-localised divergence mismatch and blew up within 6
    // steps -- the very pathology PEQSI Sec. II B documents.)
    const Switch consistencyRHS
    (
        pimple.dict().lookupOrDefault<Switch>("peqsiTrapezoidRHS", false)
    );

    // Temporal order of the Helmholtz equation (WKK App. E, Eq. E.3):
    // theta = 0.5 -> second order, theta = 1.0 -> first order.  The
    // references use FIRST order everywhere except the smoothed 1-D
    // case C: 1-D cases A/B "diverged with the second order" (PEQSI
    // Sec. III A 2), and for 2-D/3-D "a converged solution in the
    // pressure correction step cannot be obtained" at second order
    // (PEQSI Sec. IV).  Default therefore 1.0; set 0.5 only for the
    // case-C temporal-accuracy study.  The mass-conservation telescoping
    // is theta-independent: int(coef dp) = dt int(sComp) either way.
    const scalar theta
    (
        pimple.dict().lookupOrDefault<scalar>("peqsiHelmholtzTheta", 1.0)
    );

    tmp<volScalarField> tRhs;
    if (consistencyRHS)
    {
        // Literal theta-weighted trapezoid (E.3 with the PEQSI Eq. 18/19
        // consistency substitution applied at each end point)
        tRhs =
            (2.0/(theta*dt))
           *(
                (1.0 - theta)*rhoN_()*fvc::div(UN_())
              + theta*rhoStar*fvc::div(UStar)
            );
    }
    else
    {
        tRhs = (2.0/(theta*dt))*sComp_();
    }

    // Non-orthogonal correction loop, the standard pressure-equation
    // pattern.  With the 'corrected' laplacian scheme the implicit
    // operator carries only the ORTHOGONAL part; the cross-diffusion
    // correction is explicit, evaluated with the current dp -- which on
    // a single solve is the INITIAL dp = 0, so on a non-orthogonal mesh
    // the correction was simply missing while the RHS's explicit
    // laplacian of (p* + p^n) kept its full corrected form.  On an
    // orthogonal mesh the correction is identically zero and one pass
    // is bit-identical to the old single solve; on an unstructured mesh
    // set nNonOrthogonalCorrectors as for any OpenFOAM p equation.
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix dpEqn
        (
            fvm::laplacian(dp)
          + fvm::div(Fdp, dp)
          + fvm::Sp((2.0/theta)*coef/sqr(dt), dp)
         ==
          - fvc::laplacian(pStar + pN_())
          + tRhs()
        );

        mark(tPhase_[4]);

        dpEqn.solve();

        mark(tPhase_[5]);
    }

    // Interval-gated: two global reductions and a log line per step add
    // nothing the interval's last report does not, and at high rank
    // counts the reductions are synchronisation points.
    {
        const label diagN =
            pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
        static label nDp = 0;
        if (diagN > 0 && (nDp++ % diagN) == 0)
        {
            Info<< "PEQSI: dp min/max = "
                << gMin(dp.primitiveField()) << " / "
                << gMax(dp.primitiveField()) << " Pa (theta = " << theta
                << ")" << endl;
        }
    }

    // S_Y chain instrumentation (peqsiTraceSY).  The composition source
    // delivers a pressure impulse of the right size -- summed over the
    // burn it matches what the closed case actually gains -- yet almost
    // no velocity comes out.  This prints every link so the loss can be
    // located instead of guessed at:
    //
    //   S_Y            the source itself                       [Pa/s]
    //   dt*iL*S_Y      what it adds to p* in one step           [Pa]
    //   spread(p*)     how much of that is NON-uniform, i.e. the only
    //                  part a laplacian can see                 [Pa]
    //   dp             what the Helmholtz actually returned     [Pa]
    //   max|U|         what reached the velocity field          [m/s]
    //
    // A uniform p* increment has zero gradient, so if spread(p*) is tiny
    // next to dt*iL*S_Y the answer is geometric, not a coupling defect:
    // in a uniform burn the expansion can only be communicated from a
    // pressure-setting boundary, at the speed of sound.
    if (pimple.dict().lookupOrDefault<Switch>("peqsiTraceSY", false))
    {
        const scalar dtv = mesh.time().deltaTValue();
        scalar sMax = 0, kick = 0;
        if (sourceP_.valid())
        {
            const scalarField& SY = sourceP_().primitiveField();
            const scalarField& aF = alpha_.primitiveField();
            forAll(SY, i)
            {
                sMax = max(sMax, mag(SY[i]));
                kick = max(kick, mag(dtv*SY[i]/max(1.0 - aF[i], small)));
            }
            reduce(sMax, maxOp<scalar>());
            reduce(kick, maxOp<scalar>());
        }
        const scalar pLo = gMin(p_.primitiveField());
        const scalar pHi = gMax(p_.primitiveField());
        Info<< "PEQSI S_Y chain: S_Y = " << sMax << " Pa/s"
            << ", dt*iL*S_Y = " << kick << " Pa"
            << ", spread(p) = " << (pHi - pLo) << " Pa"
            << ", dp = " << (gMax(dp.primitiveField())
                           - gMin(dp.primitiveField())) << " Pa"
            << ", max|U| = " << gMax(mag(U_.primitiveField()))
            << " m/s" << endl;
    }

    // Blow-up forensics: when dp leaves the physical band, report the
    // extreme cell's location and coefficient state (each rank reports
    // its own extreme -- cheap, only fires when already abnormal)
    {
        const scalarField& dpi = dp.primitiveField();
        label iMax = -1; scalar vMax = 0;
        forAll(dpi, i)
        {
            if (mag(dpi[i]) > vMax) { vMax = mag(dpi[i]); iMax = i; }
        }
        // Per-rank print budget: both cirius MPI_ERR_TRUNCATE deaths
        // followed dpExtreme Pout bursts -- unbounded per-step Pout from
        // many ranks stresses the collated-output communicators.  The
        // forensics value is in the FIRST events; cap the rest.
        static label nDpExtremePrints = 0;
        if (iMax >= 0 && vMax > 1e6 && nDpExtremePrints < 20)
        {
            if (++nDpExtremePrints == 20)
            {
                Pout<< "PEQSI dpExtreme: print budget reached,"
                    << " further events suppressed on this rank" << endl;
            }
            Pout<< "PEQSI dpExtreme: |dp|=" << dpi[iMax]
                << " at " << mesh.C()[iMax]
                << " T=" << thermo_.T()[iMax]
                << " rho=" << rho_[iMax]
                << " rhoEOS-drift=" << (thermo_.rho()()[iMax] - rho_[iMax])
                << " alpha=" << alpha_[iMax]
                << " 1-alpha=" << (1.0 - alpha_[iMax])
                << " beta=" << beta_[iMax]
                << " coef=" << coef[iMax]
                << " sComp=" << sComp_()[iMax]
                << " h=" << h_[iMax]
                << endl;
        }
    }

    // ------------------------------------------------------------------
    // Coupled updates with dp (WKK Eqs. 22-24; PEQSI Eqs. 11-13)
    // ------------------------------------------------------------------

    // Eq. (22): rho^{n+1} = rho^* - rho^n (1-alpha)/beta dp
    rho_ = rhoStar - coef*dp;
    rho_.correctBoundaryConditions();

    // Eq. (13): p^{n+1} = p^* + dp
    p_ = pStar + dp;
    p_.correctBoundaryConditions();

    // Eq. (23): (rho u)^{n+1} =
    //   rho^* u^* - rho^n u^n (1-alpha)/beta dp
    //   - dt grad( (p^{n+1} + p^n)/2 )
    // Face-consistent pressure gradient: reconstruct from the face-normal
    // snGrad so the momentum update sees the SAME compact stencil as the
    // dp Laplacian.  The plain cell-centred fvc::grad decouples odd/even
    // cells on this colocated arrangement (the reference is a staggered-
    // type FD code) and was observed to drive a checkerboard blow-up of
    // dp at the advected interface (peqsi1d_A, ~step 60).
    const volVectorField gradPmid
    (
        "PEQSI:gradPmid",
        fvc::reconstruct(fvc::snGrad(0.5*(p_ + pN_()))*mesh.magSf())
    );

    const volVectorField rhoUNew
    (
        "PEQSI:rhoUNew",
        rhoStar*UStar
      - coef*UN_()*dp
      - dt*gradPmid
    );

    // A lip drift episode can carry the transported rho of single cells
    // through zero while advection restores them; dividing by the raw
    // field there injects inf into u and h and the NEXT step's kappa
    // laplacian FPEs (rd0110, guarded run, step 203).  Divide by a
    // magnitude-floored rho instead: healthy cells are untouched, a
    // sick cell gets a finite (sign-consistent, wrong-anyway) value the
    // closure clamps, instead of inf.
    volScalarField rhoDen("PEQSI:rhoDen", rho_);
    {
        scalarField& rd = rhoDen.primitiveFieldRef();
        forAll(rd, celli)
        {
            if (mag(rd[celli]) < 0.01)
            {
                rd[celli] = rd[celli] < 0 ? -0.01 : 0.01;
            }
        }
        rhoDen.correctBoundaryConditions();
    }

    U_ = rhoUNew/rhoDen;
    U_.correctBoundaryConditions();

    // Eq. (24): rho^{n+1} h^{n+1} =
    //   rho^* h^* - ( rho^n h^n (1-alpha)/beta - 1 ) dp
    h_ = (rhoStar*hStar - (coef*hN_() - 1.0)*dp)/rhoDen;
    h_.correctBoundaryConditions();

    // h-budget probe, acoustic side (see advectSubstep counterpart)
    if (pimple.dict().lookupOrDefault<Switch>("peqsiStiffCensus", false))
    {
        const scalar rhAfter =
            gSum(rho_.primitiveField()*h_.primitiveField()
                *mesh.V().primitiveField());
        const scalar rhStar =
            gSum(rhoStar.primitiveField()*hStar.primitiveField()
                *mesh.V().primitiveField());
        Info<< "PEQSI h-budget: acoustic update d(int rho h) = "
            << rhAfter - rhStar << " J" << endl;
    }

    // Where the acoustic step's energy goes.  Subtracting p^{n+1} from
    // Eq. (24) gives d(rho e) = -coef h^n dp cell by cell, i.e. every
    // cell's density change carries that cell's OWN enthalpy.  For a
    // change of state that is right -- h drho is the isentropic answer.
    // For mass MOVED between cells it is not: the mass leaving A arrives
    // at B carrying h_B instead of h_A, so a redistribution that is
    // exactly mass-neutral still books (h_B - h_A) dm of energy.  The
    // two sums below separate those.  If dM is ~0 while dE is not, the
    // leak is redistribution and no time-step refinement removes it --
    // which is what the dt sweep showed (halving dt bought only 0.7x).
    if (pimple.dict().lookupOrDefault<Switch>("peqsiStiffCensus", false))
    {
        const scalarField& dV = mesh.V();
        const scalarField& cf = coef.primitiveField();
        const scalarField& dpf = dp.primitiveField();
        const scalarField& hf = hN_().primitiveField();

        if (hLedgerRef_ == GREAT)
        {
            hLedgerRef_ = gAverage(h_.primitiveField());
        }

        // Book against the manifold's own datum where it is available.
        // What the acoustic step does to the manifold is measured by
        // dh = h - h_mix(Z), because Z is untouched here: writing the
        // update as a defect density D = rho h - rho h_mix(Z) gives
        // D^{n+1} - D* = dh dm + dp, so the first term is the part
        // that moving mass creates and the second is compression, which
        // is physical but sits outside a constant-pressure envelope.
        // A single global datum is not enough -- it leaves the fuel
        // cells carrying an offset of order 1.9e6 J/kg, which is the
        // same disease in weaker form.
        // Linear, unclamped mixing line -- deliberately NOT dhF_.
        //
        // The manifold's own dh subtracts dhRef(Z, gZ, c), a tabulated
        // nonlinear function, so int rho dhF dV moves whenever gZ moves
        // even if no defect is created: mixing two parcels gives
        // mean[dhRef(Z)] != dhRef(mean[Z]), and gZ grows as the
        // interface develops.  A ledger built on it reports the
        // manifold's reference surface shifting as though the substep
        // had made defect.  The advective probe already uses the linear
        // form and the reachability envelope is defined on it, so use
        // it here too and keep all three commensurable.
        tmp<scalarField> tDhLin;
        if (fgmActive_ && fgmTable_.valid())
        {
            const scalar hOx = fgmTable_().hOx();
            const scalar hFu = fgmTable_().hFuel();
            const scalarField& Zf = Z_().primitiveField();
            const scalarField& hcell = h_.primitiveField();
            tDhLin = new scalarField(Zf.size());
            scalarField& dl = tDhLin.ref();
            forAll(dl, i)
            {
                dl[i] = hcell[i] - ((1.0 - Zf[i])*hOx + Zf[i]*hFu);
            }
        }
        const scalarField* dhp = tDhLin.valid() ? &tDhLin() : nullptr;

        const bool leakField =
            pimple.dict().lookupOrDefault<Switch>("peqsiLeakField", false);

        if (leakField && !dhLeak_.valid())
        {
            dhLeak_.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        "PEQSI:dhLeak",
                        runTime.name(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(dimEnergy, 0)
                )
            );
        }

        scalar dM = 0, dE = 0, dEabs = 0, dPV = 0;
        scalarField* lk =
            leakField ? &dhLeak_().primitiveFieldRef() : nullptr;

        forAll(dpf, i)
        {
            const scalar hb = dhp ? (*dhp)[i] : (hf[i] - hLedgerRef_);
            const scalar e = -cf[i]*hb*dpf[i]*dV[i];
            dM -= cf[i]*dpf[i]*dV[i];
            dE += e;
            dEabs += mag(cf[i]*hb*dpf[i])*dV[i];
            dPV += dpf[i]*dV[i];
            if (lk) (*lk)[i] += e;
        }
        reduce(dPV, sumOp<scalar>());
        reduce(dM, sumOp<scalar>());
        reduce(dE, sumOp<scalar>());
        reduce(dEabs, sumOp<scalar>());

        const label every =
            pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
        static label nAc = 0;
        if (every > 0 && (nAc++ % every) == 0)
        {
            // The acoustic contribution to the same defect the
            // advective probe reports, so the two are addable and the
            // source can be told from the sink.
            Info<< "PEQSI defect budget: acoustic d(int rho dh) = "
                << dE + dPV << " J" << endl;

            // And to the envelope excess, which is the quantity the
            // reachability census counts.  The advective probe leaves
            // its post-substep value in envExcessPrev_, so the
            // difference here is the acoustic step's own contribution.
            if (fgmActive_ && fgmTable_.valid() && envExcessPrev_ >= 0)
            {
                const scalar hFu = fgmTable_().hFuel();
                const scalarField& Zf = Z_().primitiveField();
                const scalarField& hc = h_.primitiveField();
                const scalarField& rc = rho_.primitiveField();
                const scalarField& Vc = mesh.V();

                scalar eNow = 0;
                label nE = 0;
                forAll(Vc, i)
                {
                    const scalar zc = min(max(Zf[i], 0.0), 1.0);
                    const scalar e =
                        hc[i] - ((1.0 - zc)*hChamber_ + zc*hFu);
                    if (e > 0) { eNow += rc[i]*e*Vc[i]; nE++; }
                }
                reduce(eNow, sumOp<scalar>());
                reduce(nE, sumOp<label>());

                Info<< "PEQSI envelope budget: acoustic dE = "
                    << eNow - envExcessPrev_ << " J (E = " << eNow
                    << " J on " << nE << " cells)" << endl;
            }

            Info<< "PEQSI acoustic ledger (datum "
                << (dhp ? "dh" : "global") << "): d(mass) = " << dM
                << " kg, dh*dm = " << dE
                << " J, dp*dV = " << dPV
                << " J, cancellation = " << dE/max(dEabs, vSmall)
                << endl;
        }
    }

    // Explicit selective filter on p and U (TK Sec. 2.7 family):
    // kills the 2-cell dp/velocity noise whose Eq. (23) kick otherwise
    // exceeds the local CFL in mature vortex fields (measured: dp 1e6
    // -> 117 m/s/step at dy = 8.2 um).  sigma = 0 disables.
    const scalar filterSigma
    (
        pimple.dict().lookupOrDefault<scalar>("peqsiFilterSigma", 0.0)
    );

    // Sensor-gated shock-capturing filter (Bogey, de Cacqueray &
    // Bailly, JCP 228 (2009) 1447): a CONSERVATIVE 2nd-order filter
    // whose local strength follows a high-pass pressure sensor, so it
    // acts only on grid-scale oscillations (the dp-kick rho/p spike)
    // and is the identity on monotone transcritical fronts -- the
    // failure mode of the un-gated variants.  sigmaSC = 0 disables.
    const scalar filterSC
    (
        pimple.dict().lookupOrDefault<scalar>("peqsiFilterSC", 0.0)
    );
    if (filterSC > 0)
    {
        applySCFilter
        (
            filterSC,
            pimple.dict().lookupOrDefault<scalar>
            (
                "peqsiFilterSCThreshold", 1e-5
            )
        );
    }

    if (filterSigma > 0)
    {
        // rho in the filter set is EXPERIMENTAL (off by default): the
        // flat mid-band response of the explicit 8th-order filter
        // damages legitimate 87:1 transcritical fronts (measured strip
        // mass leak 8e-5/0.1s through the positivity cap); the
        // production-proven set is p and U.
        applyFilter
        (
            filterSigma,
            pimple.dict().lookupOrDefault<Switch>("peqsiFilterRho", false)
        );
    }

    // End-of-step mass flux
    phi_ = fvc::flux(rho_*U_);

    mark(tPhase_[6]);

    // ------------------------------------------------------------------
    // Conservation audit, the papers' metric: cumulative relative drift
    // of int(rho) dV and int(rho h) dV from the initial state (an
    // instantaneous ddt+div residual cannot measure the error of this
    // non-conservative-form scheme -- measured on case A: residual -4e-9
    // while the mass actually grew +5.6%/period with the consistency RHS).
    // ------------------------------------------------------------------
    {
        // The audit is a diagnostic, not a control: its two gSum
        // reductions (and the closure's drift reductions) are global
        // synchronisation points, and at 256 ranks those are where the
        // imbalance wait surfaces.  Reporting every step bought nothing
        // the last report of an interval does not.
        const label diagN =
            pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);

        const bool report =
            initialMass_ < 0 || (runTime.timeIndex() % max(diagN, 1) == 0);

        if (report)
        {
        const scalar M =
            gSum(rho_.primitiveField()*mesh.V().primitiveField());
        const scalar E =
            gSum
            (
                rho_.primitiveField()*h_.primitiveField()
               *mesh.V().primitiveField()
            );

        const scalar Zm =
            fgmActive_
          ? gSum
            (
                rho_.primitiveField()*Z_().primitiveField()
               *mesh.V().primitiveField()
            )
          : 0.0;

        if (initialMass_ < 0)
        {
            initialMass_ = M;
            initialRhoH_ = E;
            initialRhoZ_ = Zm;
        }

        Info<< "PEQSI conservation: mass rel = "
            << (M - initialMass_)/initialMass_
            << ", rho*h rel = "
            << (E - initialRhoH_)/max(mag(initialRhoH_), small)
            << ", rho*Z rel = "
            << (Zm - initialRhoZ_)/max(mag(initialRhoZ_), small)
            << ", rho*Z abs = " << Zm - initialRhoZ_ << " kg"
            << " | abs: int(rho) = " << M
            << " kg, int(rho h) = " << E
            << " J, int(rho Z) = " << Zm << " kg"
            << endl;
        }
    }

    // ------------------------------------------------------------------
    // Thermodynamic closure on the end-of-step (rho, h, p) state:
    // T from h(T,v) Newton inversion, then transport properties and the
    // alpha/beta coefficient fields for the next step.
    // ------------------------------------------------------------------
    invertTemperature();
    mark(tPhase_[7]);
    updateCoefficients();
    mark(tPhase_[8]);

    if (timers && (++tPhaseSteps_ % 50 == 0))
    {
        static const char* nm[9] =
        {
            "setup", "sgs", "lad", "rk3(WENO)",
            "helmAsm", "helmSolve", "updates", "T-Newton", "coeffs"
        };
        scalar tot = 0;
        for (int i = 0; i < 9; i++) tot += tPhase_[i];
        Info<< "PEQSI timers after " << tPhaseSteps_ << " steps (cpu s, "
            << "this rank):" << nl;
        // Master-rank numbers alone were misleading on the weighted
        // decomposition: rank 0 is far-field, and its rk3 entry was
        // mostly WAIT for the jet-corridor ranks at the RK-stage halo
        // syncs.  min/max across ranks separates compute imbalance
        // (spread in the busy phase) from where the wait surfaces.
        for (int i = 0; i < 9; i++)
        {
            scalar tMin = tPhase_[i], tMax = tPhase_[i];
            reduce(tMin, minOp<scalar>());
            reduce(tMax, maxOp<scalar>());
            Info<< "    " << nm[i] << " = " << tPhase_[i]
                << " (" << 100*tPhase_[i]/max(tot, small)
                << "%)  ranks[" << tMin << ", " << tMax << "]" << nl;
        }
        extern scalar peqsiOuterAcc;
        Info<< "    accounted total = " << tot
            << "  (outer machinery between steps = " << peqsiOuterAcc
            << ")" << endl;
    }

    if (timers)
    {
        extern scalar peqsiOuterMark;
        peqsiOuterMark = runTime.elapsedCpuTime();
    }

    // Release the substep state
    rhoN_.clear();
    UN_.clear();
    pN_.clear();
    hN_.clear();
    sComp_.clear();
}


// ************************************************************************* //


void Foam::solvers::peqsiFluid::ensureDirGeometry() const
{
    // Direction weights/spacings depend only on the (static) mesh:
    // built once, cached as members
    if (!ladWDir_.valid())
    {
        const surfaceScalarField deltaF(mag(mesh.delta()));

        ladWDir_.set(new PtrList<surfaceScalarField>(3));
        ladDeltaDir_.set(new PtrList<volScalarField>(3));

        for (direction cmpt = 0; cmpt < 3; cmpt++)
        {
            vector e(Zero);
            e[cmpt] = 1;

            ladWDir_().set
            (
                cmpt,
                new surfaceScalarField
                (
                sqr((mesh.Sf()/mesh.magSf()) & e)
                )
            );

            const surfaceScalarField wA(ladWDir_()[cmpt]*mesh.magSf());

            ladDeltaDir_().set
            (
                cmpt,
                new volScalarField
                (
                [&]() -> volScalarField
                {
                    // internal-only quotient lifted to a zero-gradient
                    // boundary-complete field
                    const volScalarField::Internal q
                    (
                        fvc::surfaceSum(wA*deltaF)()
                       /max
                        (
                            fvc::surfaceSum(wA)(),
                            dimensionedScalar(dimArea, vSmall)
                        )
                    );
                    volScalarField f
                    (
                        IOobject("PEQSI:DeltaDir", mesh.time().name(), mesh),
                        mesh,
                        dimensionedScalar(q.dimensions(), 0),
                        zeroGradientFvPatchScalarField::typeName
                    );
                    f.primitiveFieldRef() = q;
                    f.correctBoundaryConditions();
                    return f;
                }()
                )
            );
        }
    }
    if (!ladDeltaMin_.valid())
    {
        // smallest INTERNAL-face centre-to-centre spacing per cell
        // (empty/boundary-only directions excluded)
        ladDeltaMin_.set(new scalarField(mesh.nCells(), great));
        scalarField& dm = ladDeltaMin_();

        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        const surfaceScalarField deltaFI(mag(mesh.delta()));
        const scalarField& df = deltaFI.primitiveField();

        forAll(own, facei)
        {
            dm[own[facei]] = min(dm[own[facei]], df[facei]);
            dm[nei[facei]] = min(dm[nei[facei]], df[facei]);
        }
    }

}


void Foam::solvers::peqsiFluid::applyFilter
(
    const scalar sigma,
    const bool filterRho
)
{
    // Explicit 8th-order selective low-pass filter (see peqsiFluid.H):
    // per direction l, q -= sigma/256 * (Delta_l^2 lap_l)^4 q, applied
    // sequentially in each direction to the TK 2012 Sec. 2.7 variable
    // set: the conservative rho and rho*u plus the pressure p (their
    // pressure-evolution form filters p, not E).  The p,U-only variant
    // failed to contain the mature-vortex dp-kick spiral at the jet-lip
    // shear layer (141449/141450, x=9.2 cm): the unfiltered density
    // spike (rho -> 4614) kept re-seeding the pressure kick.
    //
    // Conservation and boundary treatment: the face weights are zeroed
    // on all physical boundaries AND on every face of a cell that
    // touches one, so each fvc::laplacian pass telescopes to exactly
    // zero over the domain (global mass invariant to machine precision;
    // momentum is NOT -- the current set filters the primitive U, and
    // rho*(dU) does not telescope when rho varies.  That is the
    // documented Visbal-Gaitonde primitive-variable practice, its
    // momentum footprint is O(sigma/256) of grid-scale content, but the
    // old claim of momentum invariance belonged to the rho*u variable
    // set this filter no longer uses) and the filter is the IDENTITY in
    // the whole
    // boundary-adjacent cell layer.  The one-sided asymmetric stencil
    // the plain boundary-zeroing left at the inlet lip drove a
    // conservative-recovery feedback there (141451: lip cell T
    // 603 -> 1356 K, rho -> 6558 within ~500 steps); TK handle the
    // same problem with dedicated one-sided boundary filter formulas
    // -- degrading to the identity is the conservative substitute.
    // All operators are standard coupled FV ops -- processor-boundary
    // consistent by construction (coupled patches keep their exchange
    // values).
    ensureDirGeometry();

    const PtrList<surfaceScalarField>& wDir = ladWDir_();
    const PtrList<volScalarField>& DeltaDir = ladDeltaDir_();

    if (!filterMask_.valid())
    {
        filterMask_.set(filterBoundaryMask().ptr());
    }
    const surfaceScalarField& maskF = filterMask_();

    // Variable set: the PRIMITIVES p, rho, U.  Filtering the
    // conservative pair rho*u / rho*h and recovering u, h by division
    // blows up at sharp fronts (measured: strip lip T 603 -> 1356 K in
    // 141451, interior front T -> 2161 K in the strip reproducer) --
    // the transcritical density ratio makes the ratio recovery
    // ill-conditioned exactly where the filter acts.  rho must be in
    // the set: the mature-vortex dp-kick spike lives in rho (141450:
    // rho -> 4614 at x = 9.2 cm) and regenerates the pressure spike
    // through sComp every step if only p and U are filtered.  rho is
    // filtered in flux form, so global mass is preserved to machine
    // precision; h is untouched (T inversion stays on the transported
    // (h, p) state).  Primitive-variable selective filtering is
    // standard published practice (Visbal & Gaitonde).
    for (direction cmpt = 0; cmpt < 3; cmpt++)
    {
        const surfaceScalarField w(wDir[cmpt]*maskF);

        const volScalarField d2l(sqr(DeltaDir[cmpt]));

        // p
        {
            volScalarField d(p_);
            for (label pass = 0; pass < 4; pass++)
            {
                d = d2l*fvc::laplacian(w, d);
            }
            p_ -= (sigma/256.0)*d;
        }

        // rho (flux form: global mass preserved -- exactly while the
        // cap below is inactive; a capped cell breaks the pairing by
        // its excess, which the mass ledger audit makes visible).
        // Positivity guard:
        if (filterRho)
        // at an 87:1 transcritical front the 8th-order increment can
        // undershoot the gas-side cell below zero (measured: strip
        // rho -> ~0, drift 3e18, SIGFPE).  Cap the per-cell decrement
        // at 20% of the local density -- inactive in normal operation
        // (increments are O(1e-6 rho)), it only disarms the
        // pathological cell.
        {
            volScalarField d(rho_);
            for (label pass = 0; pass < 4; pass++)
            {
                d = d2l*fvc::laplacian(w, d);
            }

            scalarField& rf = rho_.primitiveFieldRef();
            const scalarField& df = d.primitiveField();
            const scalar s256 = sigma/256.0;
            forAll(rf, i)
            {
                const scalar incr = s256*df[i];
                const scalar cap = 0.2*rf[i];
                rf[i] -= min(max(incr, -cap), cap);
            }
        }

        // U
        {
            volVectorField d(U_);
            for (label pass = 0; pass < 4; pass++)
            {
                d = d2l*fvc::laplacian(w, d);
            }
            U_ -= (sigma/256.0)*d;
        }

        p_.correctBoundaryConditions();
        rho_.correctBoundaryConditions();
        U_.correctBoundaryConditions();
    }
}


Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::peqsiFluid::filterBoundaryMask() const
{
    // 0 on every face of a boundary-adjacent cell, 1 elsewhere.  The
    // indicator is exchanged across coupled patches so both ranks zero
    // a shared processor face identically.
    volScalarField bInd
    (
        IOobject("peqsiFilterBInd", mesh.time().name(), mesh),
        mesh,
        dimensionedScalar(dimless, 0)
    );
    forAll(mesh.boundary(), patchi)
    {
        if (!mesh.boundary()[patchi].coupled())
        {
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            forAll(fc, i) bInd[fc[i]] = 1;
        }
    }
    bInd.correctBoundaryConditions();

    tmp<surfaceScalarField> tmaskF
    (
        new surfaceScalarField
        (
            IOobject("peqsiFilterMask", mesh.time().name(), mesh),
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    surfaceScalarField& maskF = tmaskF.ref();

    scalarField& mIn = maskF.primitiveFieldRef();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    forAll(mIn, facei)
    {
        mIn[facei] = 1 - max(bInd[own[facei]], bInd[nei[facei]]);
    }
    forAll(maskF.boundaryFieldRef(), patchi)
    {
        scalarField& mp = maskF.boundaryFieldRef()[patchi];
        if (maskF.boundaryField()[patchi].coupled())
        {
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            const scalarField bNei
            (
                bInd.boundaryField()[patchi].patchNeighbourField()
            );
            forAll(fc, i)
            {
                mp[i] = 1 - max(bInd[fc[i]], bNei[i]);
            }
        }
        else
        {
            mp = 0;
        }
    }

    return tmaskF;
}


void Foam::solvers::peqsiFluid::applySCFilter
(
    const scalar sigmaMax,
    const scalar rTh
)
{
    // Sensor-gated conservative shock-capturing filter after Bogey,
    // de Cacqueray & Bailly, JCP 228 (2009) 1447 ("adaptative spatial
    // filtering"), FV translation:
    //
    //   per direction l:
    //     Dp      = -1/4 lap_l(p)                  (high-pass pressure)
    //     Dh      = -1/4 lap_l(h)                  (high-pass enthalpy)
    //     r       = max((Dp/p)^2, (Dh/|h|)^2)      (dimensionless sensor)
    //     sig_i   = max(0, 1 - rTh/r)              (their Eq. for sigma)
    //     q      += sigmaMax/4 div(sig_f w_l D_l^2 grad q),  q in {p, rho, U}
    //
    // Every coefficient (directional weight, boundary mask, Delta^2,
    // face sensor) is folded INTO the flux, so each pass is a pure
    // flux divergence: rho's global integral is invariant to machine
    // precision on arbitrarily graded meshes.  The 2nd-order kernel is
    // positivity-safe (increment moves a cell toward its neighbours,
    // bounded by sigmaMax/4 * sum w <= 1/2 of the local difference).
    // Distinct from the background selective filter: this one is the
    // IDENTITY wherever the sensor is quiet, so monotone transcritical
    // fronts (the failure mode of un-gated rho filtering) are never
    // touched, while the grid-scale dp-kick rho/p spike -- which
    // regenerates dp through sComp every step -- is removed at a
    // strength no CFL-limited artificial diffusivity can reach.
    ensureDirGeometry();

    const PtrList<surfaceScalarField>& wDir = ladWDir_();
    const PtrList<volScalarField>& DeltaDir = ladDeltaDir_();

    if (!filterMask_.valid())
    {
        filterMask_.set(filterBoundaryMask().ptr());
    }
    const surfaceScalarField& maskF = filterMask_();

    const dimensionedScalar rThD(dimless, rTh);
    const dimensionedScalar rFloor(dimless, vSmall);

    scalar sigMaxAll = 0;   // diagnostic: strongest sensor this step
    scalar rMaxAll = 0;     // raw sensor value, reported even when quiet

    for (direction cmpt = 0; cmpt < 3; cmpt++)
    {
        const surfaceScalarField wd2
        (
            wDir[cmpt]*maskF
           *fvc::interpolate(sqr(DeltaDir[cmpt]))
        );

        const volScalarField Dp(-0.25*fvc::laplacian(wd2, p_));
        const volScalarField Dh(-0.25*fvc::laplacian(wd2, h_));

        // The sensor triggers on EITHER carrier.  BCB gate on pressure
        // alone, which is right for the acoustic spikes they target; it
        // is blind to the mode that killed 141472 at t = 6.63 ms, where
        // one cell at the jet head went T 298 -> 1032 K in four steps
        // while the pressure field stayed smooth enough that (Dp/p)^2
        // never crossed rTh.  Adding h to the filtered SET (552889f)
        // could not help on its own: with sigma = 0 the filter is the
        // identity on every variable in the set, so the trigger has to
        // see h too.
        //
        // h is not sign-definite (it runs around -1.3e5 J/kg here), so
        // the local magnitude does the normalising, floored against a
        // global scale -- otherwise a cell whose h passes through zero
        // manufactures an unbounded sensor out of a finite wiggle.
        const dimensionedScalar hScale
        (
            "hScale",
            h_.dimensions(),
            max(gMax(mag(h_.primitiveField())), vSmall)
        );

        // Enthalpy term: OFF by default, and it should stay off.
        //
        // BCB gate on (Dp/p)^2, and that works because of pressure
        // EQUILIBRIUM: a physical contact carries no pressure jump, so
        // whatever Dp survives the high-pass is numerical.  Enthalpy has
        // no such property here -- a 125 K core in a 298 K ambient is a
        // large PHYSICAL jump in h -- so (Dh/|h|)^2 does not detect
        // wiggles, it detects the interface.  Measured on the 3-D jet:
        // the sensor went from 6 firings in 3652 steps to 58 in 59, at
        // ~0.97 strength, and the jet diffused from a 20 mm head to a
        // 60 mm smear in one millisecond.  That is the same failure the
        // |grad rho| LAD sensor showed on the Mayer case: a first-
        // derivative measure cannot tell a transcritical interface from
        // grid-scale noise.
        //
        // The mode that killed 141472 was in dp, not in h, so the
        // pressure sensor is the right instrument for it -- it just
        // fired too late.  Calibrating rTh needs the sensor's margin
        // during the growth, which is why rMax is now reported every
        // step whether or not the filter fires.
        const Switch hSensor
        (
            pimple.dict().lookupOrDefault<Switch>("peqsiFilterSCEnthalpy", false)
        );

        const volScalarField r
        (
            hSensor
          ? volScalarField
            (
                max
                (
                    sqr(Dp/p_),
                    sqr(Dh/max(mag(h_), 1e-3*hScale))
                )
            )
          : volScalarField(sqr(Dp/p_))
        );

        const scalar rMaxDir = gMax(r.primitiveField());
        rMaxAll = max(rMaxAll, rMaxDir);

        // Quiet direction: sigma = max(0, 1 - rTh/r) is identically zero
        // wherever r <= rTh, so when even the maximum is below threshold
        // every apply-laplacian below would add exactly zero -- at the
        // cost of four laplacian sweeps, two full-field products and
        // twelve boundary corrections.  Profiled on the 3-D 256-rank
        // run: the updates phase carried 41% of the step with the
        // sensor quiet at rMax 3.5e-7 against rTh 1e-5.  Skipping is
        // the identity, not an approximation.
        if (rMaxDir <= rTh)
        {
            continue;
        }

        const volScalarField sigSC
        (
            max
            (
                dimensionedScalar(dimless, 0),
                1 - rThD/max(r, rFloor)
            )
        );
        sigMaxAll = max(sigMaxAll, gMax(sigSC.primitiveField()));

        const surfaceScalarField coeff
        (
            (sigmaMax/4.0)*fvc::interpolate(sigSC)*wd2
        );

        // Filter the CONSERVED quantities, not the primitives.
        //
        // Each pass is a flux divergence, so it leaves the volume
        // integral of whatever it is applied to invariant.  Applied to
        // rho that is exactly right -- mass comes out at 1e-14.  Applied
        // to U and h it conserves the integral of U and of h, which are
        // not the conserved quantities: momentum is rho*U and energy is
        // rho*h, and in a variable-density field smoothing the primitive
        // does not conserve the product.  Measured on 1-D case A, that
        // cost 22% of the energy (rho*h rel 0.223 with the filter on
        // against 7.7e-9 with it off) -- against the paper's <5% budget.
        //
        // Filtering rho*U and rho*h and dividing by the filtered rho
        // afterwards conserves all three integrals and leaves the
        // sensor, and therefore the shock capturing, untouched.
        const volScalarField rhoh("PEQSI:rhoh", rho_*h_);
        const volVectorField rhoU("PEQSI:rhoU", rho_*U_);

        p_ += fvc::laplacian(coeff, p_);
        rho_ += fvc::laplacian(coeff, rho_);

        const volScalarField rhohNew(rhoh + fvc::laplacian(coeff, rhoh));
        const volVectorField rhoUNew(rhoU + fvc::laplacian(coeff, rhoU));

        U_ = rhoUNew/max(rho_, dimensionedScalar(rho_.dimensions(), small));
        h_ = rhohNew/max(rho_, dimensionedScalar(rho_.dimensions(), small));
        // h is in the SC set (BCB apply their shock-capturing filter to
        // the full conservative set including energy): the 141455
        // runaway at t=9.30 ms rode on h/T (rho stable at 84 while
        // T 239 -> 1864 K) -- untouchable by a p,rho,U-only set.  The
        // sensor gating keeps monotone transcritical fronts untouched
        // (the failure mode of un-gated enthalpy filtering).

        p_.correctBoundaryConditions();
        rho_.correctBoundaryConditions();
        U_.correctBoundaryConditions();
        h_.correctBoundaryConditions();
    }

    // rMax every step, sigma only when it bites: the threshold can only
    // be calibrated against how close the sensor ran to it while the
    // instability was still growing.
    Info<< "PEQSI SC sensor: rMax = " << rMaxAll
        << " (rTh = " << rTh << ")";
    if (sigMaxAll > small)
    {
        Info<< ", max sigma = " << sigMaxAll;
    }
    Info<< endl;
}
