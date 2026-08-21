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
    const volScalarField coef
    (
        "PEQSI:coef",
        rhoN_()*(1.0 - alpha_)/beta_
    );

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

    Info<< "PEQSI: dp min/max = "
        << gMin(dp.primitiveField()) << " / "
        << gMax(dp.primitiveField()) << " Pa (theta = " << theta << ")"
        << endl;

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

    U_ = rhoUNew/rho_;
    U_.correctBoundaryConditions();

    // Eq. (24): rho^{n+1} h^{n+1} =
    //   rho^* h^* - ( rho^n h^n (1-alpha)/beta - 1 ) dp
    h_ = (rhoStar*hStar - (coef*hN_() - 1.0)*dp)/rho_;
    h_.correctBoundaryConditions();

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
        const scalar M =
            gSum(rho_.primitiveField()*mesh.V().primitiveField());
        const scalar E =
            gSum
            (
                rho_.primitiveField()*h_.primitiveField()
               *mesh.V().primitiveField()
            );

        if (initialMass_ < 0)
        {
            initialMass_ = M;
            initialRhoH_ = E;
        }

        Info<< "PEQSI conservation: mass rel = "
            << (M - initialMass_)/initialMass_
            << ", rho*h rel = "
            << (E - initialRhoH_)/max(mag(initialRhoH_), small)
            << endl;
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
        for (int i = 0; i < 9; i++)
        {
            Info<< "    " << nm[i] << " = " << tPhase_[i]
                << " (" << 100*tPhase_[i]/max(tot, small) << "%)" << nl;
        }
        Info<< "    accounted total = " << tot << endl;
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
    // zero over the domain (global mass and momentum invariant to
    // machine precision) and the filter is the IDENTITY in the whole
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

    const tmp<surfaceScalarField> tmaskF(filterBoundaryMask());
    const surfaceScalarField& maskF = tmaskF();

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

    const tmp<surfaceScalarField> tmaskF(filterBoundaryMask());
    const surfaceScalarField& maskF = tmaskF();

    const dimensionedScalar rThD(dimless, rTh);
    const dimensionedScalar rFloor(dimless, vSmall);

    scalar sigMaxAll = 0;   // diagnostic: strongest sensor this step

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

        const volScalarField r
        (
            max
            (
                sqr(Dp/p_),
                sqr(Dh/max(mag(h_), 1e-3*hScale))
            )
        );

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

        p_ += fvc::laplacian(coeff, p_);
        rho_ += fvc::laplacian(coeff, rho_);
        U_ += fvc::laplacian(coeff, U_);
        // h is in the SC set (BCB apply their shock-capturing filter to
        // the full conservative set including energy): the 141455
        // runaway at t=9.30 ms rode on h/T (rho stable at 84 while
        // T 239 -> 1864 K) -- untouchable by a p,rho,U-only set.  The
        // sensor gating keeps monotone transcritical fronts untouched
        // (the failure mode of un-gated enthalpy filtering).
        h_ += fvc::laplacian(coeff, h_);

        p_.correctBoundaryConditions();
        rho_.correctBoundaryConditions();
        U_.correctBoundaryConditions();
        h_.correctBoundaryConditions();
    }

    if (sigMaxAll > small)
    {
        Info<< "PEQSI SC filter: max sensor strength = "
            << sigMaxAll << endl;
    }
}
