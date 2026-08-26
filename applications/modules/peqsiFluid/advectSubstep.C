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
#include "zeroGradientFvPatchFields.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcFlux.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcSurfaceIntegrate.H"
#include "localMax.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Between-substep (outer machinery) CPU time: set at the end of the
// acoustic substep, read at the start of the next advective substep.
// Defined in Foam::solvers -- the extern declarations sit inside member
// functions of that namespace and resolve there.
namespace solvers
{
    scalar peqsiOuterMark = -1;
    scalar peqsiOuterAcc = 0;
}

// Advective (non-conservative) operator  u . grad q  in flux form
// div(phiv q) - q div(phiv), with the frozen n-state volumetric flux.
// The scheme name keys the case fvSchemes entry (e.g. "div(phiv,rho)").
template<class Type>
static tmp<VolField<Type>> uGrad
(
    const surfaceScalarField& phiv,
    const volScalarField& divPhiv,
    const VolField<Type>& q,
    const word& scheme
)
{
    return fvc::div(phiv, q, scheme) - q*divPhiv;
}

// Wrap an internal-only property field (OF-13 thermo mu()/kappa() return
// DimensionedField) into a boundary-complete volScalarField with
// zero-gradient patches -- adequate for the smooth transport properties
// used in the explicit substep operators.
static volScalarField completeField
(
    const word& name,
    const fvMesh& mesh,
    const tmp<volScalarField::Internal>& tif
)
{
    volScalarField f
    (
        IOobject(name, mesh.time().name(), mesh),
        mesh,
        dimensionedScalar(tif().dimensions(), 0),
        zeroGradientFvPatchScalarField::typeName
    );
    f.primitiveFieldRef() = tif();
    f.correctBoundaryConditions();
    return f;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::peqsiFluid::momentumPredictor()
{
    // ------------------------------------------------------------------
    // Advective substep (PEQSI Eqs. 8-10, WKK Eqs. 15-17): explicit
    // SSP-RK3 update of rho, rho u, p, rho h with the FROZEN n-state
    // velocity (characteristic split: only the u-eigenvalue family is
    // advanced here; the acoustic family is handled implicitly in
    // pressureCorrector).  Non-conservative advective form throughout.
    // ------------------------------------------------------------------

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

    // The phase timers cover the two substeps; everything BETWEEN a
    // step's acoustic end and the next step's advective start (Courant,
    // deltaT logic, writes, function objects, log I/O) was invisible --
    // and it measured 5.2 of 17.8 s on the 2-D case, too big to leave
    // unattributed.
    if (timers)
    {
        extern scalar peqsiOuterMark, peqsiOuterAcc;
        if (peqsiOuterMark > 0)
        {
            peqsiOuterAcc += runTime.elapsedCpuTime() - peqsiOuterMark;
        }
    }

    const dimensionedScalar& dt = runTime.deltaT();

    if (rhoN_.valid())
    {
        FatalErrorInFunction
            << "advective substep entered twice without an acoustic substep "
            << "in between -- PEQSI is a fractional step, not an iterated "
            << "corrector: set nOuterCorrectors 1"
            << exit(FatalError);
    }

    // Snapshot the n-state (consumed by the acoustic substep)
    rhoN_.set(new volScalarField("PEQSI:rhoN", rho_));
    UN_.set(new volVectorField("PEQSI:UN", U_));
    if (fgmActive_)
    {
        ZN_.reset(new volScalarField("PEQSI:ZN", Z_()));
        ZvarN_.reset(new volScalarField("PEQSI:ZvarN", Zvar_()));
        YcN_.reset(new volScalarField("PEQSI:YcN", Yc_()));
    }
    pN_.set(new volScalarField("PEQSI:pN", p_));
    hN_.set(new volScalarField("PEQSI:hN", h_));

    // Frozen advecting volumetric flux and its divergence
    const surfaceScalarField phiv("PEQSI:phiv", fvc::flux(UN_()));
    const volScalarField divPhiv("PEQSI:divPhiv", fvc::div(phiv));

    mark(tPhase_[0]);

    // Frozen sources at the n-state:
    //
    // L_h: enthalpy-equation RHS without Dp/Dt (WKK Eq. 3) -- here the
    // conductive part div(kappa grad T); viscous heating is neglected as
    // in the reference low-Mach cryogenic-jet applications.
    volScalarField kappa
    (
        completeField("PEQSI:kappa", mesh, thermo_.kappa())
    );
    volScalarField mu
    (
        completeField("PEQSI:mu", mesh, thermo_.mu())
    );

    // Subgrid contribution (V3, paper 3-D LES: SGS model with SGS
    // Prandtl number 0.7 for the subgrid heat flux).  The base class's
    // momentumTransport is corrected at the n-state; with
    // simulationType laminar nut == 0 and this is a no-op.
    if (sgsActive_ < 0)
    {
        sgsActive_ =
            word(momentumTransport->lookup("simulationType")) != "laminar";
    }

    if (sgsActive_)
    {
        momentumTransport->correct();

        const scalar prSgs
        (
            pimple.dict().lookupOrDefault<scalar>("peqsiPrSGS", 0.7)
        );

        const volScalarField muSgs
        (
            "PEQSI:muSgs",
            rhoN_()*momentumTransport->nut()
        );

        mu += muSgs;
        kappa += muSgs*thermo_.Cp()/prSgs;
    }

    mark(tPhase_[1]);

    // (L_h is assembled after the LAD block: TK Eq. 33 augments the
    // conductivity with kappa_art)

    // Viscous stress divergence for the momentum substep (WKK Eq. 16),
    // explicit at the n-state.  Molecular part only for now; the LES
    // subgrid contribution is added in the V3 stage.
    const volVectorField divTau
    (
        "PEQSI:divTau",
        fvc::laplacian(mu, UN_())
      + fvc::div(mu*dev2(T(fvc::grad(UN_()))))
    );

    // Real-gas coefficient combinations, frozen at the n-state
    // (alpha_, beta_ were evaluated at the end of the previous step)
    const volScalarField aL("PEQSI:aL", alpha_/(1.0 - alpha_));
    const volScalarField iL("PEQSI:iL", 1.0/(1.0 - alpha_));

    // Terashima-Koshi LAD, TAKEN VERBATIM from TK JCP 231 (2012) 6907,
    // Sec. 2.5 (the method PEQSI Sec. III B runs at C_rho = 0.002 -- the
    // smallest value of TK's own 1-D sweep {0.002, 0.01, 0.05}; TK's 2-D
    // jet uses 0.02, minimum-required 0.005):
    //
    //   mass     (Eq. 22-23):  + div( rho_art grad(rho) ),
    //       rho_art = C_rho (c/rho) sum_l |d^4 rho_bar/dx_l^4| Dl^5
    //   momentum (Eq. 24)   :  + div( rho_art (u (x) g) grad(rho) )
    //   energy   (Eq. 33)   :  kappa += kappa_art,
    //       kappa_art = C_kappa (rho c^3/T^2) sum_l |d^4 T_bar/dx_l^4| Dl^5
    //       with C_kappa = 0.01 (TK's standard)
    //   pressure (Eq. 30)   :  NO artificial term (pressure-equilibrium
    //                          preservation) -- p* stays untouched.
    //
    // FV translation (documented deviations, isotropic on hex meshes):
    //   - directional 4th derivatives -> |lap(lap(q_bar))| * Delta^5,
    //     Delta = V^(1/3)
    //   - TK's truncated Gaussian filter -> 2nd-order equivalent
    //     q_bar = q + (Delta^2/24) lap(q)
    //   - Eq. 24 discretised as u^n * [mass-diffusion increment]: exact
    //     under velocity equilibrium, which is the property Eq. 24 exists
    //     to preserve.
    //
    // NOTE: heap-allocate INTO the tmp -- assigning a stack-local field
    // to a tmp creates a non-owning reference that dangles out of scope
    // (measured: freed name() read -> std::length_error in fvc name
    // concatenation).
    const scalar ladCoeff
    (
        pimple.dict().lookupOrDefault<scalar>("peqsiLADCoeff", scalar(0))
    );
    const scalar ladKappaCoeff
    (
        pimple.dict().lookupOrDefault<scalar>
        (
            "peqsiLADKappaCoeff",
            ladCoeff > 0 ? scalar(0.01) : scalar(0)
        )
    );

    // peqsiLADKappaCoeff defaults to 0.01 whenever the mass LAD is on,
    // which is harmless for a single-species case but not for a
    // manifold one: the artificial conductivity and the artificial
    // mass diffusivity both act on the enthalpy side (kappa, and the
    // laplacians of rho and rho*h) while the mixture fraction gets
    // neither.  A manifold assumes h and Z move together, and every
    // route that diffuses one without the other pushes cells to states
    // no stream combination can reach -- the failure that ended the
    // rd0110 run at t = 404 us.  Enabling the mass LAD here therefore
    // switches on a third such route silently, so say so.
    if
    (
        fgmActive_ && (ladKappaCoeff > 0 || ladCoeff > 0)
     && pimple.dict().lookupOrDefault<Switch>("peqsiLADWarn", true)
    )
    {
        static bool warned = false;
        if (!warned)
        {
            warned = true;
            WarningInFunction
                << "LAD active with the manifold (peqsiLADCoeff = "
                << ladCoeff << ", peqsiLADKappaCoeff = "
                << ladKappaCoeff << "): both terms diffuse the enthalpy"
                << " side and neither touches Z -- the conductivity"
                << " through kappa, and the mass coefficient through"
                << " the laplacians of rho AND rho*h, so zeroing"
                << " peqsiLADKappaCoeff alone does not remove the"
                << " asymmetry.  It drives cells off the manifold"
                << " (peqsiLADWarn no silences this)." << endl;
        }
    }

    tmp<volScalarField> tDart;      // rho_art [m2/s]
    tmp<volScalarField> tKappaArt;  // kappa_art [W/(m K)]

    if (ladCoeff > 0 || ladKappaCoeff > 0)
    {
        // TK's length scale is DIRECTIONAL (Eqs. 22-24, 33: the
        // diffusivity is a SCALAR whose magnitude sums the per-
        // direction terms |d4 q/dx_l^4| Delta_l^5).  Two isotropic
        // shortcuts both fail on anisotropic grids (measured):
        //   - |lap(lap q))| * V^(5/3): overestimates by
        //     (V^(1/3)/dy_min)^5 = 243x on the paper's 17.5 x 8.2 um
        //     jet grid -> LAD diffusion number 6.4, blow-up in 4 steps;
        //   - the isotropic biharmonic times the FACE spacing^5:
        //     charges the y-front's 4th derivative to the coarse
        //     x-faces (38x overshoot on the smoke grid).
        // Faithful translation on a Cartesian grid: directional second
        // derivatives via direction-weighted Laplacians,
        //   d2_l(q) = lap(w_l, q),  w_l = (n_f . e_l)^2,
        // applied twice for |d4_l q|, each weighted by the per-cell
        // spacing Delta_l^5 (face-normal spacing averaged per cell).

        ensureDirGeometry();

        const scalarField& deltaMin = ladDeltaMin_();
        const PtrList<surfaceScalarField>& wDir = ladWDir_();
        const PtrList<volScalarField>& DeltaDir = ladDeltaDir_();

        // TK truncated-Gaussian filter, DIRECTIONAL form:
        //   qBar = q + sum_l (Delta_l^2/24) d2q/dx_l^2.
        // The isotropic V^(1/3) filter length re-introduces the
        // anisotropy bug through the back door: on the paper grid its
        // (V^(1/3)/dy)^2 = 8.8x overestimated correction plants a fake
        // wiggle at the 3-cell erf inlet band whose 4th derivative then
        // dominates kappa_art (measured 1.4 W/mK where the physical
        // scale is ~5e-6) -> p* spike -> dp -8.7e7 Pa in ONE step.
        auto dirFilter = [&](const volScalarField& q) -> tmp<volScalarField>
        {
            tmp<volScalarField> tq(new volScalarField(q));
            for (direction cmpt = 0; cmpt < 3; cmpt++)
            {
                tq.ref() +=
                    sqr(DeltaDir[cmpt])/24.0
                   *fvc::laplacian(wDir[cmpt], q);
            }
            return tq;
        };

        // Directional |d4 q/dx_l^4| Delta_l^5 sum for a filtered field
        auto d4Sum = [&](const volScalarField& qBar) -> tmp<volScalarField>
        {
            tmp<volScalarField> tSum;
            for (direction cmpt = 0; cmpt < 3; cmpt++)
            {
                const volScalarField d2
                (
                    fvc::laplacian(wDir[cmpt], qBar)
                );
                tmp<volScalarField> term
                (
                    mag(fvc::laplacian(wDir[cmpt], d2))
                   *pow(DeltaDir[cmpt], 5)
                );
                if (tSum.valid()) tSum.ref() += term();
                else tSum = term;
            }
            return tSum;
        };

        // sound speed from the coefficient fields already in hand:
        // c^2 = -beta/((1-alpha) rho)  (WKK App. D identity, verified to
        // 5e-16 at every step) -- thermo.gamma() would re-run BOTH the
        // Cp and Cv SRK sweeps every step just for this
        const volScalarField c
        (
            "PEQSI:c",
            sqrt(max(-beta_/((1.0 - alpha_)*rhoN_()),
                     dimensionedScalar(sqr(dimVelocity), 0)))
        );

        if (ladCoeff > 0)
        {
            // filtered density (TK truncated Gaussian, directional)
            const volScalarField rhoBar(dirFilter(rhoN_()));

            tDart = tmp<volScalarField>
            (
                new volScalarField
                (
                    "PEQSI:Dart",
                    ladCoeff*(c/rhoN_())*d4Sum(rhoBar)
                )
            );

            // The artificial diffusivity is an INTERIOR discontinuity
            // regularisation: zero it on domain boundary faces.  An
            // inflow face carries a rho jump that is a boundary
            // condition, not a captured discontinuity; letting the mass
            // diffusion act through it injects boundary mass whose
            // Eq. (24) momentum consistency (u_cell x massArt) is at its
            // least valid point.  Measured on the 2-D jet: inlet-cell
            // pressure collapse (4 -> 1.7 MPa by step 5) and T
            // undershoot to the Newton floor, blow-up by step ~15; with
            // the boundary flux removed the startup is smooth.
            // Domain boundaries only.  A processor or cyclic patch is an
            // INTERIOR face of the global domain: zeroing the exchanged
            // neighbour values there makes the face diffusivity depend
            // on the decomposition -- serial and parallel runs diverge
            // at every processor boundary the LAD crosses.
            forAll(tDart.ref().boundaryFieldRef(), patchi)
            {
                if (tDart.ref().boundaryField()[patchi].coupled()) continue;
                tDart.ref().boundaryFieldRef()[patchi] == 0.0;
            }

            // Explicit-stability cap on the ARTIFICIAL diffusivity
            // (numerics hygiene on an artificial device): at clamp-
            // saturated front cells the TK formula can spiral (measured
            // rho_art 8e9 m^2/s in the 2-D vortex blow-up).  The cap
            // (diffusion number 0.45 on the smallest internal spacing)
            // sits ABOVE normal LAD operation (~0.3 measured) and binds
            // only at pathological cells.
            {
                scalarField& D = tDart.ref().primitiveFieldRef();
                const scalar rdt = 0.45/dt.value();

                // Time-step limit from the UNCAPPED demand -- it must be
                // read before the min() below, or the capped value feeds
                // back into dt and the bound degenerates to dt <= dt.
                //
                // A single clamp-saturated cell can spiral the TK demand
                // to 8e9 m^2/s (measured, 2-D vortex blow-up); such a
                // cell must not crush dt to nothing, so the demand used
                // for the limit is itself bounded (peqsiLadDtDemandMax,
                // default 5e-3 m^2/s) and anything above that is left to
                // the 0.45 cap to handle, exactly as today.  The safety
                // factor covers the one-step lag: the limit computed
                // here steers the NEXT step's dt.
                const scalar Dsane =
                    pimple.dict().lookupOrDefault<scalar>
                    (
                        "peqsiLadDtDemandMax", 5e-3
                    );
                const scalar fSafe =
                    pimple.dict().lookupOrDefault<scalar>
                    (
                        "peqsiLadDtSafety", 0.8
                    );

                scalar dtLim = great;
                forAll(D, i)
                {
                    const scalar Duse = min(D[i], Dsane);
                    if (Duse > vSmall)
                    {
                        dtLim =
                            min(dtLim, 0.45*sqr(deltaMin[i])/Duse);
                    }
                    D[i] = min(D[i], rdt*sqr(deltaMin[i]));
                }
                reduce(dtLim, minOp<scalar>());
                ladDtLimit_ = fSafe*dtLim;
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
                    Info<< "PEQSI LAD: rho_art max = "
                        << gMax(tDart().primitiveField()) << " m^2/s"
                        << endl;
                }
            }
        }

        if (ladKappaCoeff > 0)
        {
            const volScalarField& T = thermo_.T();
            const volScalarField TBar(dirFilter(T));

            tKappaArt = tmp<volScalarField>
            (
                new volScalarField
                (
                    "PEQSI:kappaArt",
                    ladKappaCoeff*rhoN_()*pow3(c)/sqr(T)*d4Sum(TBar)
                )
            );

            // Interior regularisation only (same rationale as Dart):
            // the molecular kappa keeps its boundary flux
            // Domain boundaries only (see the Dart loop above).
            forAll(tKappaArt.ref().boundaryFieldRef(), patchi)
            {
                if (tKappaArt.ref().boundaryField()[patchi].coupled())
                {
                    continue;
                }
                tKappaArt.ref().boundaryFieldRef()[patchi] == 0.0;
            }

            // Same explicit-stability cap on kappa_art/(rho cp)
            // (cp from the previous step's coefficient evaluation;
            // measured kappa_art 1e47 W/(m K) in the blow-up)
            if (cpPrev_.size() == mesh.nCells())
            {
                scalarField& K = tKappaArt.ref().primitiveFieldRef();
                const scalarField& rhoNi = rhoN_().primitiveField();
                const scalar rdt = 0.45/dt.value();
                forAll(K, i)
                {
                    K[i] =
                        min(K[i], rdt*sqr(deltaMin[i])*rhoNi[i]*cpPrev_[i]);
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
                    Info<< "PEQSI LAD: kappa_art max = "
                        << gMax(tKappaArt().primitiveField()) << " W/(m K)"
                        << endl;
                }
            }
        }
    }

    // L_h: enthalpy-equation RHS without Dp/Dt (WKK Eq. 3) -- conduction
    // with the TK Eq. (33) artificial conductivity when active.
    volScalarField Lh
    (
        "PEQSI:Lh",
        tKappaArt.valid()
      ? fvc::laplacian((kappa + tKappaArt())(), thermo_.T())()
      : fvc::laplacian(kappa, thermo_.T())()
    );

    // The other half of the multicomponent heat flux.
    //
    // q = -kappa gradT + sum_k h_k j_k, and only the first term was
    // here.  The species fluxes j_k = -rho D grad Y_k are not
    // hypothetical: the mixture fraction is diffused with exactly that
    // coefficient a few hundred lines below, so fuel mass leaves a cell
    // carrying h_fu and oxidiser arrives carrying h_ox while the
    // enthalpy equation is told nothing about it.  With h_fu - h_ox of
    // order 1.8 MJ/kg across these two streams, moving 5% of a cell's
    // mixture fraction is worth ~90 kJ/kg of enthalpy -- which is the
    // size of the defect the reachability census has been reporting.
    //
    // Written from the transported Y rather than from stream
    // compositions, so it stays correct once the progress variable is
    // active and the composition is no longer linear in Z.  Species
    // that are absent are skipped, which on the non-reacting case
    // leaves four of thirty.
    if
    (
        pimple.dict().lookupOrDefault<Switch>
        (
            "peqsiSpeciesEnthalpyFlux", true
        )
    )
    {
        const scalar leZh =
            pimple.dict().lookupOrDefault<scalar>("peqsiLe", 1.0);
        const scalar scSgsh =
            pimple.dict().lookupOrDefault<scalar>("peqsiScSGS", 0.7);

        tmp<volScalarField> tRhoD;
        if (leZh > 0)
        {
            const volScalarField kmol
            (
                completeField("PEQSI:kappaMolH", mesh, thermo_.kappa())
            );
            tRhoD = (kmol/(thermo_.Cp()*leZh)).ptr();
        }
        if (sgsActive_ && scSgsh > 0)
        {
            const volScalarField mus
            (
                "PEQSI:muSgsH", rhoN_()*momentumTransport->nut()
            );
            if (tRhoD.valid()) tRhoD.ref() += mus/scSgsh;
            else tRhoD = mus/scSgsh;
        }

        if (tRhoD.valid())
        {
            const surfaceScalarField rhoDf
            (
                "PEQSI:rhoDf", fvc::interpolate(tRhoD())
            );
            surfaceScalarField jh
            (
                "PEQSI:jh",
                0.0*rhoDf*mesh.magSf()
               *dimensionedScalar(dimEnergy/dimMass/dimLength, 1.0)
            );

            // Referenced to the mixture enthalpy, not used raw.  The
            // sum is physically invariant to a common datum because
            // sum_k grad Y_k = 0, but only if that identity holds in
            // the discrete field -- and it does not, both because the
            // manifold's mass fractions carry interpolation error and
            // because the loop below skips absent species.  With
            // absolute enthalpies of order 1e6 J/kg, a residual of
            // 1e-6 in the sum of gradients is already worth joules per
            // kilogram, and the first attempt at this term drove the
            // temperature to its floor within thirty steps.
            // Subtracting h makes the cancellation exact by
            // construction rather than by hope.
            const PtrList<volScalarField>& Yk = thermo_.Y();
            forAll(Yk, k)
            {
                if (gMax(Yk[k].primitiveField()) < 1e-6) continue;

                const volScalarField hRel
                (
                    "PEQSI:hRel",
                    thermo_.hai(k, p_, thermo_.T()) - h_
                );

                jh +=
                    fvc::interpolate(hRel)
                   *rhoDf*fvc::snGrad(Yk[k])()*mesh.magSf();
            }

            Lh += fvc::div(jh);
        }
    }

    mark(tPhase_[2]);

    // RK working fields (rho h transported as a product)
    volScalarField r("PEQSI:rkRho", rhoN_());
    volVectorField ru("PEQSI:rkRhoU", rhoN_()*UN_());
    volScalarField pw("PEQSI:rkP", pN_());
    volScalarField rh("PEQSI:rkRhoH", rhoN_()*hN_());

    // Transport carrier: ALL transported quantities (rho, p, rho h,
    // and the three components of rho u) packed as one symmTensor field
    // so a SINGLE WENO reconstruction pass serves every equation.  The
    // WENO weights are computed per component independently (WENOCoeff
    // loops nComp inside the stencil loop), so the math is identical to
    // per-field passes -- but the per-cell stencil matrices are loaded
    // from memory ONCE instead of four times.  Profiling: the WENO pass
    // is memory-bandwidth bound (86% of step time; an AVX2 rebuild
    // gained nothing), so cutting the matrix sweeps 4 -> 1 per stage
    // attacks the actual bottleneck (measured 1.33x total already at
    // 3-scalars-only packing; bitwise-identical results).
    //
    // Component map: XX=rho, XY=p, XZ=rho*h, YY/YZ/ZZ=rho*u (x,y,z).
    //
    // Boundary handling: physical patches are 'calculated' and receive
    // copies of the component fields' own evaluated boundary values each
    // stage (identical face values to what the per-field div would use);
    // constraint patches (processor/empty/cyclic) keep their coupled
    // type so parallel interpolation is untouched.
    wordList qBcTypes(mesh.boundary().size());
    forAll(mesh.boundary(), patchi)
    {
        const word& pt = mesh.boundaryMesh()[patchi].type();
        qBcTypes[patchi] =
            (pt == "patch" || pt == "wall")
          ? word("calculated")
          : p_.boundaryField()[patchi].type();  // processor/empty/cyclic
    }

    volSymmTensorField Q
    (
        IOobject("PEQSI:Q", runTime.name(), mesh),
        mesh,
        dimensionedSymmTensor(dimless, Zero),
        qBcTypes
    );

    auto packQ = [&]()
    {
        symmTensorField& Qi = Q.primitiveFieldRef();
        const scalarField& ri = r.primitiveField();
        const scalarField& pi = pw.primitiveField();
        const scalarField& rhi = rh.primitiveField();
        const vectorField& rui = ru.primitiveField();
        forAll(Qi, i)
        {
            Qi[i] = symmTensor
            (
                ri[i], pi[i], rhi[i],
                rui[i].x(), rui[i].y(),
                rui[i].z()
            );
        }
        forAll(Q.boundaryFieldRef(), patchi)
        {
            fvPatchSymmTensorField& Qp = Q.boundaryFieldRef()[patchi];
            if (Qp.coupled()) continue;
            const fvPatchScalarField& rp = r.boundaryField()[patchi];
            const fvPatchScalarField& pp = pw.boundaryField()[patchi];
            const fvPatchScalarField& rhp = rh.boundaryField()[patchi];
            const fvPatchVectorField& rup = ru.boundaryField()[patchi];
            forAll(Qp, i)
            {
                Qp[i] = symmTensor
                (
                    rp[i], pp[i], rhp[i],
                    rup[i].x(), rup[i].y(),
                    rup[i].z()
                );
            }
        }
        Q.correctBoundaryConditions();
    };

    // One SSP-RK3 stage: q <- cOld*qn + cNew*(q + dt*L(q)).
    // The updates are FUSED into primitive loops (internal + physical
    // boundaries): the field-expression form allocated ~25 full-field
    // temporaries per stage, ~30% of the non-WENO substep cost.  The
    // per-element arithmetic and its order are identical.
    const scalar dtv = dt.value();

    auto stage = [&](const scalar cOld, const scalar cNew)
    {
        packQ();
        const volSymmTensorField LQ
        (
            -uGrad(phiv, divPhiv, Q, "div(phiv,rho)")
        );

        // LAD mass/energy diffusion increments (TK Eqs. 22, 24, 31 --
        // see the comment block above)
        tmp<volScalarField> tMassArt, tRhArt;
        if (tDart.valid())
        {
            tMassArt = fvc::laplacian(tDart(), r);
            tRhArt = fvc::laplacian(tDart(), rh);
        }

        // ---- internal cells ----
        {
            const symmTensorField& LQi = LQ.primitiveField();
            const scalarField& Lhi = Lh.primitiveField();
            const scalarField& aLi = aL.primitiveField();
            const scalarField& iLi = iL.primitiveField();
            // FPV composition source S_Y: alpha and beta are
            // fixed-composition coefficients, so without this the
            // pressure equation never learns that the manifold moved Y
            // and heat release produces no dilatation at all.  Null
            // whenever the manifold is inactive (single-species), and
            // the branch below then adds nothing -- byte-for-byte the
            // previous arithmetic.
            const scalarField* sYi =
                sourceP_.valid() ? &sourceP_().primitiveField() : nullptr;
            const vectorField& dTi = divTau.primitiveField();
            const scalarField& rhoNi = rhoN_().primitiveField();
            const vectorField& UNi = UN_().primitiveField();
            const scalarField& pNi = pN_().primitiveField();
            const scalarField& hNi = hN_().primitiveField();
            const scalarField* mAi =
                tMassArt.valid() ? &tMassArt().primitiveField() : nullptr;
            const scalarField* rAi =
                tRhArt.valid() ? &tRhArt().primitiveField() : nullptr;

            scalarField& ri = r.primitiveFieldRef();
            vectorField& rui = ru.primitiveFieldRef();
            scalarField& pwi = pw.primitiveFieldRef();
            scalarField& rhi = rh.primitiveFieldRef();

            forAll(ri, i)
            {
                const symmTensor& l = LQi[i];

                scalar Lr = l.xx();
                scalar Lp = l.xy() + aLi[i]*Lhi[i];
                scalar Lrh = l.xz() + iLi[i]*Lhi[i];
                if (sYi)
                {
                    // S_Y enters BOTH.  rho Dh/Dt = Dp/Dt + Lh, which is
                    // where Lrh's iL comes from in the first place:
                    // aL + 1 = alpha/(1-alpha) + 1 = 1/(1-alpha) = iL.
                    // So any term added to Dp/Dt appears unchanged in the
                    // enthalpy equation, and adding it to Lp alone
                    // silently drains h -- measured at -363 kJ/kg, 31%,
                    // on the constant-pressure reacting case, against a
                    // pressure impulse of ~2 MPa over rho ~ 4, i.e. the
                    // same 500 kJ/kg to within the transient.
                    const scalar sy = iLi[i]*(*sYi)[i];
                    Lp += sy;
                    Lrh += sy;
                }
                vector Lru
                (
                    l.yy() + dTi[i].x(),
                    l.yz() + dTi[i].y(),
                    l.zz() + dTi[i].z()
                );

                if (mAi)
                {
                    Lr += (*mAi)[i];
                    Lru += UNi[i]*(*mAi)[i];
                    Lrh += (*rAi)[i];
                }

                ri[i] = cOld*rhoNi[i] + cNew*(ri[i] + dtv*Lr);
                rui[i] = cOld*(rhoNi[i]*UNi[i]) + cNew*(rui[i] + dtv*Lru);
                pwi[i] = cOld*pNi[i] + cNew*(pwi[i] + dtv*Lp);
                rhi[i] = cOld*(rhoNi[i]*hNi[i]) + cNew*(rhi[i] + dtv*Lrh);
            }
        }

        // ---- physical boundary values (read back by packQ and by the
        //      boundary contributions of div/laplacian next stage) ----
        forAll(r.boundaryFieldRef(), patchi)
        {
            if (r.boundaryField()[patchi].coupled()) continue;

            const symmTensorField& LQb = LQ.boundaryField()[patchi];
            const scalarField& Lhb = Lh.boundaryField()[patchi];
            const scalarField& aLb = aL.boundaryField()[patchi];
            const scalarField& iLb = iL.boundaryField()[patchi];
            const vectorField& dTb = divTau.boundaryField()[patchi];
            const scalarField& rhoNb = rhoN_().boundaryField()[patchi];
            const vectorField& UNb = UN_().boundaryField()[patchi];
            const scalarField& pNb = pN_().boundaryField()[patchi];
            const scalarField& hNb = hN_().boundaryField()[patchi];
            const scalarField* mAb =
                tMassArt.valid()
              ? &tMassArt().boundaryField()[patchi] : nullptr;
            const scalarField* rAb =
                tRhArt.valid()
              ? &tRhArt().boundaryField()[patchi] : nullptr;

            const scalarField* sYb =
                sourceP_.valid()
              ? &sourceP_().boundaryField()[patchi] : nullptr;

            fvPatchScalarField& rb = r.boundaryFieldRef()[patchi];
            fvPatchVectorField& rub = ru.boundaryFieldRef()[patchi];
            fvPatchScalarField& pwb = pw.boundaryFieldRef()[patchi];
            fvPatchScalarField& rhb = rh.boundaryFieldRef()[patchi];

            forAll(rb, i)
            {
                const symmTensor& l = LQb[i];

                scalar Lr = l.xx();
                scalar Lp = l.xy() + aLb[i]*Lhb[i];
                scalar Lrh = l.xz() + iLb[i]*Lhb[i];
                if (sYb)
                {
                    const scalar sy = iLb[i]*(*sYb)[i];
                    Lp += sy;
                    Lrh += sy;
                }
                vector Lru
                (
                    l.yy() + dTb[i].x(),
                    l.yz() + dTb[i].y(),
                    l.zz() + dTb[i].z()
                );

                if (mAb)
                {
                    Lr += (*mAb)[i];
                    Lru += UNb[i]*(*mAb)[i];
                    Lrh += (*rAb)[i];
                }

                rb[i] = cOld*rhoNb[i] + cNew*(rb[i] + dtv*Lr);
                rub[i] = cOld*(rhoNb[i]*UNb[i]) + cNew*(rub[i] + dtv*Lru);
                pwb[i] = cOld*pNb[i] + cNew*(pwb[i] + dtv*Lp);
                rhb[i] = cOld*(rhoNb[i]*hNb[i]) + cNew*(rhb[i] + dtv*Lrh);
            }
        }

        r.correctBoundaryConditions();
        ru.correctBoundaryConditions();
        pw.correctBoundaryConditions();
        rh.correctBoundaryConditions();

        if (fgmActive_)
        {
            advectManifoldStage(phiv, divPhiv, cOld, cNew, dtv);
        }
    };

    stage(0.0, 1.0);            // q1 = qn + dt L(qn)
    const volScalarField r1("PEQSI:r1", r);
    stage(0.75, 0.25);          // q2 = 3/4 qn + 1/4 (q1 + dt L(q1))
    const volScalarField r2("PEQSI:r2", r);
    stage(1.0/3.0, 2.0/3.0);    // q* = 1/3 qn + 2/3 (q2 + dt L(q2))

    // The substep's own compression bookkeeping: SSP-RK3's effective
    // quadrature q* = q^n + dt (1/6 L0 + 1/6 L1 + 2/3 L2) applied to the
    // +rho div(phiv) part of the advective operator.  The acoustic
    // substep uses this as the Helmholtz source so that the Eq. (22)
    // density increment cancels the substep's mass change identically.
    sComp_.set
    (
        new volScalarField
        (
            "PEQSI:sComp",
            (rhoN_()/6.0 + r1/6.0 + 2.0*r2/3.0)*divPhiv
        )
    );

    // Publish the starred state into the solver fields (BCs from the
    // registered fields' own conditions)
    rho_ = r;
    rho_.correctBoundaryConditions();

    // Magnitude-floored divisor: a drift-episode cell crossing rho = 0
    // must not inject inf into u and h (see the Eq. 23/24 counterpart
    // in correctPressurePEQSI.C)
    volScalarField rDen("PEQSI:rDen", r);
    {
        scalarField& rd = rDen.primitiveFieldRef();
        forAll(rd, celli)
        {
            if (mag(rd[celli]) < 0.01)
            {
                rd[celli] = rd[celli] < 0 ? -0.01 : 0.01;
            }
        }
        rDen.correctBoundaryConditions();
    }

    U_ = ru/rDen;
    U_.correctBoundaryConditions();

    p_ = pw;
    p_.correctBoundaryConditions();

    h_ = rh/rDen;
    h_.correctBoundaryConditions();

    // ------------------------------------------------------------------
    // Mesh-general front regularisation (peqsiAVCoeff > 0 to enable).
    //
    // Jameson-type 2nd-order artificial dissipation with a density
    // front sensor -- the unstructured-FV standard (JST; the AVBP
    // real-gas LES line applies exactly this sensor/operator family to
    // transcritical coaxial injection).  Replaces the retired
    // DIRECTIONAL devices (peqsiLADCoeff shear, peqsiFilterSigma/SC,
    // peqsiBoundZ), whose grid-aligned stencils were measured hostile
    // on the 42%-skewed rd0110 mesh: here the operator is the standard
    // face-based conservative laplacian (non-orthogonal correction and
    // coupled patches included), so the discretisation is mesh-general
    // by construction.
    //
    //   sensor  psi_P = |sum_f (rho_N - rho_P)| / sum_f (rho_N + rho_P)
    //           (normalised undivided second difference: ~O(1) across a
    //           2-cell front, ~0 on smooth fields AND on linear shear,
    //           so resolved K-H at the lip is not damped)
    //   coeff   Gamma_f = k2 psi_f (|u_f| + c_f) |d_f|,  c = 1/sqrt(psi_T)
    //           (isothermal sound speed from the thermo compressibility:
    //           within sqrt(gamma) of c and positive for any state)
    //   update  q <- q + dt div(Gamma grad q) on the CONSERVATIVE set
    //           {rho, rho u, rho h, rho Z, rho Zvar, rho Yc}: at a
    //           uniform-(u, h, p) contact the same operator on rho and
    //           rho q leaves u, h, p untouched -- the substep's
    //           pressure-equilibrium property survives, and the flux
    //           form conserves every integral (the Helmholtz mass
    //           telescoping is unaffected).
    //
    // p is NOT smoothed: the acoustic substep owns it, and at a front
    // p is uniform anyway.
    const scalar k2
    (
        pimple.dict().lookupOrDefault<scalar>("peqsiAVCoeff", 0.0)
    );
    if (k2 > 0)
    {
        // Undivided second difference of rho, accumulated by a plain
        // face loop (processor/cyclic neighbours through
        // patchNeighbourField; one-sided physical patches contribute
        // nothing -- the front never needs sensing ON a wall face)
        scalarField num(mesh.nCells(), 0.0);
        scalarField den(mesh.nCells(), 0.0);
        {
            const labelUList& own = mesh.owner();
            const labelUList& nei = mesh.neighbour();
            const scalarField& rf = rho_.primitiveField();
            forAll(own, facei)
            {
                const scalar d = rf[nei[facei]] - rf[own[facei]];
                const scalar a = rf[nei[facei]] + rf[own[facei]];
                num[own[facei]] += d;
                num[nei[facei]] -= d;
                den[own[facei]] += a;
                den[nei[facei]] += a;
            }
            forAll(rho_.boundaryField(), patchi)
            {
                const fvPatchScalarField& prho =
                    rho_.boundaryField()[patchi];
                if (!prho.coupled()) continue;
                const labelUList& fc = prho.patch().faceCells();
                const tmp<scalarField> tnf(prho.patchNeighbourField());
                const scalarField& nf = tnf();
                forAll(fc, facei)
                {
                    num[fc[facei]] += nf[facei] - rf[fc[facei]];
                    den[fc[facei]] += nf[facei] + rf[fc[facei]];
                }
            }
        }
        volScalarField psiS
        (
            IOobject("PEQSI:avSensor", runTime.name(), mesh),
            mesh,
            dimensionedScalar(dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );
        {
            scalarField& pf = psiS.primitiveFieldRef();
            forAll(pf, celli)
            {
                pf[celli] =
                    min(mag(num[celli])/max(den[celli], small), 1.0);
            }
        }
        psiS.correctBoundaryConditions();

        const volScalarField lambda
        (
            mag(U_) + sqrt(1.0/max(thermo_.psi(), dimensionedScalar(dimDensity/dimPressure, vSmall)))
        );

        // Face value = MAX of the two sides (the JST/AVBP convention for
        // both the sensor and the spectral radius), not the linear
        // average: at a liquid/gas mixed face the liquid side's c ~ 1300
        // then sets the scale, so the dissipation is strongest exactly at
        // the transcritical faces where the dp/c^2 feedback lives --
        // averaging halved it there.
        const localMax<scalar> maxInterp(mesh);
        surfaceScalarField Gamma
        (
            "PEQSI:avGamma",
            k2*maxInterp.interpolate(psiS)()*maxInterp.interpolate(lambda)()
           /mesh.nonOrthDeltaCoeffs()
        );
        // explicit-diffusion stability ceiling: dt Gamma delta^2 <= 1/4
        Gamma = min
        (
            Gamma,
            0.25/(dt*sqr(mesh.nonOrthDeltaCoeffs()))
        );

        const volScalarField rho0(rho_);
        rho_ += dt*fvc::laplacian(Gamma, rho0);
        rho_.correctBoundaryConditions();
        volScalarField rhoSafe("PEQSI:avRhoDen", rho_);
        {
            scalarField& rd = rhoSafe.primitiveFieldRef();
            forAll(rd, celli)
            {
                if (mag(rd[celli]) < 0.01)
                {
                    rd[celli] = rd[celli] < 0 ? -0.01 : 0.01;
                }
            }
            rhoSafe.correctBoundaryConditions();
        }

        auto smooth = [&](volScalarField& q)
        {
            const volScalarField rq(rho0*q);
            q = (rq + dt*fvc::laplacian(Gamma, rq))/rhoSafe;
            q.correctBoundaryConditions();
        };

        {
            const volVectorField rU(rho0*U_);
            U_ = (rU + dt*fvc::laplacian(Gamma, rU))/rhoSafe;
            U_.correctBoundaryConditions();
        }
        smooth(h_);
        if (fgmActive_)
        {
            smooth(Z_());
            smooth(Zvar_());
            smooth(Yc_());
        }

        if (pimple.dict().lookupOrDefault<Switch>("peqsiStiffCensus", false))
        {
            const label every =
                pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
            static label n = 0;
            if (every > 0 && (n++ % every) == 0)
            {
                Info<< "PEQSI AV: sensor max = " << gMax(psiS)
                    << ", Gamma max = " << gMax(Gamma.primitiveField())
                    << " m2/s" << endl;
            }
        }
    }

    hStarSaved_.reset(new volScalarField("PEQSI:hStarSaved", h_));

    // Algebraic subgrid mixture-fraction variance (peqsiZvarCoeff > 0):
    // Pierce & Moin (Phys. Fluids 10:3041, 1998), constant-coefficient
    // form Zvar = C_v Delta^2 |grad Z|^2 with Delta = V^(1/3), the
    // standard closure of flamelet/FGM LES.  The solver's "Zvar" field
    // is the SEGREGATION factor g = Zvar/(Z(1-Z)) -- the closure hands
    // it to the table's gZ axis unscaled -- so the model value is
    // normalised here and clamped to the axis range.  Cells of nearly
    // pure stream (Z(1-Z) < 1e-6) carry no meaningful variance and get
    // g = 0.  Overwrites the advected field: with the model armed the
    // transport is a placeholder (revisit if a transported-variance
    // model is ever wanted).  Default 0 keeps the old behaviour
    // bit-identically.
    const scalar zvarCv
    (
        pimple.dict().lookupOrDefault<scalar>("peqsiZvarCoeff", 0.0)
    );
    if (fgmActive_ && zvarCv > 0)
    {
        const tmp<volVectorField> tgZ(fvc::grad(Z_()));
        const vectorField& gZ = tgZ().primitiveField();
        const scalarField& Zf = Z_().primitiveField();
        const scalarField& Vf = mesh.V();
        scalarField& gv = Zvar_().primitiveFieldRef();
        // Track the RAW model output as well as the stored (clamped)
        // one.  g = Zvar/(Z(1-Z)) is a segregation factor and cannot
        // exceed 1 for any realisable state, so a raw value above the
        // clamp is the model reporting that its input |grad Z| is not
        // the gradient of a filtered field.  The stored field alone
        // cannot show this -- it saturates at 0.99 and looks healthy.
        scalar gRawMax = 0;
        label nOver = 0;
        forAll(gv, i)
        {
            const scalar zz = Zf[i]*(1.0 - Zf[i]);
            if (zz > 1e-6)
            {
                const scalar d2 = cbrt(sqr(Vf[i]));
                const scalar raw = zvarCv*d2*magSqr(gZ[i])/zz;
                gRawMax = max(gRawMax, raw);
                if (raw > 0.99) nOver++;
                gv[i] = min(raw, scalar(0.99));
            }
            else
            {
                gv[i] = 0.0;
            }
        }
        Zvar_().correctBoundaryConditions();

        if (pimple.dict().lookupOrDefault<Switch>("peqsiStiffCensus", false))
        {
            const label every =
                pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
            static label n = 0;
            if (every > 0 && (n++ % every) == 0)
            {
                reduce(gRawMax, maxOp<scalar>());
                reduce(nOver, sumOp<label>());
                Info<< "PEQSI Zvar (Pierce-Moin): g max = "
                    << gMax(Zvar_()) << ", mean = "
                    << gAverage(Zvar_())
                    << ", raw max = " << gRawMax
                    << ", above clamp on " << nOver << " cells" << endl;
            }
        }
    }

    // h-budget probe (peqsiStiffCensus): who is heating the field?  The
    // rd0110 restart shows a LINEAR global T rise (~1.7 K/step over 88%
    // of cells) that survives outlet-BC and LES/laminar bisection; the
    // remaining suspects are this substep's h equation and the acoustic
    // Eq. 24 update.  Print the substep's own d(int rho h) so the two
    // can be told apart.
    if (pimple.dict().lookupOrDefault<Switch>("peqsiStiffCensus", false))
    {
        const scalar rhAfter =
            gSum(rho_.primitiveField()*h_.primitiveField()
                *mesh.V().primitiveField());
        const scalar rhBefore =
            gSum(rhoN_().primitiveField()*hN_().primitiveField()
                *mesh.V().primitiveField());
        Info<< "PEQSI h-budget: advective substep d(int rho h) = "
            << rhAfter - rhBefore << " J" << endl;

        // The same substep measured against the manifold rather than
        // against zero.  d(int rho h) answers "is the field heating",
        // which a uniform compression also does; the defect
        // int rho (h - h_mix(Z)) answers "is h leaving Z behind",
        // which is the thing the manifold cannot tolerate.  Printing
        // it here in the same units the acoustic ledger uses makes the
        // two substeps' contributions addable, and that sum is what
        // decides which of them creates the defect and which absorbs
        // it -- a distinction we cannot make from either number alone.
        if (fgmActive_ && fgmTable_.valid())
        {
            const scalar hOx = fgmTable_().hOx();
            const scalar hFu = fgmTable_().hFuel();
            const scalarField& V = mesh.V();
            const scalarField& Za = Z_().primitiveField();
            const scalarField& Zb = ZN_().primitiveField();
            const scalarField& ra = rho_.primitiveField();
            const scalarField& rb = rhoN_().primitiveField();
            const scalarField& ha = h_.primitiveField();
            const scalarField& hb = hN_().primitiveField();

            // Unclamped on purpose.  The manifold lookup clamps Z to
            // [0, 1] before forming its mixing line, so a cell sitting
            // at 1 + eps has h_mix evaluated at 1 and books
            // eps*(h_fu - h_ox) of defect every step -- 1.6 J/kg per
            // 1e-6 of overshoot, from neither scheme nor diffusion.
            // Using the clamped value here would hide that in the
            // budget instead of showing it, so the excursion is
            // counted separately below.
            scalar dAfter = 0, dBefore = 0, zOver = 0, zOverE = 0;
            forAll(V, i)
            {
                dAfter +=
                    ra[i]*(ha[i] - ((1.0 - Za[i])*hOx + Za[i]*hFu))*V[i];
                dBefore +=
                    rb[i]*(hb[i] - ((1.0 - Zb[i])*hOx + Zb[i]*hFu))*V[i];

                const scalar zc = min(max(Za[i], 0.0), 1.0);
                if (zc != Za[i])
                {
                    zOver += mag(Za[i] - zc)*V[i];
                    zOverE += ra[i]*(Za[i] - zc)*(hFu - hOx)*V[i];
                }
            }
            reduce(zOver, sumOp<scalar>());
            reduce(zOverE, sumOp<scalar>());
            reduce(dAfter, sumOp<scalar>());
            reduce(dBefore, sumOp<scalar>());

            Info<< "PEQSI defect budget: advective d(int rho dh) = "
                << dAfter - dBefore << " J (level " << dAfter
                << " J), Z outside [0,1] on " << zOver
                << " m3 worth " << zOverE << " J" << endl;

            // The same accounting restricted to the population the
            // reachability census actually counts.  int rho dh is
            // 98.7% chamber gas sitting above the cold LOX anchor --
            // correct, expected, and six orders larger than the
            // interface excess, so nothing done at the interface can
            // move it.  Excess above the upper mixing line is exactly
            // zero wherever the state is reachable, so what is left is
            // the violation itself.
            if (hChamber_ == GREAT)
            {
                // The WARMEST oxidiser-side state, not their mean.
                // Both the chamber gas and the LOX sit at Z = 0, so a
                // mass-weighted mean over Z < 1e-6 lands between them
                // and puts the chamber itself above its own envelope --
                // every cell in the domain then counts as a violation,
                // which is the same population error one level down.
                // The upper edge is defined by mixing fuel with the
                // hottest oxidiser present, so take the maximum.
                scalar hMaxOx = -GREAT, hMaxFu = -GREAT;
                forAll(V, i)
                {
                    if (Zb[i] < 1e-6) hMaxOx = max(hMaxOx, hb[i]);
                    if (Zb[i] > 1.0 - 1e-6) hMaxFu = max(hMaxFu, hb[i]);
                }
                reduce(hMaxOx, maxOp<scalar>());
                reduce(hMaxFu, maxOp<scalar>());
                hChamber_ = hMaxOx > -GREAT ? hMaxOx : hOx;

                // Both endpoints from the field, for the same reason.
                // Taking one from the field and the other from the
                // table's constant is what put the whole fuel stream
                // outside its own envelope.
                hFuelRef_ = hMaxFu > -GREAT ? hMaxFu : hFu;
            }

            scalar eA = 0, eB = 0;
            label nE = 0;
            forAll(V, i)
            {
                const scalar za = min(max(Za[i], 0.0), 1.0);
                const scalar zb = min(max(Zb[i], 0.0), 1.0);
                const scalar ea =
                    ha[i] - ((1.0 - za)*hChamber_ + za*hFuelRef_);
                const scalar eb =
                    hb[i] - ((1.0 - zb)*hChamber_ + zb*hFuelRef_);
                if (ea > 0) { eA += ra[i]*ea*V[i]; nE++; }
                if (eb > 0) eB += rb[i]*eb*V[i];
            }
            reduce(eA, sumOp<scalar>());
            reduce(eB, sumOp<scalar>());
            reduce(nE, sumOp<label>());

            if (envBooked_ == 0 && envExcessPrev_ < 0) envE0_ = eB;
            envBooked_ += eA - eB;

            Info<< "PEQSI envelope budget: advective dE = " << eA - eB
                << " J (E = " << eA << " J on " << nE
                << " cells, ends " << hChamber_ << " / " << hFuelRef_
                << " J/kg frozen, booked " << envBooked_
                << " vs E-E0 " << eA - envE0_ << " J)" << endl;

            envExcessPrev_ = eA;
        }
    }

    mark(tPhase_[3]);
}




void Foam::solvers::peqsiFluid::advectManifoldStage
(
    const surfaceScalarField& phiv,
    const volScalarField& divPhiv,
    const scalar cOld,
    const scalar cNew,
    const scalar dtv
)
{
    // One SSP-RK3 stage of the manifold coordinates, ADVECTIVE form:
    //
    //     dq/dt + u . grad q = S,      q in {Z, Zvar, Yc}
    //
    // Non-conservative on purpose, twice over.  First, composition must
    // not respond to the acoustic projection at a contact -- the dp
    // kick changes rho, and a conservative rho*q carrier would need the
    // same Delta-M pairing the momentum needs; the advective form makes
    // the invariance exact by construction (TK 2012 transport their
    // mass fractions the same way).  Second, u . grad q obeys a maximum
    // principle, so with a bounded reconstruction the coordinates stay
    // in their physical ranges without clipping.
    //
    // The price is that int(rho q) is conserved only through rho's own
    // ledger -- acceptable for the stage-2a wiring scope, revisited
    // when the source terms arrive.
    //
    // Boundedness hygiene: the TK Eq. 27 bounding artificial
    // diffusivity on Z (C_Y = 100, active only OUTSIDE [0,1], so it is
    // the identity on a healthy field).
    auto stageOne =
        [&](volScalarField& q, const volScalarField& qN, const word& key,
            const volScalarField* Sq)
    {
        tmp<volScalarField> tLq =
            -(fvc::div(phiv, q, key) - q*divPhiv);

        if (Sq)
        {
            tLq.ref() += *Sq/max(rho_, dimensionedScalar(dimDensity, small));
        }

        const scalarField& Lq = tLq().primitiveField();
        const scalarField& qNf = qN.primitiveField();
        scalarField& qf = q.primitiveFieldRef();

        forAll(qf, i)
        {
            qf[i] = cOld*qNf[i] + cNew*(qf[i] + dtv*Lq[i]);
        }
        q.correctBoundaryConditions();
    };

    // Packed carrier (peqsiPackManifold, default on): (Z, Zvar, Yc) as
    // one vector field so the WENO reconstruction loads each cell's
    // stencil matrices ONCE per stage instead of three times -- the same
    // memory-bandwidth argument, and the same per-component weight
    // computation, as the symmTensor carrier for the conservative set
    // (WENOCoeff loops nComp inside the stencil loop; the math per
    // component is identical to the per-field passes).  One scheme
    // entry, "div(phiv,Zpack)", serves all three components; the
    // bounded 01 fit is valid for each: Z in [0,1], Zvar in [0, 1/4],
    // Yc in [0, Cnorm] -- all inside [0, 1].
    if (pimple.dict().lookupOrDefault<Switch>("peqsiPackManifold", true))
    {
        wordList qBcTypes(mesh.boundary().size());
        forAll(mesh.boundary(), patchi)
        {
            const word& pt = mesh.boundaryMesh()[patchi].type();
            qBcTypes[patchi] =
                (pt == "patch" || pt == "wall")
              ? word("calculated")
              : Z_().boundaryField()[patchi].type();
        }
        volVectorField Q3
        (
            IOobject("PEQSI:Zpack", mesh.time().name(), mesh),
            mesh,
            dimensionedVector(dimless, Zero),
            qBcTypes
        );
        {
            vectorField& Qf = Q3.primitiveFieldRef();
            const scalarField& Zf = Z_().primitiveField();
            const scalarField& Zvf = Zvar_().primitiveField();
            const scalarField& Ycf = Yc_().primitiveField();
            forAll(Qf, i)
            {
                Qf[i] = vector(Zf[i], Zvf[i], Ycf[i]);
            }
            forAll(Q3.boundaryFieldRef(), patchi)
            {
                fvPatchVectorField& Qp = Q3.boundaryFieldRef()[patchi];
                if (Qp.coupled()) continue;
                const fvPatchScalarField& Zp =
                    Z_().boundaryField()[patchi];
                const fvPatchScalarField& Zvp =
                    Zvar_().boundaryField()[patchi];
                const fvPatchScalarField& Ycp =
                    Yc_().boundaryField()[patchi];
                forAll(Qp, i)
                {
                    Qp[i] = vector(Zp[i], Zvp[i], Ycp[i]);
                }
            }
            Q3.correctBoundaryConditions();
        }

        // Operand check (peqsiPackCheck): the vector-WENO FPE names
        // calcCoeff but not its input -- scan the carrier for
        // non-finite components BEFORE the reconstruction and dump the
        // first offender's cell state, so "bad input" and "bad scheme
        // math" are told apart without burning a dt ladder (the
        // Cantera high-chi chain precedent: two plausible physical
        // stories, the real cause only visible in the operand).
        if (pimple.dict().lookupOrDefault<Switch>("peqsiPackCheck", false))
        {
            const vectorField& Qf = Q3.primitiveField();
            label bad = -1;
            forAll(Qf, i)
            {
                const vector& q = Qf[i];
                if
                (
                    !std::isfinite(q.x())
                 || !std::isfinite(q.y())
                 || !std::isfinite(q.z())
                )
                {
                    bad = i;
                    break;
                }
            }
            if (bad >= 0)
            {
                Pout<< "PEQSI pack check: NON-FINITE carrier at "
                    << mesh.C()[bad]
                    << " Q3 = " << Qf[bad]
                    << " rho = " << rho_[bad]
                    << " h = " << h_[bad]
                    << " p = " << p_[bad]
                    << " T = " << thermo_.T()[bad] << endl;
            }
        }

        volVectorField LQ3
        (
            -(fvc::div(phiv, Q3, "div(phiv,Zpack)") - Q3*divPhiv)
        );

        // Subgrid scalar flux (peqsiScSGS > 0 to enable).
        //
        // The energy equation adds muSgs*Cp/PrSgs to kappa, but the
        // mixture fraction had no subgrid term at all: h was closed by
        // an explicit SGS model and Z by the scheme's implicit
        // dissipation alone.  Measured on rd0110: muSgs reaches 2x the
        // molecular viscosity at the median chamber cell and 90x in the
        // tail, so the two are not closed alike, and the manifold
        // assumes they move together.  What that asymmetry produced is
        // documented: the resolved |grad Z| is too steep, the algebraic
        // segregation factor exceeds its own realizability bound
        // (measured g = 2.02 against g <= 1), cells leave the manifold
        // on the specific-volume gate, the EOS fallback is handed an
        // (h, v) pair no temperature satisfies, and the floor clamp
        // fills with cells until the run collapses (t = 404 us).
        //
        // Standard gradient-diffusion closure, D_sgs = nu_t/Sc_sgs, in
        // the advective form the substep uses: (1/rho) div(mu_sgs/Sc
        // grad q).  Carries the same explicit-stability cap as the
        // bounding diffusivity below -- an unstable subgrid term would
        // amplify what it exists to smooth.
        const scalar scSgs =
            pimple.dict().lookupOrDefault<scalar>("peqsiScSGS", 0.7);

        // Molecular species diffusion (peqsiLe > 0 to enable).
        //
        // Takahashi is off at Tier-0, so the mixture fraction carried
        // no molecular diffusion either -- only the subgrid term above
        // and the scheme's own dissipation, while the enthalpy has
        // always had the full thermal conductivity.  That asymmetry is
        // the larger of the two: at Le = 0.29 the molecular
        // diffusivity is 2.1x the subgrid one in the fuel and 1.5x in
        // the chamber, which is why closing only the subgrid half
        // slowed the reachability violations without stopping them.
        //
        // rho*D = kappa/(cp Le) with the MOLECULAR conductivity, not
        // the one the energy equation has already had muSgs folded
        // into -- otherwise the subgrid part is counted twice.
        //
        // Le = 1 is not an assumption here, it is the manifold's own
        // convention worked through: the table defines the mixture
        // fraction's Lewis number as Le_Z = Pr = nu/alpha, so the
        // diffusivity it implies is D_Z = nu/Le_Z = alpha exactly, and
        // the viscosity cancels.  That cancellation is worth having --
        // high-pressure Chung overpredicts mu by ~300x in the cold
        // fuel, and none of that error reaches Z through this term.
        // The progress variable is a different matter (Le_C uses a
        // Takahashi-corrected D and does not cancel), but Yc is zero
        // throughout a non-reacting run.
        const scalar leZ =
            pimple.dict().lookupOrDefault<scalar>("peqsiLe", 1.0);

        if ((sgsActive_ && scSgs > 0) || leZ > 0)
        {
            volScalarField Dsgs
            (
                "PEQSI:Dsgs",
                0.0*rhoN_()*dimensionedScalar(dimArea/dimTime, 1.0)
            );

            if (leZ > 0)
            {
                const volScalarField kappaMol
                (
                    completeField("PEQSI:kappaMol", mesh, thermo_.kappa())
                );
                Dsgs += kappaMol/(thermo_.Cp()*leZ);
            }

            if (sgsActive_ && scSgs > 0)
            {
                Dsgs += rhoN_()*momentumTransport->nut()/scSgs;
            }

            {
                ensureDirGeometry();
                const scalarField& dmin = ladDeltaMin_();
                const scalarField& rf0 = rho_.primitiveField();
                scalarField& Df = Dsgs.primitiveFieldRef();
                const scalar rdt = 0.45/runTime.deltaTValue();
                label nCap = 0;
                forAll(Df, i)
                {
                    const scalar cap = rdt*sqr(dmin[i])*rf0[i];
                    if (Df[i] > cap) { Df[i] = cap; nCap++; }
                }
                Dsgs.correctBoundaryConditions();

                // A capped subgrid diffusivity is a subgrid model that
                // is not being applied.  Without this count, "the SGS
                // flux did not help" and "the SGS flux never got in"
                // read the same in the log.
                if
                (
                    pimple.dict().lookupOrDefault<Switch>
                    (
                        "peqsiStiffCensus", false
                    )
                )
                {
                    reduce(nCap, sumOp<label>());
                    const label every =
                        pimple.dict().lookupOrDefault<label>
                        (
                            "peqsiDiagInterval", 10
                        );
                    static label nS = 0;
                    if (every > 0 && (nS++ % every) == 0 && nCap)
                    {
                        Info<< "PEQSI SGS scalar flux: diffusivity "
                            << "capped on " << nCap << " cells" << endl;
                    }
                }
            }

            const volScalarField rhoSafe
            (
                "PEQSI:rhoSgsDen",
                max(rho_, dimensionedScalar(dimDensity, 0.01))
            );

            const volScalarField dZ
            (
                "PEQSI:sgsZ", fvc::laplacian(Dsgs, Z_())/rhoSafe
            );
            const volScalarField dZv
            (
                "PEQSI:sgsZvar", fvc::laplacian(Dsgs, Zvar_())/rhoSafe
            );
            const volScalarField dYc
            (
                "PEQSI:sgsYc", fvc::laplacian(Dsgs, Yc_())/rhoSafe
            );

            vectorField& Lw = LQ3.primitiveFieldRef();
            const scalarField& a = dZ.primitiveField();
            const scalarField& b = dZv.primitiveField();
            const scalarField& c = dYc.primitiveField();
            forAll(Lw, i)
            {
                Lw[i] += vector(a[i], b[i], c[i]);
            }

            forAll(LQ3.boundaryFieldRef(), patchi)
            {
                fvPatchVectorField& Lb = LQ3.boundaryFieldRef()[patchi];
                const fvPatchScalarField& ab =
                    dZ.boundaryField()[patchi];
                const fvPatchScalarField& bb =
                    dZv.boundaryField()[patchi];
                const fvPatchScalarField& cb =
                    dYc.boundaryField()[patchi];
                forAll(Lb, i)
                {
                    Lb[i] += vector(ab[i], bb[i], cb[i]);
                }
            }
        }

        const vectorField& Lf = LQ3.primitiveField();
        const scalarField& Sf = sourceYc_().primitiveField();
        const scalarField& rf = rho_.primitiveField();
        const scalarField& ZNf = ZN_().primitiveField();
        const scalarField& ZvNf = ZvarN_().primitiveField();
        const scalarField& YcNf = YcN_().primitiveField();
        scalarField& Zf = Z_().primitiveFieldRef();
        scalarField& Zvf = Zvar_().primitiveFieldRef();
        scalarField& Ycf = Yc_().primitiveFieldRef();
        forAll(Zf, i)
        {
            Zf[i]  = cOld*ZNf[i]  + cNew*(Zf[i]  + dtv*Lf[i].x());
            Zvf[i] = cOld*ZvNf[i] + cNew*(Zvf[i] + dtv*Lf[i].y());
            const scalar LYc = Lf[i].z() + Sf[i]/max(rf[i], small);
            Ycf[i] = cOld*YcNf[i] + cNew*(Ycf[i] + dtv*LYc);
        }
        Z_().correctBoundaryConditions();
        Zvar_().correctBoundaryConditions();
        Yc_().correctBoundaryConditions();
    }
    else
    {
        stageOne(Z_(), ZN_(), "div(phiv,Z)", nullptr);
        stageOne(Zvar_(), ZvarN_(), "div(phiv,Zvar)", nullptr);
        stageOne(Yc_(), YcN_(), "div(phiv,Yc)", &sourceYc_());
    }

    // Realizability on Yc, not on c.  The table axis is c = Yc/Cnorm(Z),
    // and Cnorm falls to 1.59e-3 at Z = 0, so an absolute Yc error of
    // 1e-4 -- ordinary transport noise -- lands as dc = 0.063 there
    // against dc = 4e-4 at Z = 0.1.  Clamping c does the arithmetic in
    // that ill-conditioned variable and, worse, leaves the transported
    // Yc infeasible: measured 2026-08-23 on dp8lo, a pure-oxidiser cell
    // reached c = 4.31, i.e. Yc = 6.9e-3 of combustion products in a
    // cell holding no fuel, and the run died there.  Yc <= Cnorm(Z) is
    // the same set of states expressed in the well-conditioned
    // variable, and it leaves the c axis definition untouched so every
    // existing table stays valid.
    if
    (
        pimple.dict().lookupOrDefault<Switch>("peqsiBoundYc", true)
     && fgmTable_.valid() && fgmTable_().hasCnorm()
    )
    {
        const scalarField& Zf = Z_().primitiveField();
        scalarField& Ycf = Yc_().primitiveFieldRef();
        label nLo = 0, nHi = 0;
        scalar worst = 0;

        forAll(Ycf, i)
        {
            const scalar Cn =
                fgmTable_().interpolateCnorm(min(max(Zf[i], 0.0), 1.0));

            if (Ycf[i] < 0) nLo++;
            if (Ycf[i] > Cn)
            {
                nHi++;
                worst = max(worst, Ycf[i]/max(Cn, small));
            }
            Ycf[i] = min(max(Ycf[i], 0.0), Cn);
        }
        Yc_().correctBoundaryConditions();

        // Once this clamp is on, the c census can no longer report c > 1
        // -- it is the same bound seen from the other side -- so the
        // violation count has to be reported from here or it goes dark.
        if (pimple.dict().lookupOrDefault<Switch>("peqsiStiffCensus", false))
        {
            reduce(nLo, sumOp<label>());
            reduce(nHi, sumOp<label>());
            reduce(worst, maxOp<scalar>());

            const label every =
                pimple.dict().lookupOrDefault<label>("peqsiDiagInterval", 10);
            static label nYc = 0;
            if (every > 0 && (nYc++ % every) == 0 && (nLo || nHi))
            {
                Info<< "PEQSI Yc realizability: below 0 on " << nLo
                    << ", above Cnorm(Z) on " << nHi
                    << " cells; worst Yc/Cnorm = " << worst << endl;
            }
        }
    }

    // Bounding diffusivity on Z (explicit, inside the stage): inactive
    // on a bounded field, O(dx) strength where an excursion exists
    if (pimple.dict().lookupOrDefault<Switch>("peqsiBoundZ", true))
    {
        // Early exit on a bounded field.  D* is identically zero then, so
        // the whole block below evaluates to Z += lap(0, Z) == Z exactly
        // -- at the cost of ~6 field constructions and their halo
        // exchanges, PER RK STAGE.  One local scan + one reduction buys
        // all of that back on every healthy step (same pattern as the
        // SC filter's quiet-direction skip; bit-identical by the
        // zero-flux argument above).
        scalar violMax = 0;
        {
            const scalarField& Zf = Z_().primitiveField();
            forAll(Zf, i)
            {
                violMax = max(violMax, max(Zf[i] - 1.0, -Zf[i]));
            }
            reduce(violMax, maxOp<scalar>());
        }
        if (violMax > 0)
        {
        // speed of sound from the App-D identity beta/(1-alpha) = -rho c^2
        // (alpha_/beta_ are refreshed every step by updateCoefficients)
        const volScalarField cSound
        (
            "PEQSI:cSound",
            sqrt
            (
                max
                (
                    -beta_/((1.0 - alpha_)
                   *max(rho_, dimensionedScalar(dimDensity, small))),
                    dimensionedScalar(sqr(dimVelocity), small)
                )
            )
        );
        tmp<volScalarField> tD =
            boundingArtDiffusivity(Z_(), cSound, 100.0);

        // Explicit-stability cap, the SAME 0.45 diffusion-number guard
        // the LAD mass diffusivity carries -- this term never got one.
        // Uncapped, the demand D* = 100 c viol Delta crosses the
        // explicit limit 0.45 dmin^2/dt at a Z overshoot of only ~15%
        // on the 2-D case's grid, and an explicitly unstable bounding
        // term AMPLIFIES the violation it exists to remove.  Inactive
        // on healthy fields (bit-identical there), it only disarms the
        // pathological cell, exactly like the Dart cap.
        {
            ensureDirGeometry();
            const scalarField& dmin = ladDeltaMin_();
            scalarField& D = tD.ref().primitiveFieldRef();
            const scalar rdt = 0.45/runTime.deltaTValue();
            forAll(D, i)
            {
                D[i] = min(D[i], rdt*sqr(dmin[i]));
            }
        }

        Z_() += runTime.deltaT()*cNew*fvc::laplacian(tD(), Z_());
        Z_().correctBoundaryConditions();
        }
    }
}


// ************************************************************************* //
