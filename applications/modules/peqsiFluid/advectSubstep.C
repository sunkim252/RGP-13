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
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

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

    const dimensionedScalar& dt = runTime.deltaT();

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

    // IMEX: treat the Eq.22/24 artificial MASS diffusion implicitly in a
    // split step after the RK3 stages.  The explicit form caps dt at
    // 0.45 dmin^2/D (measured 3.0e-7 s on the 3-D case, 3.3x below the
    // Courant/maxDeltaT bound), and that cap cannot be relaxed by
    // lowering C_rho: the interface sharpens until |d4rho| ~ 1/(n dx)^4
    // restores the same D (measured C 2.5x down -> D 1.26x UP).
    const Switch ladImplicit
    (
        pimple.dict().lookupOrDefault<Switch>("peqsiLadImplicit", false)
    );

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
            forAll(tDart.ref().boundaryFieldRef(), patchi)
            {
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

                if (ladImplicit)
                {
                    // Implicit: no explicit-stability constraint, so dt is
                    // left to the Courant bound and the dt-DEPENDENT cap
                    // (0.45 dmin^2/dt) must not be used -- it scales as
                    // 1/dt, so at the larger dt it would clip normal LAD
                    // operation.  The pathological-cell guard is kept, but
                    // expressed on a dt-INDEPENDENT, grid-aware scale:
                    // an acoustic cell diffusivity c*dmin.  (An absolute
                    // ceiling such as peqsiLadDtDemandMax cannot serve --
                    // it is a fixed number, so on a 1-cm 1-D mesh where
                    // normal operation is O(1) m2/s it would clip by 500x,
                    // while on the 45-um 3-D mesh it never binds.)
                    const scalar fCap =
                        pimple.dict().lookupOrDefault<scalar>
                        (
                            "peqsiLadImplicitCap", 1.0
                        );
                    const scalarField& ci = c.primitiveField();
                    forAll(D, i)
                    {
                        D[i] = min(D[i], fCap*ci[i]*deltaMin[i]);
                    }
                    ladDtLimit_ = great;
                }
                else
                {
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
            }

            Info<< "PEQSI LAD: rho_art max = "
                << gMax(tDart().primitiveField()) << " m^2/s" << endl;
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
            forAll(tKappaArt.ref().boundaryFieldRef(), patchi)
            {
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

            Info<< "PEQSI LAD: kappa_art max = "
                << gMax(tKappaArt().primitiveField()) << " W/(m K)" << endl;
        }
    }

    // L_h: enthalpy-equation RHS without Dp/Dt (WKK Eq. 3) -- conduction
    // with the TK Eq. (33) artificial conductivity when active.
    const volScalarField Lh
    (
        "PEQSI:Lh",
        tKappaArt.valid()
      ? fvc::laplacian((kappa + tKappaArt())(), thermo_.T())()
      : fvc::laplacian(kappa, thermo_.T())()
    );

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
        if (tDart.valid() && !ladImplicit)
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
                if (sYi) Lp += iLi[i]*(*sYi)[i];
                scalar Lrh = l.xz() + iLi[i]*Lhi[i];
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
                if (sYb) Lp += iLb[i]*(*sYb)[i];
                scalar Lrh = l.xz() + iLb[i]*Lhb[i];
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

    // ---- IMEX split step: implicit Eq. 22/24 artificial mass diffusion
    //
    //     (r - r*)/dt = div(D grad r),    solved implicitly for r
    //     ru += u^n (r - r*)              <- Eq. 24, SAME increment
    //     (rh - rh*)/dt = div(D grad rh)
    //
    // Placed AFTER sComp_ on purpose: sComp is the RK3 quadrature of the
    // COMPRESSION term rho div(phiv), and the acoustic substep relies on
    // it to cancel the Eq. 22 density increment exactly.  LAD is a pure
    // divergence with Dart == 0 on every boundary, so it moves no net
    // mass (verified: the substep mass ledger is unchanged) and must not
    // enter that quadrature.
    //
    // Eq. 24 uses u^n, exactly as the explicit branch does inside the
    // stages (UNi there, UN_() here): under velocity equilibrium the
    // momentum increment is then u times the mass increment, which is
    // the invariance Eq. 24 exists to preserve.
    if (tDart.valid() && ladImplicit)
    {
        const dimensionedScalar rdtD("rdt", dimless/dimTime, 1.0/dt.value());

        // theta weighting of the split diffusion.  Default backward
        // Euler.  theta = 0.5 (Crank-Nicolson) was tried and is WORSE
        // here, which locates the residual discrepancy: measured on the
        // Mayer-scale 1-D interface at dt 1e-6, L1(implicit - explicit)
        // went 5.33e-3 (theta 1) -> 7.70e-3 (theta 0.5) and the 10-90
        // width stayed 55 against the explicit branch's 53 either way.
        // The gap is therefore not the diffusion's time order but the
        // SPLITTING: the explicit branch applies LAD inside each RK
        // stage, so it interacts with the advection within the substep,
        // while this branch applies it once at the end.  Both converge
        // to the same solution -- the width difference is gone by
        // dt 1e-7 -- so the cost is a dt-truncation effect at the
        // operating point, not a change of model.
        const scalar theta =
            pimple.dict().lookupOrDefault<scalar>("peqsiLadTheta", 1.0);

        const volScalarField rStar("PEQSI:rStarLAD", r);

        fvScalarMatrix rEqn
        (
            fvm::Sp(rdtD, r) - theta*fvm::laplacian(tDart(), r)
         ==
            rdtD*rStar + (1.0 - theta)*fvc::laplacian(tDart(), rStar)
        );
        rEqn.solve(mesh.solution().solverDict("rhoLAD"));

        const volScalarField rhStar("PEQSI:rhStarLAD", rh);

        fvScalarMatrix rhEqn
        (
            fvm::Sp(rdtD, rh) - theta*fvm::laplacian(tDart(), rh)
         ==
            rdtD*rhStar + (1.0 - theta)*fvc::laplacian(tDart(), rhStar)
        );
        rhEqn.solve(mesh.solution().solverDict("rhoLAD"));

        // Eq. 24 momentum increment, from the SAME mass increment
        {
            const scalarField& rNew = r.primitiveField();
            const scalarField& rOld = rStar.primitiveField();
            const vectorField& UNi = UN_().primitiveField();
            vectorField& rui = ru.primitiveFieldRef();
            forAll(rui, i)
            {
                rui[i] += UNi[i]*(rNew[i] - rOld[i]);
            }
        }

        r.correctBoundaryConditions();
        rh.correctBoundaryConditions();
        ru.correctBoundaryConditions();
    }

    // Publish the starred state into the solver fields (BCs from the
    // registered fields' own conditions)
    rho_ = r;
    rho_.correctBoundaryConditions();

    U_ = ru/r;
    U_.correctBoundaryConditions();

    p_ = pw;
    p_.correctBoundaryConditions();

    h_ = rh/r;
    h_.correctBoundaryConditions();

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

    stageOne(Z_(), ZN_(), "div(phiv,Z)", nullptr);
    stageOne(Zvar_(), ZvarN_(), "div(phiv,Zvar)", nullptr);
    stageOne(Yc_(), YcN_(), "div(phiv,Yc)", &sourceYc_());

    // Bounding diffusivity on Z (explicit, inside the stage): inactive
    // on a bounded field, O(dx) strength where an excursion exists
    if (pimple.dict().lookupOrDefault<Switch>("peqsiBoundZ", true))
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
        const tmp<volScalarField> tD =
            boundingArtDiffusivity(Z_(), cSound, 100.0);
        Z_() += runTime.deltaT()*cNew*fvc::laplacian(tD(), Z_());
        Z_().correctBoundaryConditions();
    }
}


// ************************************************************************* //
