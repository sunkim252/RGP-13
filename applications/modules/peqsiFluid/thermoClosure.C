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

    reduce(maxDev, maxOp<scalar>());
    Info<< "PEQSI coefficients: max |beta/(1-alpha) + rho c^2| / (rho c^2) = "
        << maxDev << endl;

    alpha_.correctBoundaryConditions();
    beta_.correctBoundaryConditions();
}


void Foam::solvers::peqsiFluid::invertTemperature()
{
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
            // Constant-v closure, driven through the thermo's own p.
            //
            // No new thermo API is needed: p_ IS thermo_.p(), so writing
            // a trial pressure into it and calling correct() evaluates
            // the equation of state exactly at that pressure.  That
            // turns the p-direction of the Newton from a linearisation
            // (which failed: the step needed is tens of MPa, far outside
            // the range where rho is linear in p) into a real one, and
            // it touches no shared thermophysical code.
            //
            // On exit the properties are left at the CONSISTENT state
            // (pEos, T) while p_ is restored to the transported
            // pressure the momentum and pressure equations need.  That
            // is the point: psi, cp, cv -- and therefore alpha and beta
            // -- stop being a mix of the transported density with
            // properties taken at a pressure that does not produce it.
            // ----------------------------------------------------------
            const label nOuter =
                pimple.dict().lookupOrDefault<label>("peqsiConstantVIters", 8);

            const volScalarField pSave("PEQSI:pSave", p_);
            const scalarField& rhoTf = rho_.primitiveField();
            scalarField& pf = p_.primitiveFieldRef();

            const scalar dpFrac = 0.5;   // per-pass |dp|/p clamp

            for (; iter < nOuter && maxRel > tol; ++iter)
            {
                maxRel = 0;

                thermo_.correct();

                const auto tCp = thermo_.Cp();
                const auto tCv = thermo_.Cv();
                const scalarField& cpf = tCp().primitiveField();
                const scalarField& cvf = tCv().primitiveField();
                const scalarField& psif = thermo_.psi().primitiveField();
                const scalarField& rhoEf = thermo_.rho()().primitiveField();

                const volScalarField hk(thermo_.he(p_, Tw));
                const scalarField& hkf = hk.primitiveField();

                forAll(Tf, i)
                {
                    const scalar rho = max(rhoEf[i], small);
                    const scalar v = 1.0/rho;
                    const scalar psi = max(psif[i], small);
                    const scalar dpdv = -sqr(rho)/psi;
                    const scalar dpdT =
                        rho*sqrt(max((cpf[i] - cvf[i])/(Tf[i]*psi), 0.0));
                    const scalar dvdT = -dpdT/dpdv;      // (dv/dT)_p

                    //   [ psi          -rho^2 (dv/dT)_p ] [dp]   [F1]
                    //   [ v - T(dv/dT)_p       cp       ] [dT] = [F2]
                    const scalar J11 = psi;
                    const scalar J12 = -sqr(rho)*dvdT;
                    const scalar J21 = v - Tf[i]*dvdT;
                    const scalar J22 = cpf[i];
                    const scalar det = J11*J22 - J12*J21;

                    const scalar F1 = rhoTf[i] - rhoEf[i];
                    const scalar F2 = hf[i] - hkf[i];

                    scalar dp, dT;
                    if (mag(det) > vSmall)
                    {
                        dp = (F1*J22 - F2*J12)/det;
                        dT = (J11*F2 - J21*F1)/det;
                    }
                    else
                    {
                        dp = 0;
                        dT = F2/max(cpf[i], small);
                    }

                    // Clamps: the first passes can ask for enormous
                    // steps out of a badly inconsistent start, and an
                    // unclamped step leaves the SRK validity range (and
                    // can drive p negative, where the root solver's
                    // floor silently takes over)
                    const scalar dpLim = dpFrac*mag(pf[i]);
                    dp = min(max(dp, -dpLim), dpLim);
                    dT = min(max(dT, -dTmax), dTmax);

                    pf[i] = max(pf[i] + dp, 1e4);

                    const scalar Tnew = Tf[i] + dT;
                    if ((Tnew <= Tmin && dT < 0) || (Tnew >= Tmax && dT > 0))
                    {
                        Tf[i] = Tnew <= Tmin ? Tmin : Tmax;
                        nSaturated_++;
                        continue;
                    }
                    Tf[i] = min(max(Tnew, Tmin), Tmax);

                    maxRel =
                        max
                        (
                            maxRel,
                            max
                            (
                                mag(F1)/max(rhoTf[i], small),
                                mag(dT)/max(Tf[i], small)
                            )
                        );
                }

                reduce(maxRel, maxOp<scalar>());
                p_.correctBoundaryConditions();
                Tw.correctBoundaryConditions();
            }

            // Refresh properties at the converged, CONSISTENT state...
            thermo_.correct();

            // ...then hand the transported pressure back to the solver.
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
                << " Pa, outer iters = " << iter
                << ", maxRel = " << maxRel << endl;
        }
        else
        {
            // cp slope evaluated ONCE per step (adequate: cp varies
            // slowly over the per-step temperature increments; profiling
            // showed the per-iteration Cp() re-evaluation doubled the
            // Newton cost)
            const auto tCp = thermo_.Cp();
            const scalarField& cpf = tCp().primitiveField();

            for (; iter < 60 && maxRel > tol; ++iter)
            {
                maxRel = 0;

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
                    maxRel = max(maxRel, mag(dT)/max(Tf[i], small));
                }

                reduce(maxRel, maxOp<scalar>());
            }
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
        thermo_.correct();
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
