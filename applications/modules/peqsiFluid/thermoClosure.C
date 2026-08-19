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
    // Deliberate deviation from WKK Fig. 3 (documented in the wiki):
    // the reference inverts h(T, v) at fixed specific volume with the
    // slope xi; this implementation inverts at fixed p with the thermo's
    // own cp-slope Newton.  The difference is of EOS-residual order;
    // revisit against the 1-D validation if needed.
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

        // cp slope evaluated ONCE per step (adequate: cp varies slowly
        // over the per-step temperature increments; profiling showed the
        // per-iteration Cp() re-evaluation doubled the Newton cost)
        const auto tCp = thermo_.Cp();
        const scalarField& cpf = tCp().primitiveField();

        scalarField& Tf = Tw.primitiveFieldRef();
        const scalarField& hf = h_.primitiveField();

        // After the first full-field pass only the cells still moving
        // (steep-front cells, typically a small fraction) are iterated,
        // via the thermo's cell-subset he(T, cells) evaluation --
        // profiling showed the full-field SRK sweep per iteration
        // dominating the closure cost.
        labelList active;
        label nSaturated_ = 0;

        for (; iter < 60 && maxRel > tol; ++iter)
        {
            maxRel = 0;

            if (iter == 0)
            {
                const volScalarField hk(thermo_.he(p_, Tw));
                const scalarField& hkf = hk.primitiveField();

                DynamicList<label> nextActive(Tf.size()/8);
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
                    if (rel > tol) nextActive.append(i);
                    maxRel = max(maxRel, rel);
                }
                active.transfer(nextActive);
            }
            else
            {
                scalarField Tsub(active.size());
                forAll(active, k) Tsub[k] = Tf[active[k]];

                const tmp<scalarField> thk(thermo_.he(Tsub, active));
                const scalarField& hkf = thk();

                DynamicList<label> nextActive(active.size());
                forAll(active, k)
                {
                    const label i = active[k];
                    scalar dT = (hf[i] - hkf[k])/max(cpf[i], small);
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
                    if (rel > tol) nextActive.append(i);
                    maxRel = max(maxRel, rel);
                }
                active.transfer(nextActive);
            }

            reduce(maxRel, maxOp<scalar>());
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
    thermo_.correct();

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

    scalar maxDrift = 0;
    {
        const scalarField& re = thermo_.rho()().primitiveField();
        const scalarField& rt = rho_.primitiveField();
        forAll(rt, i)
        {
            maxDrift = max(maxDrift, mag(re[i] - rt[i])/max(rt[i], small));
        }
        reduce(maxDrift, maxOp<scalar>());
    }

    Info<< "PEQSI thermo closure: T = ["
        << gMin(thermo_.T().primitiveField()) << ", "
        << gMax(thermo_.T().primitiveField())
        << "] K, rho drift (EOS vs transported) max = "
        << maxDrift << endl;
}


// ************************************************************************* //
