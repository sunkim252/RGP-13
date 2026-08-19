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
    const auto tGamma = thermo_.gamma();

    const scalarField& cp = tCp().primitiveField();
    const scalarField& cv = tCv().primitiveField();
    const scalarField& gam = tGamma().primitiveField();
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

        // App. D identity: beta/(1-alpha) == -rho c^2 = -gamma/psi
        const scalar lhs = b[i]/(1.0 - a[i]);
        const scalar rhs = -gam[i]/psi[i]/rho[i]*rho[i];  // -gamma/psi
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

    const volScalarField rhoTrans("PEQSI:rhoTrans", rho_);

    thermo_.he() = h_;
    thermo_.he().correctBoundaryConditions();
    thermo_.correct();

    // rho_ now holds the EOS density: measure the drift, then restore
    scalar maxDrift = 0;
    {
        const scalarField& re = rho_.primitiveField();
        const scalarField& rt = rhoTrans.primitiveField();
        forAll(re, i)
        {
            maxDrift = max(maxDrift, mag(re[i] - rt[i])/max(rt[i], small));
        }
        reduce(maxDrift, maxOp<scalar>());
    }

    rho_ = rhoTrans;
    rho_.correctBoundaryConditions();

    Info<< "PEQSI thermo closure: T = ["
        << gMin(thermo_.T().primitiveField()) << ", "
        << gMax(thermo_.T().primitiveField())
        << "] K, rho drift (EOS vs transported) max = "
        << maxDrift << endl;
}


// ************************************************************************* //
