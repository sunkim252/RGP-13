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
    // TODO(M1): exact SRK partials.  The faithful coefficients are
    //   xi    = (dh/dT)_v
    //   alpha = (1/(rho xi)) (dp/dT)_v                    (WKK Eq. 11)
    //   beta  = (1/rho)[(dp/dv)_T - (1/xi)(dp/dT)_v (dh/dv)_T]  (Eq. 12)
    // evaluated analytically from the SRK EoS (WKK App. B), with the
    // cell-wise consistency check  beta/(1-alpha) == -rho c^2 (App. D).
    //
    // M0 surrogate: the ideal-gas limit of alpha ((gamma-1)/gamma) and
    // the App. D identity for beta with the isentropic sound speed
    // c^2 = gamma/psi, psi = (drho/dp)_T (real for SRKGas).  Exact in
    // the ideal-gas limit; a controlled approximation for the scaffold.
    // ------------------------------------------------------------------

    const volScalarField gamma(thermo_.gamma());
    const volScalarField c2("PEQSI:c2", gamma/thermo_.psi());

    alpha_ = (gamma - 1.0)/gamma;
    alpha_.correctBoundaryConditions();

    beta_ = -rho_*c2*(1.0 - alpha_);
    beta_.correctBoundaryConditions();
}


void Foam::solvers::peqsiFluid::invertTemperature()
{
    // ------------------------------------------------------------------
    // TODO(M1): Newton inversion of h(T, v, Y) at fixed specific volume
    // v = 1/rho for T, with the slope xi = (dh/dT)_v (WKK Fig. 3 --
    // explicitly NOT cp), initial guess = current T (or the manifold T
    // in the reacting extension), tolerance 1e-6 as in WKK.  Then
    // refresh mu/kappa (Chung) at the new (T, rho, Y) state WITHOUT
    // calling thermo_.correct(), which would overwrite the transported
    // density with the EOS density and destroy the scheme's mass
    // conservation.
    //
    // M0 scaffold: T is left lagged (properties frozen at the previous
    // state).  The 1-D validation (M3) requires the real inversion.
    // ------------------------------------------------------------------

    static bool warned = false;
    if (!warned)
    {
        warned = true;
        WarningInFunction
            << "M0 scaffold: T inversion not yet implemented -- "
            << "temperature and transport properties are lagged"
            << endl;
    }
}


// ************************************************************************* //
