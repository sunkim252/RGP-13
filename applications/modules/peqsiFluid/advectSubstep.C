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

    const dimensionedScalar& dt = runTime.deltaT();

    // Snapshot the n-state (consumed by the acoustic substep)
    rhoN_.set(new volScalarField("PEQSI:rhoN", rho_));
    UN_.set(new volVectorField("PEQSI:UN", U_));
    pN_.set(new volScalarField("PEQSI:pN", p_));
    hN_.set(new volScalarField("PEQSI:hN", h_));

    // Frozen advecting volumetric flux and its divergence
    const surfaceScalarField phiv("PEQSI:phiv", fvc::flux(UN_()));
    const volScalarField divPhiv("PEQSI:divPhiv", fvc::div(phiv));

    // Frozen sources at the n-state:
    //
    // L_h: enthalpy-equation RHS without Dp/Dt (WKK Eq. 3) -- here the
    // conductive part div(kappa grad T); viscous heating is neglected as
    // in the reference low-Mach cryogenic-jet applications.
    const volScalarField kappa
    (
        completeField("PEQSI:kappa", mesh, thermo_.kappa())
    );
    const volScalarField mu
    (
        completeField("PEQSI:mu", mesh, thermo_.mu())
    );

    const volScalarField Lh
    (
        "PEQSI:Lh",
        fvc::laplacian(kappa, thermo_.T())
    );

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

    // RK working fields (rho h transported as a product)
    volScalarField r("PEQSI:rkRho", rhoN_());
    volVectorField ru("PEQSI:rkRhoU", rhoN_()*UN_());
    volScalarField pw("PEQSI:rkP", pN_());
    volScalarField rh("PEQSI:rkRhoH", rhoN_()*hN_());

    // One SSP-RK3 stage: q <- cOld*qn + cNew*(q + dt*L(q))
    auto stage = [&](const scalar cOld, const scalar cNew)
    {
        const volScalarField Lr(-uGrad(phiv, divPhiv, r, "div(phiv,rho)"));
        const volVectorField Lru
        (
            -uGrad(phiv, divPhiv, ru, "div(phiv,rhoU)") + divTau
        );
        const volScalarField Lp
        (
            -uGrad(phiv, divPhiv, pw, "div(phiv,p)") + aL*Lh
        );
        const volScalarField Lrh
        (
            -uGrad(phiv, divPhiv, rh, "div(phiv,rhoh)") + iL*Lh
        );

        r == cOld*rhoN_() + cNew*(r + dt*Lr);
        ru == cOld*(rhoN_()*UN_()) + cNew*(ru + dt*Lru);
        pw == cOld*pN_() + cNew*(pw + dt*Lp);
        rh == cOld*(rhoN_()*hN_()) + cNew*(rh + dt*Lrh);

        r.correctBoundaryConditions();
        ru.correctBoundaryConditions();
        pw.correctBoundaryConditions();
        rh.correctBoundaryConditions();
    };

    stage(0.0, 1.0);            // q1 = qn + dt L(qn)
    stage(0.75, 0.25);          // q2 = 3/4 qn + 1/4 (q1 + dt L(q1))
    stage(1.0/3.0, 2.0/3.0);    // q* = 1/3 qn + 2/3 (q2 + dt L(q2))

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

    // Intermediate mass flux (diagnostics; rebuilt after the acoustic
    // substep from the end-of-step state)
    phi_ = fvc::flux(rho_*U_);
}


// ************************************************************************* //
