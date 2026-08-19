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

    tmp<volScalarField> tRhs;
    if (consistencyRHS)
    {
        tRhs =
            (2.0/dt)
           *(
                rhoN_()*fvc::div(UN_())
              + rhoStar*fvc::div(UStar)
            );
    }
    else
    {
        tRhs = (4.0/dt)*sComp_();
    }

    fvScalarMatrix dpEqn
    (
        fvm::laplacian(dp)
      + fvm::div(Fdp, dp)
      + fvm::Sp(4.0*coef/sqr(dt), dp)
     ==
      - fvc::laplacian(pStar + pN_())
      + tRhs()
    );

    dpEqn.solve();

    Info<< "PEQSI: dp min/max = "
        << gMin(dp.primitiveField()) << " / "
        << gMax(dp.primitiveField()) << " Pa" << endl;

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

    // End-of-step mass flux
    phi_ = fvc::flux(rho_*U_);

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
    updateCoefficients();

    // Release the substep state
    rhoN_.clear();
    UN_.clear();
    pN_.clear();
    hN_.clear();
    sComp_.clear();
}


// ************************************************************************* //
