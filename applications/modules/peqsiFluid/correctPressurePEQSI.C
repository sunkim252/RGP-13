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

    // dp with zero-gradient boundaries (constraint patches override);
    // outflow/non-reflecting treatment is introduced with the 2-D jet
    // stage -- the 1-D validation cases are periodic.
    volScalarField dp
    (
        IOobject("PEQSI:dp", runTime.name(), mesh),
        mesh,
        dimensionedScalar(dimPressure, 0),
        zeroGradientFvPatchScalarField::typeName
    );

    // Implicit convective-coefficient face flux:
    //   d/dx_i ( 2 rho^n u^n_i (1-alpha)/beta dp / dt )
    const surfaceScalarField Fdp
    (
        "PEQSI:Fdp",
        fvc::flux(2.0*coef*UN_())/dt
    );

    // Modified Helmholtz, PEQSI Eq. (19): the RHS carries the
    // consistency-substituted volumetric source
    //   (4/dt) * (1/2) ( rho^n div u^n + rho^* div u^* )
    // in place of WKK Eq. (28)'s discrete continuity residual.
    fvScalarMatrix dpEqn
    (
        fvm::laplacian(dp)
      + fvm::div(Fdp, dp)
      + fvm::Sp(4.0*coef/sqr(dt), dp)
     ==
      - fvc::laplacian(pStar + pN_())
      + (2.0/dt)
       *(
            rhoN_()*fvc::div(UN_())
          + rhoStar*fvc::div(UStar)
        )
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
    // Mass audit: the fractional step conserves mass by construction;
    // report the residual of the discrete budget it satisfies.
    // ------------------------------------------------------------------
    {
        const volScalarField contErr
        (
            (rho_ - rhoN_())/dt + fvc::div(phi_)
        );
        Info<< "PEQSI mass balance: int(contErr) dV = "
            << gSum
               (
                   contErr.primitiveField()*mesh.V().primitiveField()
               )
            << " kg/s" << endl;
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
}


// ************************************************************************* //
