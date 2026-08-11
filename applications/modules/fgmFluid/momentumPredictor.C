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

#include "fgmFluid.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "safeReconstruct.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::fgmFluid::momentumPredictor()
{
    volVectorField& U(U_);

    // --- Localized Artificial (shear) Viscosity on the momentum -------------
    // Companion to the scalar mass-diffusivity LAD: damps at source the
    // spurious VELOCITY overshoot at the transcritical density interface
    // (diagnosed recess-tip |U| spike, limitMag-capped at 500, that throttles
    // the time step). A density-gradient-sensed artificial dynamic viscosity
    //   muArt = LADUCoeff * V^(2/3) * |U| * |grad(rho)|   [kg/(m s)]
    // is added to the viscous stress (-fvm::laplacian(muArt, U)); cell-sized so
    // its own viscous time step stays of order the convective one and smooth
    // regions (|grad(rho)| ~ 0) are untouched -- the LES subgrid stress is
    // unaffected away from the interface. Cook & Cabot, J. Comput. Phys. 195
    // (2004); Kawai & Lele, J. Comput. Phys. 227 (2008); Kawai, Terashima &
    // Negishi, J. Comput. Phys. 300 (2015). LADUCoeff read from the PIMPLE
    // dict each step (runTimeModifiable); default 0 = off.
    const scalar LADUCoeff
    (
        pimple.dict().lookupOrDefault<scalar>("LADUCoeff", scalar(0))
    );
    volScalarField muArt
    (
        IOobject
        (
            "muArt",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    );
    if (LADUCoeff > 0)
    {
        const scalarField V23(pow(scalarField(mesh.V()), 2.0/3.0));
        muArt.primitiveFieldRef() =
            LADUCoeff*V23
           *mag(U)().primitiveField()
           *mag(fvc::grad(rho))().primitiveField();
        muArt.correctBoundaryConditions();
        Info<< "LAD-U: muArt max = " << gMax(muArt.primitiveField())
            << " kg/(m s)" << endl;
    }

    // Artificial BULK viscosity (Cook-Cabot) -- damps the DILATATIONAL /
    // compressive part directly, targeting the injector pressure oscillation
    // (p +/- spike and |U| overshoot at the fine tangential-hole cells) that
    // the shear muArt and the Rhie-Chow fixes do not reach. Dilatation-sensed
    //   betaArt = LADbulkCoeff * rho * V^(2/3) * |div U|   [kg/(m s)]
    // added as -grad(betaArt*div(U)) to the momentum. Cook & Cabot, J. Comput.
    // Phys. 195 (2004) 594; Kawai, Terashima & Negishi, J. Comput. Phys. 300
    // (2015) 116. Read each step (runTimeModifiable); default 0 = off.
    const scalar LADbulkCoeff
    (
        pimple.dict().lookupOrDefault<scalar>("LADbulkCoeff", scalar(0))
    );
    const volScalarField divU(fvc::div(U));
    volScalarField betaArt
    (
        IOobject
        (
            "betaArt",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    );
    if (LADbulkCoeff > 0)
    {
        const scalarField V23(pow(scalarField(mesh.V()), 2.0/3.0));
        betaArt.primitiveFieldRef() =
            LADbulkCoeff*rho.primitiveField()*V23
           *mag(divU.primitiveField());
        betaArt.correctBoundaryConditions();
        Info<< "LAD-bulk: betaArt max = " << gMax(betaArt.primitiveField())
            << " kg/(m s)" << endl;
    }

    // --- Korteweg capillary stress (transcritical density interface) --------
    // NOBLE (physics-restoring) alternative to the ad-hoc LAD/limitU regulari-
    // sation: restores the surface-tension-like force that a fully conservative
    // real-fluid code OMITS at the transcritical contact. Jofre & Urzay, Prog.
    // Energy Combust. Sci. 82 (2021) 100877: the SRK thermodynamic pressure
    // genuinely dips O(10-100 bar) across the diffuse LOX/gas interface (Maxwell
    // loop); in reality this is balanced by the Korteweg (van der Waals square-
    // gradient) stress
    //   tau_K = [lambda*rho*lap(rho) + 0.5*lambda*|grad(rho)|^2] I
    //         - lambda*grad(rho) (x) grad(rho)
    // whose divergence reduces by the standard identity to the body force
    //   f_K = div(tau_K) = lambda * rho * grad(lap(rho))          [N/m^3]
    // Added as a momentum source (RHS). It self-localises to the interface
    // (grad(lap(rho)) ~ 0 in smooth regions) so, unlike LAD, it does not diffuse
    // the resolved turbulence, and it OPPOSES the spurious force at the root
    // instead of damping the resulting velocity -- so limitU can be relaxed.
    // lambda [m^7/(kg s^2)] calibrated to the resolved interface thickness
    // (regularised: the true 5 nm interface is unresolvable). Default 0 = off.
    //
    // BALANCED-FORCE (2026-07-22): grad(lap(rho)) is reconstructed FROM the face
    // snGrad -- reconstruct(snGrad(lap(rho))*magSf) -- NOT the Green-Gauss cell
    // grad. This is the SAME discrete operator the faceGradP cell-velocity update
    // uses for -grad(p) (reconstruct(snGrad(p)*magSf), see pressureCorrector.C).
    // With the naive cell grad, the Korteweg body force and the face-consistent
    // pressure gradient live on DIFFERENT stencils, so at the transcritical
    // interface they do not cancel discretely and a residual parasitic velocity
    // survives (the classic surface-tension spurious-current defect; Francois
    // et al., JCP 213 (2006) 141; interFoam). Matching the operators makes
    // -snGrad(p) + lambda*rho*snGrad(lap(rho)) cancel face-by-face at mechanical
    // equilibrium, so the reconstructed cell velocity sees zero net force at the
    // interface -- the pressure dip is balanced, not merely opposed. Requires
    // faceGradP=true (asserted below) for the operators to match.
    const scalar kortewegLambda
    (
        pimple.dict().lookupOrDefault<scalar>("kortewegLambda", scalar(0))
    );
    volVectorField fKorteweg
    (
        IOobject
        (
            "fKorteweg",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimForce/dimVolume, Zero),
        zeroGradientFvPatchVectorField::typeName
    );
    if (kortewegLambda != 0)
    {
        const dimensionedScalar lambdaK
        (
            "lambdaK",
            dimensionSet(-1, 7, -2, 0, 0, 0, 0),
            kortewegLambda
        );
        // balanced-force: reconstruct grad(lap(rho)) from the face snGrad so it
        // shares the faceGradP -grad(p) stencil (see header note above).
        // safeReconstruct (not fvc::reconstruct): AMR hanging-node cells can
        // singularize the reconstruction tensor -- see safeReconstruct.H.
        const volScalarField lapRho(fvc::laplacian(rho));
        fKorteweg =
            lambdaK*rho
           *safeReconstruct(fvc::snGrad(lapRho)*mesh.magSf(), fvc::grad(lapRho));
        fKorteweg.correctBoundaryConditions();
        Info<< "Korteweg: |fK| max = "
            << gMax(mag(fKorteweg)().primitiveField()) << " N/m^3" << endl;
    }

    tUEqn =
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + momentumTransport->divDevTau(U)
      - fvm::laplacian(muArt, U)
      - fvc::grad(betaArt*divU)
     ==
        fvModels().source(rho, U)
      + fKorteweg
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    // --- ladRhoMomentum: Kawai Eq. (32) right-hand side ----------------------
    // The artificial mass diffusion A_rho = Dr*grad(rho) that pressureCorrector
    // adds to the (volumetric) continuity must also carry momentum:
    //     d(rho u)/dt + div(rho u (x) u + p d - tau) = div(A_rho (x) u)
    // With a uniform u both div(rho u (x) u) and div(A_rho (x) u) collapse to
    // u times the corresponding continuity terms, so the momentum equation
    // reduces to u*[continuity] and the uniform field is EXACTLY preserved
    // (Kawai, Terashima & Negishi, JCP 300 (2015) 116, Sec. 2.2). Omitting it
    // -- as this solver did -- means the mass smoothing injects momentum at
    // every density gradient. phiA_amd = Sf & A_rho is published by
    // pressureCorrector and is identically zero when LADrhoCoeff = 0.
    // NOTE the caveat: this solver solves no discrete continuity equation at
    // all (the pressure equation is the volumetric PEP form), so the
    // constant-velocity property is recovered only to the extent that the
    // pressure equation's -div(Dr grad rho)/rho mirrors continuity's +div(A_rho).
    const Switch ladRhoMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("ladRhoMomentum", false)
    );
    if (ladRhoMomentum && mesh.foundObject<surfaceScalarField>("phiA_amd"))
    {
        const surfaceScalarField& phiA =
            mesh.lookupObject<surfaceScalarField>("phiA_amd");
        UEqn -= fvm::div(phiA, U, "div(phiA,U)");
        Info<< "LAD-rho momentum: |phiA| max = "
            << gMax(mag(phiA.primitiveField())) << " kg/s" << endl;
    }

    UEqn.relax();

    fvConstraints().constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        if (buoyancy.valid())
        {
            solve
            (
                UEqn
             ==
                netForce()
            );
        }
        else
        {
            // faceGradP: face-consistent pressure gradient in the momentum
            // predictor, matching the pressureCorrector's balanced-force cell
            // velocity update (see pressureCorrector.C). Default off.
            const Switch faceGradP
            (
                pimple.dict().lookupOrDefault<Switch>("faceGradP", false)
            );
            if (faceGradP)
            {
                // safeReconstruct: see safeReconstruct.H (AMR hanging-node
                // singular-tensor SIGFPE guard, 2026-07-23). Falls back to
                // the same fvc::grad(p) the non-faceGradP branch below uses.
                solve
                (
                    UEqn
                 ==
                   -safeReconstruct(fvc::snGrad(p)*mesh.magSf(), fvc::grad(p))
                );
            }
            else
            {
                solve(UEqn == -fvc::grad(p));
            }
        }

        fvConstraints().constrain(U);
        K = 0.5*magSqr(U);
    }
}


// ************************************************************************* //
