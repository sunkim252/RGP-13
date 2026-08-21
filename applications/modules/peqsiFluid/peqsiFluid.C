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
#include "localEulerDdtScheme.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(peqsiFluid, 0);
    addToRunTimeSelectionTable(solver, peqsiFluid, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::peqsiFluid::peqsiFluid(fvMesh& mesh)
:
    isothermalFluid
    (
        mesh,
        autoPtr<fluidThermo>(fluidMulticomponentThermo::New(mesh).ptr())
    ),

    thermo_(refCast<fluidMulticomponentThermo>(isothermalFluid::thermo_)),

    h_
    (
        IOobject
        (
            "h",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        // Fallback for fresh starts: seed from the thermo enthalpy at the
        // initial (p, T) state; restarts read the transported field.
        volScalarField
        (
            IOobject("hInit", runTime.name(), mesh),
            isothermalFluid::thermo_.he()
        )
    ),

    fgmActive_
    (
        pimple.dict().lookupOrDefault<Switch>("peqsiFGM", false)
    ),

    alpha_
    (
        IOobject("PEQSI:alpha", runTime.name(), mesh),
        mesh,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    beta_
    (
        IOobject("PEQSI:beta", runTime.name(), mesh),
        mesh,
        // beta = (1/rho)(dp/dv)_T + ... : pressure dimensions;
        // beta/(1-alpha) == -rho c^2 [Pa]
        dimensionedScalar(dimPressure, -1),
        zeroGradientFvPatchScalarField::typeName
    ),

    acousticTimeIndex_(-1),

    initialMass_(-1),

    initialRhoH_(0)
{
    if (fv::localEulerDdt::enabled(mesh))
    {
        FatalErrorInFunction
            << "peqsiFluid is a transient fractional-step solver: "
            << "LTS (localEuler) is not supported"
            << exit(FatalError);
    }

    // Initial coefficient evaluation from the restart state
    updateCoefficients();

    if (fgmActive_)
    {
        // Stage 1: the state is declared and read, the manifold closure
        // and the multi-species pressure source are not wired yet --
        // refuse rather than run a half-coupled system silently.
        FatalErrorInFunction
            << "peqsiFGM is set but the FGM coupling is still at stage 1 "
            << "(state declaration only).  See the wiki section "
            << "'FGM coupling design' for the staged plan; run with "
            << "peqsiFGM off for the inert path."
            << exit(FatalError);
    }

    Info<< "peqsiFluid: PEQSI fractional-step solver "
        << "(Wada et al., Phys. Fluids 36 (2024) 116104)" << nl
        << "    advective substep: SSP-RK3, p advected (PEQSI Eq. 10)" << nl
        << "    acoustic substep: consistency-form Helmholtz (Eq. 19)" << nl
        << "    thermo closure: T from h(T,v) Newton inversion (WKK Fig. 3)"
        << nl
        << "    composition: single-species (FGM coupling: stage 1/6)"
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::peqsiFluid::~peqsiFluid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::peqsiFluid::preSolve()
{
    // Base preSolve handles the Courant number (transient branch) and
    // fvModels preUpdateMesh; the static mesh makes the rest a no-op.
    isothermalFluid::preSolve();
}


void Foam::solvers::peqsiFluid::prePredictor()
{
    // No-op: the base isothermalFluid::prePredictor() solves a continuity
    // predictor for rho ("rhoFinal" solve), which must NOT run -- density
    // is advanced exclusively by the fractional step (advective substep +
    // dp increment), which is what makes the scheme mass conservative.
}


void Foam::solvers::peqsiFluid::thermophysicalPredictor()
{
    // No-op: the thermodynamic closure runs after the acoustic substep
    // (correctPressurePEQSI), on the fully updated (rho, h, p) state.
}


void Foam::solvers::peqsiFluid::motionCorrector()
{
    // Static-mesh no-op: mesh motion is not supported.
}


void Foam::solvers::peqsiFluid::postCorrector()
{
    // No-op: the acoustic substep finalised the state; the base
    // implementation must not re-sync density or transport here.
    // TODO(V3): momentumTransport correct() for LES.
}


void Foam::solvers::peqsiFluid::postSolve()
{
    // Base postSolve does  rho_ = thermo.rho()  every step end -- the
    // solver's continuity density is SEPARATE storage from the thermo's
    // EOS density, and that re-sync silently replaced the transported
    // (mass-conserving) field with the EOS one (+0.8% mass/step measured
    // on case A before this override).  Keep only the cleanup.
    divrhoU.clear();
}


// ************************************************************************* //
