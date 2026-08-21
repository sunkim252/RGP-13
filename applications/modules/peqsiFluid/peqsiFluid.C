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
#include "thermodynamicConstants.H"
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
    ladDtLimit_(great),
    srkReplicaValid_(false),
    srkChecked_(false),
    srkW_(0), srkB_(0), srkCoef1_(0), srkCoef2_(0), srkCoef3_(0),
    srkC_(0), srkCq0_(0), srkCq1_(0), srkCq2_(0), srkCTlo_(0), srkCThi_(0),

    initialMass_(-1),

    initialRhoH_(0)
{
    // ------------------------------------------------------------------
    // Single-species SRK replica coefficients for the constant-v closure
    // (thermoClosure.C).  Replicates the SRKGas dictionary constructor:
    // same formulas, same dictionary, cross-checked against the library
    // state at first use.
    // ------------------------------------------------------------------
    {
        const wordList& sp = thermo_.species();
        if (sp.size() == 1)
        {
            const dictionary& props =
                mesh.lookupObject<IOdictionary>("physicalProperties");

            if (props.found(sp[0]))
            {
                const dictionary& spDict = props.subDict(sp[0]);
                const dictionary& rf = spDict.subDict("rfProperties");

                srkW_ = spDict.subDict("specie").lookup<scalar>("molWeight");

                const scalar Tc = rf.lookup<scalar>("Tc");
                const scalar Pc = rf.lookup<scalar>("Pc");
                const scalar omega = rf.lookup<scalar>("omega");
                const scalar RR = constant::thermodynamic::RR;

                srkB_ = 0.08664*RR*Tc/Pc;
                const scalar a = 0.42747*sqr(RR*Tc)/Pc;
                const scalar S = 0.48508 + 1.5517*omega - 0.15613*sqr(omega);
                srkCoef1_ = a*sqr(1.0 + S);
                srkCoef2_ = a*2*S*(1 + S)/sqrt(Tc);
                srkCoef3_ = a*sqr(S)/Tc;

                if (rf.found("c"))
                {
                    srkC_ = rf.lookup<scalar>("c");
                }
                else if (rf.lookupOrDefault<Switch>("penelouxShift", false))
                {
                    const scalar Zra = 0.29056 - 0.08775*omega;
                    srkC_ = 0.40768*(0.29441 - Zra)*RR*Tc/Pc;
                }

                if (rf.found("penelouxCoeffs"))
                {
                    const scalarList pc
                    (
                        rf.lookup<scalarList>("penelouxCoeffs")
                    );
                    srkCq0_ = pc[0]; srkCq1_ = pc[1]; srkCq2_ = pc[2];
                    srkCTlo_ = pc[3]; srkCThi_ = pc[4];
                }

                srkReplicaValid_ = true;
            }
        }
    }

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

Foam::scalar Foam::solvers::peqsiFluid::maxDeltaT() const
{
    const scalar coDt = isothermalFluid::maxDeltaT();

    if (ladDtLimit_ < coDt)
    {
        Info<< "PEQSI dt limit: LAD explicit-diffusion bound "
            << ladDtLimit_ << " s (Courant bound " << coDt << " s)"
            << endl;
        return ladDtLimit_;
    }

    return coDt;
}


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
