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
#include "fvcFlux.H"
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

    // Mesh-alignment audit.  The directional devices -- LAD, the
    // selective and shock-capturing filters, the Z bounding
    // diffusivity's overbar -- decompose operators with per-face
    // direction weights (n . e_l)^2, which is a Cartesian-grid
    // construction: on faces not aligned with a coordinate axis the
    // "directional" second derivative smears across directions and the
    // device magnitudes are no longer the TK/Visbal-Gaitonde formulas.
    // Measure the alignment once and say so, rather than letting an
    // unstructured run discover it as a mis-scaled artificial term.
    {
        const surfaceVectorField nf(mesh.Sf()/mesh.magSf());
        scalar worst = 0;   // deviation of max |n.e_l| from 1
        const vectorField& nfi = nf.primitiveField();
        forAll(nfi, facei)
        {
            const vector& n = nfi[facei];
            const scalar amax =
                max(mag(n.x()), max(mag(n.y()), mag(n.z())));
            worst = max(worst, 1.0 - amax);
        }
        reduce(worst, maxOp<scalar>());
        if (worst > 0.01)
        {
            WarningInFunction
                << "mesh faces deviate from axis alignment by up to "
                << 100*worst << "% (1 - max|n.e_l|).  The directional "
                << "devices (peqsiLADCoeff, peqsiFilterSigma, "
                << "peqsiFilterSC, peqsiBoundZ) assume grid-aligned "
                << "hexes; on this mesh their magnitudes are not the "
                << "published formulas.  The core scheme (RK3 + "
                << "Helmholtz + updates) is mesh-general -- set "
                << "nNonOrthogonalCorrectors as for any p equation."
                << endl;
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
        // ------------------------------------------------------------
        // Stage 2a: manifold coordinates transported, composition and
        // the real-gas coefficient blocks looked up from the baked
        // table.  Scope: wiring, conservation and stability -- the
        // underlying table is the verification bake (dpm2), not
        // production physics.
        // ------------------------------------------------------------
        fgmTable_.reset(new FGMTable(mesh, "fgmProperties"));

        if (fgmTable_().hasChi())
        {
            Info<< "peqsiFGM: native chi-axis table (nChi = "
                << fgmTable_().nChi() << ") -- 4th coordinate is chi_st "
                << "from 2 (Deff/rho)|grad Z|^2 with the Pitsch-Steiner "
                << "Z(1-Z) mapping (fgmFluid convention).  ADIABATIC: "
                << "no dh correction is available from this table."
                << endl;
        }

        // The table must carry the PEQSI coefficient blocks and W;
        // refuse a plain combustion table -- the closure would then
        // silently run without the real-gas coefficients.
        for (const char* nm :
             {"PEQSI_xi", "PEQSI_alpha", "PEQSI_betan", "W"})
        {
            if (!fgmTable_().hasOptTable(nm))
            {
                FatalErrorInFunction
                    << "peqsiFGM: table lacks block '" << nm
                    << "' -- bake with bake_peqsi_coeffs.py"
                    << exit(FatalError);
            }
        }

        // Every table species must exist in the thermo so the looked-up
        // composition can be written into it
        forAll(fgmTable_().speciesNames(), i)
        {
            const word& sp = fgmTable_().speciesNames()[i];
            if (findIndex(thermo_.species(), sp) < 0)
            {
                FatalErrorInFunction
                    << "peqsiFGM: table species " << sp
                    << " missing from the thermo" << exit(FatalError);
            }
        }

        auto readField = [&](const word& nm) -> volScalarField*
        {
            return new volScalarField
            (
                IOobject
                (
                    nm, runTime.name(), mesh,
                    IOobject::MUST_READ, IOobject::AUTO_WRITE
                ),
                mesh
            );
        };
        Z_.reset(readField("Z"));
        Zvar_.reset(readField("Zvar"));
        Yc_.reset(readField("Yc"));

        sourceYc_.reset
        (
            new volScalarField
            (
                IOobject("PEQSI:sourceYc", runTime.name(), mesh),
                mesh,
                dimensionedScalar(dimMass/dimVolume/dimTime, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }

    // Restart consistency (peqsiRestartH, default "off"): a foreign
    // restart -- fgmFluid state, or a table-lineage change -- carries an
    // h that was NOT built against the composition/dh convention THIS
    // table uses.  Two failure modes, two reseed modes.
    //
    // "eos": run the closure once so Y is the manifold's, then reseed
    // h = he(p, T, Y_manifold) -- the case KEEPS its own T, only the
    // composition changes.  Fixes a Y-lineage mismatch (measured on the
    // rd0110 3M bring-up: ~486 injector-lip cells pinned at the
    // [50, 4000] K clamps, maxRel 0.49, flat over 28 steps).
    //
    // "manifold": the case's h encodes a NON-ADIABATIC dh -- e.g. the
    // per-inlet temperature spread (LOX/fuel streams enter below/above
    // the adiabatic mixing-line value) -- under a DIFFERENT table's
    // dhRef/hOx/hFuel convention.  "eos" alone does not fix this: it
    // keeps the case T, so the reseeded h still encodes the foreign
    // convention's inlet-temperature offset at THIS table's hOx/hFuel
    // anchors, which is not the same physical dh at all (measured:
    // 71.9% of cells sit above this table's dh axis by step 30, and the
    // upper-slice clamp is self-sustaining -- +517 K, flat).  "manifold"
    // discards the case dh entirely and lands exactly on THIS table's
    // own dh = 0 mixing-line state, h = hMix + dhRef(Z, gz, c), at the
    // case's transported (Z, gz, c).  The true per-inlet dh re-forms
    // from the boundary conditions themselves (each inlet patch is
    // fixedValue T) over the inlet residence time; only the INTERIOR
    // initial state changes.
    const word restartH
    (
        pimple.dict().lookupOrDefault<word>("peqsiRestartH", "off")
    );
    if (restartH != "off" && restartH != "eos" && restartH != "manifold")
    {
        FatalErrorInFunction
            << "peqsiRestartH must be 'off', 'eos' or 'manifold', got '"
            << restartH << "'" << exit(FatalError);
    }
    if (fgmActive_ && restartH != "off")
    {
        const word zvarMode
        (
            pimple.dict().lookupOrDefault<word>("peqsiZvar", "algebraic")
        );
        if (zvarMode == "algebraic")
        {
            updateSegregation();
        }
        fgmClosure();

        if (restartH == "eos")
        {
            h_ = thermo_.he();
            Info<< "PEQSI restart: h reseeded from he(p, T, Y_manifold) "
                << "after one closure pass" << endl;
        }
        else
        {
            FGMTable& tbl = fgmTable_();
            const scalar hOx = tbl.hOx();
            const scalar hFuel = tbl.hFuel();
            const bool haveDhRef = tbl.hasDhRef();

            const scalarField& Zf = Z_().primitiveField();
            const scalarField& gZf = Zvar_().primitiveField();
            const scalarField& cF = cNormF_();
            scalarField& hf = h_.primitiveFieldRef();

            forAll(hf, celli)
            {
                const scalar Zcl = min(max(Zf[celli], 0.0), 1.0);
                const scalar gz = max(gZf[celli], 0.0);
                const scalar Ccl = cF[celli];
                const scalar hMix = (1.0 - Zcl)*hOx + Zcl*hFuel;
                hf[celli] =
                    hMix
                  + (haveDhRef
                      ? tbl.interpolateDhRef(Zcl, gz, Ccl) : 0.0);
            }
            Info<< "PEQSI restart: h reseeded to the table's own "
                << "adiabatic mixing-line state (dh = 0) at the case's "
                << "(Z, gz, c) -- the case's non-adiabatic dh is "
                << "discarded, not carried across the table lineage"
                << endl;
        }
        h_.correctBoundaryConditions();

        // Full-state reseed (manifold mode only).  Measured on the
        // rd0110 3M restart: the foreign initial field was CONVERGED
        // under wrong inlet BCs (LOX 100 K instead of the canonical
        // 143.3 K), leaving a 340k-cell core at rho ~1120 kg/m3 --
        // 1.47x the table's entire oxidiser branch (761 at 143.3 K), a
        // state the manifold cannot represent at all.  Reseeding h
        // alone leaves that transported rho as a permanent 47% EOS
        // mismatch feeding the injector-lip dp/rho regeneration loop
        // (the thing that forced dt down to 1e-10).  Flushing it by
        // advection would take ~ms of physical time; reseed instead:
        // run the closure a SECOND time (the first h reseed moved the
        // dh coordinate, so composition and the temperature guess must
        // be re-looked-up at the new point), take the manifold
        // temperature outright, and rebuild rho and phi from the EOS at
        // that state.  What survives from the case: p, U, Z, c -- the
        // dynamical memory.  What is discarded: every thermodynamic
        // quantity of the off-anchor state.
        if (restartH == "manifold")
        {
            fgmClosure();

            volScalarField& Tw = const_cast<volScalarField&>(thermo_.T());
            scalarField& Tf = Tw.primitiveFieldRef();
            const scalarField& Tg = Tguess_();
            forAll(Tf, celli)
            {
                if (Tg[celli] > 0)
                {
                    Tf[celli] =
                        min
                        (
                            max
                            (
                                Tg[celli],
                                pimple.dict().lookupOrDefault<scalar>
                                (
                                    "peqsiTmin", 100
                                )
                            ),
                            scalar(4000)
                        );
                }
            }
            Tw.correctBoundaryConditions();

            rho_ = thermo_.rho();
            rho_.correctBoundaryConditions();
            phi_ = fvc::flux(rho_*U_);

            Info<< "PEQSI restart: T from the manifold and rho from the "
                << "EOS at the reseeded state; phi rebuilt.  Case p, U, "
                << "Z, c kept." << endl;
        }
    }

    Info<< "peqsiFluid: PEQSI fractional-step solver "
        << "(Wada et al., Phys. Fluids 36 (2024) 116104)" << nl
        << "    advective substep: SSP-RK3, p advected (PEQSI Eq. 10)" << nl
        << "    acoustic substep: consistency-form Helmholtz (Eq. 19)" << nl
        << "    thermo closure: T from h(p, T) inversion at the "
        << "transported p" << nl
        << "        (WKK Fig. 3 constant-v form: peqsiConstantV, "
        << "single-species only)" << nl
        << "    composition: "
        << (fgmActive_
            ? "manifold (FGM stage 2a: wiring verification)"
            : "single-species (FGM off)")
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
