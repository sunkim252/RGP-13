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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
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
        dimensionedScalar(sqr(dimVelocity), -1),
        zeroGradientFvPatchScalarField::typeName
    ),

    acousticTimeIndex_(-1)
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

    Info<< "peqsiFluid: PEQSI fractional-step solver "
        << "(Wada et al., Phys. Fluids 36 (2024) 116104)" << nl
        << "    advective substep: SSP-RK3, p advected (PEQSI Eq. 10)" << nl
        << "    acoustic substep: consistency-form Helmholtz (Eq. 19)" << nl
        << "    thermo closure: T from h(T,v) Newton inversion (WKK Fig. 3)"
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


void Foam::solvers::peqsiFluid::thermophysicalPredictor()
{
    // No-op: the thermodynamic closure runs after the acoustic substep
    // (correctPressurePEQSI), on the fully updated (rho, h, p) state.
}


void Foam::solvers::peqsiFluid::motionCorrector()
{
    // Static-mesh no-op: mesh motion is not supported.
}


// ************************************************************************* //
