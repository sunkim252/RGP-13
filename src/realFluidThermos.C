/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 RGP-13
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Explicit template instantiation for SRKchungTakaMixture with
    psiMulticomponentThermo and fluidMulticomponentThermo.
\*---------------------------------------------------------------------------*/

#include "psiMulticomponentThermo.H"
#include "rhoFluidMulticomponentThermo.H"
#include "SRKchungTakaMixture.H"
#include "forRealFluidGases.H"
#include "forRealFluidGasesElyHanley.H"
#include "makeFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // ---------------- chungTransport backend ----------------
    // psi-based (compressible)
    forRealFluidGases
    (
        makeFluidMulticomponentThermos,
        psiThermo,
        psiMulticomponentThermo,
        SRKchungTakaMixture
    );

    // rho-based (density-direct)
    forRealFluidGases
    (
        makeFluidMulticomponentThermos,
        rhoFluidThermo,
        rhoFluidMulticomponentThermo,
        SRKchungTakaMixture
    );

    // ---------------- elyHanleyTransport backend (Phase 2) ----------------
    // psi-based (compressible)
    forRealFluidGasesEH
    (
        makeFluidMulticomponentThermos,
        psiThermo,
        psiMulticomponentThermo,
        SRKchungTakaMixture
    );

    // rho-based (density-direct)
    forRealFluidGasesEH
    (
        makeFluidMulticomponentThermos,
        rhoFluidThermo,
        rhoFluidMulticomponentThermo,
        SRKchungTakaMixture
    );
}

// ************************************************************************* //
