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
#include "SRKchungTakaMixture.H"
#include "forRealFluidGases.H"
#include "makeFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forRealFluidGases
    (
        makeFluidMulticomponentThermos,
        psiThermo,
        psiMulticomponentThermo,
        SRKchungTakaMixture
    );
}

// ************************************************************************* //
