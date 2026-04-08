/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "FGM.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(FGM, 0);
    addToRunTimeSelectionTable(combustionModel, FGM, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::FGM::FGM
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    combustionModel(modelType, thermo, turb),
    table_(this->mesh()),
    pvName_
    (
        this->coeffs().lookupOrDefault<word>("progressVariable", "PV")
    ),
    pvIndex_(-1),
    sourcePV_
    (
        IOobject
        (
            "FGM:sourcePV",
            this->mesh().time().name(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimDensity/dimTime, 0)
    ),
    Qdot_
    (
        IOobject
        (
            "FGM:Qdot",
            this->mesh().time().name(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    )
{
    // Find the progress variable index in the species list
    const speciesTable& species = thermo.species();

    pvIndex_ = species[pvName_];

    Info<< "FGM combustion model:" << nl
        << "    Progress variable: " << pvName_ << nl
        << "    Species index: " << pvIndex_ << nl
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::FGM::~FGM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::FGM::correct()
{
    const volScalarField& PV = this->thermo().Y(pvIndex_);

    forAll(sourcePV_, celli)
    {
        const scalar pvVal = PV[celli];

        sourcePV_[celli] =
            table_.interpolate(table_.sourcePV(), pvVal);
    }

    // Step 1: Qdot = 0. Temperature evolution relies on species
    // transport and thermodynamic coupling. A proper Qdot based on
    // enthalpy differences can be added in a later step.
    Qdot_ = dimensionedScalar(dimEnergy/dimVolume/dimTime, 0);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::combustionModels::FGM::R(const label speciei) const
{
    if (speciei == pvIndex_)
    {
        return sourcePV_;
    }

    // For all other species, return zero source
    return
        volScalarField::Internal::New
        (
            typedName("R_" + this->thermo().Y()[speciei].name()),
            this->mesh(),
            dimensionedScalar(dimDensity/dimTime, 0)
        );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::FGM::R(volScalarField& Y) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));

    if (Y.name() == pvName_)
    {
        // Add explicit source: the fvScalarMatrix source convention is
        // that source() is on the RHS with a negative sign, i.e.
        //   [A][x] = [source]  =>  need source = rate * V
        tSu.ref().source() -= sourcePV_ * this->mesh().V();
    }

    return tSu;
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::FGM::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typedName("Qdot")),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    );
}


bool Foam::combustionModels::FGM::read()
{
    if (combustionModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
