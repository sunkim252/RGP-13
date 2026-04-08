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

#include "FGMTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FGMTable::FGMTable
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "fgmProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    pvAxis_(lookup("PV")),
    sourcePV_(lookup("sourcePV")),
    T_table_(lookup("T")),
    rho_table_(lookup("rho"))
{
    Info<< "\nFGM Table initialisation" << endl;
    Info<< "Table length: " << pvAxis_.size() << endl;
    Info<< "Table contents:" << endl;
    Info<< "{" << endl;
    Info<< "    PV (progress variable axis)" << endl;
    Info<< "    sourcePV (reaction source term)" << endl;
    Info<< "    T (temperature)" << endl;
    Info<< "    rho (density)" << endl;
    Info<< "}" << nl << endl;

    if (pvAxis_.size() < 2)
    {
        FatalErrorInFunction
            << "FGM table must have at least 2 entries, but has "
            << pvAxis_.size()
            << exit(FatalError);
    }

    if
    (
        sourcePV_.size() != pvAxis_.size()
     || T_table_.size() != pvAxis_.size()
     || rho_table_.size() != pvAxis_.size()
    )
    {
        FatalErrorInFunction
            << "FGM table arrays must all have the same size." << nl
            << "PV size: " << pvAxis_.size() << nl
            << "sourcePV size: " << sourcePV_.size() << nl
            << "T size: " << T_table_.size() << nl
            << "rho size: " << rho_table_.size()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FGMTable::~FGMTable()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::FGMTable::interpolate
(
    const List<scalar>& table,
    scalar pvValue
) const
{
    // Clamp PV to [0, 1]
    pvValue = min(pvValue, scalar(1));
    pvValue = max(pvValue, scalar(0));

    // Handle exact zero
    if (pvValue == 0)
    {
        return table[0];
    }

    // Small number to prevent divide by zero
    const scalar smallValue = 1e-5;

    scalar lower_pv = 0;
    scalar upper_pv = 0;
    scalar lower_val = 0;
    scalar upper_val = 0;

    // Find the bracketing interval in the PV axis
    for (label j = 0; j < pvAxis_.size(); j++)
    {
        if (pvAxis_[j] >= pvValue)
        {
            // Guard against j == 0 (should not happen since pvValue > 0
            // and pvAxis_[0] is typically 0, but be safe)
            label jLow = max(j - 1, label(0));

            lower_pv = pvAxis_[jLow];
            upper_pv = pvAxis_[j];

            lower_val = table[jLow];
            upper_val = table[j];

            break;
        }
    }

    const scalar rate =
        (upper_val - lower_val) / max(upper_pv - lower_pv, smallValue);

    scalar interpolatedValue = (pvValue - lower_pv) * rate + lower_val;

    // Clamp to table minimum
    interpolatedValue = max(interpolatedValue, min(table));

    return interpolatedValue;
}


// ************************************************************************* //
