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

Class
    Foam::FGMTable

Description
    Lookup table for Flamelet-Generated Manifold (FGM) combustion model.
    Reads 1D tables (PV, sourcePV, T, rho) from constant/fgmProperties
    and provides linear interpolation as a function of progress variable.

SourceFiles
    FGMTable.C

\*---------------------------------------------------------------------------*/

#ifndef FGMTable_H
#define FGMTable_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "scalar.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FGMTable Declaration
\*---------------------------------------------------------------------------*/

class FGMTable
:
    public IOdictionary
{
    // Private data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Progress variable axis
        List<scalar> pvAxis_;

        //- Source term for progress variable
        List<scalar> sourcePV_;

        //- Temperature table
        List<scalar> T_table_;

        //- Density table
        List<scalar> rho_table_;


public:

    // Constructors

        //- Construct from mesh
        FGMTable(const fvMesh& mesh);


    //- Destructor
    ~FGMTable();


    // Member Functions

        //- Interpolate a value from a table given a progress variable value
        scalar interpolate
        (
            const List<scalar>& table,
            scalar pv
        ) const;

        // Access

            //- Return const access to the PV axis
            inline const List<scalar>& pvAxis() const
            {
                return pvAxis_;
            }

            //- Return const access to the source PV table
            inline const List<scalar>& sourcePV() const
            {
                return sourcePV_;
            }

            //- Return const access to the temperature table
            inline const List<scalar>& T() const
            {
                return T_table_;
            }

            //- Return const access to the density table
            inline const List<scalar>& rho() const
            {
                return rho_table_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
