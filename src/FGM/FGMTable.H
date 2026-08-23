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
    Lookup table for the diffusion-flame (FPV + presumed beta-PDF) FGM
    combustion model.

    Reads a 3-D or 4-D table of the Favre-filtered progress-variable source
    term and the tabulated state from constant/fgmProperties.

    Axes:
        Z~   = Favre mean mixture fraction,
        gZ   = Zvar/(Z~(1-Z~)) = mixture-fraction segregation factor [0,1),
        C~   = Favre mean progress variable,
        chi  = scalar dissipation rate evaluated at Z = Z_st [1/s]
               (4-D version only; legacy 3-D tables omit this axis).

    The dictionary provides the axis sizes (nZ, nGz, nC, [nChi]), the axes
    (Z, gZ, C, [chi]) and the flat tables. C-order indexing:

        4-D:  idx = ((iZ*nGz + iGz)*nC + iC)*nChi + iChi
        3-D:  idx =  (iZ*nGz + iGz)*nC + iC        (legacy)

    Interpolation is quadrilinear in 4-D, trilinear in 3-D. The legacy
    interpolate(Z, gZ, C) entry points still work and evaluate the table
    at chi = chi_axis_[0] when the table is 4-D.

SourceFiles
    FGMTable.C

\*---------------------------------------------------------------------------*/

#ifndef FGMTable_H
#define FGMTable_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "scalar.H"
#include "List.H"
#include "HashTable.H"
#include "wordList.H"

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

        //- Axis sizes
        label nZ_;
        label nGz_;
        label nC_;
        label nChi_;            //!< 1 if the table is legacy 3-D

        //- True if a chi axis was provided in fgmProperties
        bool hasChi_;

        //- True if the 4th axis is an enthalpy-defect axis (non-adiabatic FPV,
        //  method b) rather than scalar dissipation. When set, the 4th
        //  coordinate passed to interpolate() is the enthalpy defect
        //  dh = h_total - ((1-Z) hOx + Z hFuel) and the same quadrilinear
        //  machinery is reused (chi_axis_ then holds the dh axis [J/kg]).
        bool useEnthalpy_;

        //- True if the 4th axis is a steam-dilution axis W = steam mole
        //  fraction in the oxidiser stream (power-generation H2/O2/H2O line).
        //  When set, the 4th coordinate passed to interpolate() is a
        //  transported conserved dilution scalar W (not a defect): the same
        //  quadrilinear machinery is reused (chi_axis_ then holds the W axis).
        bool useDilution_;

        //- Pure-stream total enthalpies at the flamelet boundary T [J/kg],
        //  used to form the adiabatic reference mixing line h_ad(Z).
        scalar hOx_;
        scalar hFuel_;

        //- Mixture-fraction axis (Z~)
        List<scalar> Z_axis_;

        //- Segregation-factor axis (gZ = Zvar/(Z~(1-Z~)))
        List<scalar> gZ_axis_;

        //- Progress-variable axis (C~)
        List<scalar> C_axis_;

        //- Scalar-dissipation axis (chi_st [1/s]); single entry {0} for
        //  legacy 3-D tables.
        List<scalar> chi_axis_;

        //- Flat PV source table.  Length nZ*nGz*nC*nChi.
        List<scalar> sourcePV_;

        //- Optional tabulated temperature [K] (empty if not present).
        List<scalar> T_table_;

        //- Optional pre-tabulated isentropic compressibility psis
        //  (empty if not present, 2026-07-24).
        List<scalar> psis_table_;

        //- Optional adiabatic-enthalpy offset of the manifold from the
        //  conserved-scalar mixing line (empty if not present, 2026-08-06):
        //      dhRef(Z,gZ,c) = h(T_ad, Y)|manifold - ((1-Z) hOx + Z hFuel)
        //  With differential diffusion the adiabatic flamelet does NOT lie on
        //  the mixing line (measured: -0.67 MJ/kg at Z_st for SRK+HP-Chung,
        //  +1.40 MJ/kg for mixture-averaged -- 43-91% of the whole dh axis),
        //  while the table's dh=0 slice IS that adiabatic state. Without this
        //  offset the solver's dh = h - h_mix reads a non-zero defect for a
        //  perfectly adiabatic cell and slides down the axis. Replicated along
        //  the 4th axis so the standard quadrilinear machinery applies.
        List<scalar> dhRef_table_;

        //- Equilibrium value of the RAW progress variable at each Z axis node,
        //  Yc_eq(Z) = sum_k w_k Y_k at chemical equilibrium (length nZ). This
        //  is the normalisation the table's c = Yc/Yc_eq(Z) uses. Optional;
        //  needed only when the solver transports the UNNORMALISED Yc.
        List<scalar> Cnorm_table_;

        //- Optional per-node tables read if present (PEQSI real-gas
        //  coefficient blocks and the mean molecular weight).  Absent
        //  in the fgmFluid production tables, so reading them is inert
        //  there; the peqsiFluid FGM stage consumes them.
        HashTable<List<scalar>> optTables_;

        //- Table reference pressure [Pa] (scalar entry "pressure",
        //  optional): the state the pressure-sensitive coefficients
        //  were baked at.  The dimensionless forms rescale with the
        //  cell pressure; see the PEQSI closure.
        scalar pRef_;

        //- Tabulated species mass fractions, keyed by species name.
        //  Each entry is a flat list of length nZ*nGz*nC*nChi.
        HashTable<List<scalar>> Y_tables_;

        //- Names of the tabulated species (order as in the 'species' entry)
        wordList speciesNames_;

        //- Tier-2: optional pre-tabulated real-gas mixture coefficients, in
        //  tabulatedRealGasMixture::coeffNames() order (table key
        //  "RG_"+name). Each entry is a flat list of length nZ*nGz*nC*nChi.
        //  Empty when the table does not carry the RG_* fields.
        List<List<scalar>> RGcoeff_;

        //- Tier-4: optional differential-diffusion Lewis-number fields, keyed
        //  by control-variable name ("Z", "C"); dictionary key "Le_"+name.
        //  Each entry is a flat list of length nZ*nGz*nC*nChi. When present the
        //  solver applies a per-cell Le(Z,gZ,c[,chi]) to that variable's
        //  laminar diffusivity (rho*D = mu/Le), generalising the constant 'Le'
        //  sub-dict to the standard differential-diffusion FGM closure.
        HashTable<List<scalar>> Le_tables_;


    // Binary sidecar cache (constructor cost).  Parsing the ASCII dict is
    // the whole startup: 7.7 s at 388 MB, ~30 s at 1.5 GB, paid by EVERY
    // rank on EVERY run.  After one full parse the constructor writes two
    // sidecars next to the RESOLVED dict file (so cases that symlink the
    // same table share one cache):
    //   <dict>.stub  every entry EXCEPT the node-count-length arrays --
    //                axes, scalars, species, anything a consumer may
    //                lookup() on the dictionary -- so the dictionary
    //                interface stays complete;
    //   <dict>.bin   the big arrays, raw doubles with a header recording
    //                the source file's size and mtime.
    // A later construction whose probe matches loads stub + bin instead
    // of parsing.  A stale cache cannot be read silently: the header is
    // checked against the CURRENT source file (realpath'd, so replacing
    // the symlink target invalidates), and a mismatch falls back to the
    // full parse and rewrites both sidecars.  Delete the .bin to opt out.

        bool cacheActive_ = false;
        fileName cacheSrc_;
        HashTable<List<scalar>> cacheBlocks_;

        //- The dict file the case actually points at (symlinks resolved)
        static fileName resolvedDictPath
        (
            const fvMesh& mesh,
            const word& dictName
        );

        //- True if <src>.stub and <src>.bin exist and the bin header
        //  matches the current source file
        static bool probeCache(const fileName& src);

        //- Read every block of <src>.bin into cacheBlocks_
        void loadCacheBin();

        //- Write the sidecars (master rank only; best-effort)
        void writeCache() const;

        //- One big block: from the cache when active, else lookup()
        List<scalar> readBig(const word& key);

        //- Presence of a big block: cache when active, else found()
        bool foundBig(const word& key) const;


    // Packed hot-table block.  The closure reads ~50 tables per cell at
    // one shared stencil; as separate arrays that is ~424 scattered
    // cache lines per DISTINCT stencil (16 corners x 53 arrays, h-pairs
    // sharing lines).  Interleaving the hot tables node-major makes each
    // corner one contiguous 53-value read: 16 x 7 = 112 lines, ~3.8x
    // fewer.  Invisible on a small case whose visited stencils fit L3;
    // the target is the production-scale run whose visited-node set does
    // not (the table itself is 231 MB at 506k nodes).  Costs a full
    // duplicate of the packed tables in memory.

        //- Interleaved copy: packedTables_[node*nPackedT_ + t]
        List<scalar> packedTables_;

        //- Keys in slot order (t of the layout above)
        wordList packedKeys_;

        label nPackedT_ = 0;

        //- Resolve a big-table key to member storage (nullptr if absent)
        const List<scalar>* tablePtr(const word& key) const;

public:

        //- Build the interleaved copy of the named tables.  A name that
        //  resolves to no table is a FatalError -- the caller's slot
        //  indices are positions in THIS list, so a silent skip would
        //  shift every later slot.
        void armPacked(const UList<word>& keys);

        inline bool packedArmed() const { return nPackedT_ > 0; }

        inline label nPackedTables() const { return nPackedT_; }

private:


    // Private Member Functions

        //- Quadrilinear interpolation of an arbitrary flat 4-D table.
        scalar interpolateTable
        (
            const List<scalar>& table,
            scalar Z, scalar gZ, scalar C, scalar chi
        ) const;

        //- Locate value v in a sorted axis; returns lower index i and sets
        //  the linear weight w in [0,1] such that
        //  result = (1-w)*axis[i] + w*axis[i+1]. Clamps to the ends. For
        //  axes of length 1 the call is degenerate: returns i=0, w=0 so
        //  the multilinear interpolation collapses to that slice.
        //- 비균등 축 가속용 균등 버킷 인덱스.
        //  buck[k] = axis 에서 (a0 + span*k/M) 이하인 최대 인덱스.
        //  ★왜: bracket 의 O(1) 빠른 경로는 축이 균등이라는 가정에 기댄다.
        //  Z 축을 |dT/dZ| 등분포로 잡으면 밀도비가 70배까지 벌어져 초기추정이
        //  크게 빗나가고 선형 걷기가 최대 44회로 늘어난다(실측). 버킷을 두면
        //  0.10회(최대 1회)로 돌아온다 -- 균등 축과 사실상 동일. 이분탐색(6.5회)
        //  보다 낫고 캐시 지역성도 유지된다. 메모리는 축당 ~4 KB.
        //  균등 축에서도 같은 셀을 찾으므로 결과는 비트 동일하다.
        List<label> Z_buck_;
        List<label> gZ_buck_;
        List<label> C_buck_;
        List<label> chi_buck_;

        //- 버킷 생성. M = ceil(span/min(dx)) 로 축 셀당 버킷 1개 이상 보장.
        static void buildBucket(const List<scalar>& axis, List<label>& buck);

        void bracket
        (
            const List<scalar>& axis,
            const List<label>& buck,
            scalar v,
            label& i,
            scalar& w
        ) const;

        //- Flat index into the 4-D table.
        inline label flatIndex
        (
            label iZ, label iGz, label iC, label iChi
        ) const
        {
            return ((iZ*nGz_ + iGz)*nC_ + iC)*nChi_ + iChi;
        }


public:

    //- Shared interpolation stencil: the 4-axis bracket + 16 corner flat
    //  indices of ONE query point (Z,gZ,C,chi), reusable across every
    //  tabulated field. updateManifold evaluates ~120+ fields (106 Y + 13
    //  RG + Le + T + source) at the SAME point per cell; sharing the
    //  bracket search + index computation removes that per-field
    //  duplication. Corner order matches interpolateTable's c0000..c1111,
    //  and interpolate(table, stencil) reproduces its nested-blend
    //  arithmetic verbatim -- results are bit-identical.
    struct FGMStencil
    {
        label idx[16];
        scalar wZ, wG, wC, wK;
    };


    // Constructors

        //- Construct from mesh (reads constant/fgmProperties)
        FGMTable(const fvMesh& mesh);

        //- Construct from mesh, reading a NAMED dictionary. Lets the solver
        //  load a SECOND manifold alongside the main one -- e.g. a
        //  (Z,gZ,C,chi_st) source-only UFPV table (constant/fgmPropertiesChi)
        //  next to the (Z,gZ,C,dh) non-adiabatic thermo table (2026-07-28).
        //  Same reader/validation path as the main constructor; only the
        //  dictionary name differs.
        FGMTable(const fvMesh& mesh, const word& dictName);


    //- Destructor
    ~FGMTable();


    // Member Functions

        // 4-D entry points (preferred for new clients)

            //- Quadrilinear interpolation of the PV source term.
            scalar interpolate
            (
                scalar Z, scalar gZ, scalar C, scalar chi
            ) const;

            //- Quadrilinear interpolation of tabulated temperature [K].
            scalar interpolateT
            (
                scalar Z, scalar gZ, scalar C, scalar chi
            ) const;

            //- Quadrilinear interpolation of a tabulated species mass
            //  fraction.
            scalar interpolateY
            (
                const word& specie,
                scalar Z, scalar gZ, scalar C, scalar chi
            ) const;

            //- Quadrilinear interpolation of the Tier-2 real-gas mixture
            //  coefficients. Fills coeffs (resized to the number of RG_*
            //  fields) in coeffNames() order. No-op if hasRealGasCoeffs() is
            //  false.
            void interpolateRealGasCoeffs
            (
                scalar Z, scalar gZ, scalar C, scalar chi,
                List<scalar>& coeffs
            ) const;

            //- Quadrilinear interpolation of the Tier-4 tabulated Lewis number
            //  Le_<var> for a control variable. Caller must check
            //  hasLeField(var) first.
            scalar interpolateLe
            (
                const word& var,
                scalar Z, scalar gZ, scalar C, scalar chi
            ) const;

            //- Quadrilinear interpolation of the pre-tabulated isentropic
            //  compressibility psis=psi/(rho*gamma) (offline SRK+JANAF
            //  re-derivation, gamma-clipped at build time -- see
            //  build_psis_table_v2.py, 2026-07-24). Caller must check
            //  hasPsis() first.
            scalar interpolatePsis
            (
                scalar Z, scalar gZ, scalar C, scalar chi
            ) const;

        // Shared-stencil entry points (fast path for many-field queries)

            //- Build the interpolation stencil for one query point.
            void makeStencil
            (
                scalar Z, scalar gZ, scalar C, scalar chi,
                FGMStencil& st
            ) const;

            //- Evaluate an arbitrary flat table on a prebuilt stencil
            //  (bit-identical to interpolateTable at the same point).
            scalar interpolate
            (
                const List<scalar>& table,
                const FGMStencil& st
            ) const;

            //- Evaluate ALL packed tables at one stencil into
            //  vals[0..nPackedTables()-1], slot order = armPacked key
            //  order.  Same nested blend, same floating-point order as
            //  interpolate(): bit-identical values.
            void interpolatePacked(const FGMStencil& st, scalar* vals) const;

            //- Direct table access for hoisting out of cell loops.
            inline const List<scalar>& sourcePVTable() const
            { return sourcePV_; }

            inline const List<scalar>& Ttable() const
            { return T_table_; }

            inline const List<scalar>& psisTable() const
            { return psis_table_; }

            inline const List<scalar>& Ytable(const word& sp) const
            { return Y_tables_[sp]; }

            inline label nRGcoeffs() const { return RGcoeff_.size(); }

            inline const List<scalar>& RGtable(const label k) const
            { return RGcoeff_[k]; }

            inline const List<scalar>& LeTable(const word& var) const
            { return Le_tables_[var]; }

            //- The 4th-axis coordinate a legacy 3-D query collapses to.
            inline scalar chi0() const { return chi_axis_[0]; }

            //- The c-axis node interval bracketing C: [cLo, cHi].
            //  The manifold is piecewise-linear in c between two nodes, so a
            //  difference taken across this interval IS the interpolant's
            //  local slope -- exactly, not approximately.  A free-parameter
            //  step is not equivalent: where the table is steep in c the
            //  slope jumps at the nodes, and a step small enough to sit
            //  inside one interval measures a different number than one that
            //  spans several.  Measured 2026-08-23 at a matched state,
            //  shrinking a fixed step from 0.02 to 0.00125 moved the
            //  derivative by ~8% -- small, but it is 8% of arbitrariness
            //  that no one can defend a value for.  (An earlier reading of
            //  the same sweep claimed the derivative diverged as 1/dc; that
            //  compared max|S_Y| over whole runs which had already drifted
            //  to different states, so it compared different flows, not
            //  different steps.)  The node interval is the finest structure
            //  the table actually carries, so asking for a slope below it
            //  asks for information that is not there.
            void cInterval(scalar C, scalar& cLo, scalar& cHi) const;

            //- Opt-1 (base-blend node interpolation): return the 16 surrounding
            //  manifold node flat indices and their quadrilinear weights for
            //  (Z,gZ,C,chi). nodes[m] and weights[m], m=0..15; sum(weights)=1.
            //  For a 3-D table the chi axis collapses (the chi-high corners
            //  duplicate the chi-low ones with weight 0). Lets a client blend
            //  pre-tabulated node quantities instead of summing over species.
            void interpStencil
            (
                scalar Z, scalar gZ, scalar C, scalar chi,
                label nodes[16], scalar weights[16]
            ) const;

            //- Opt-1: total number of manifold nodes (nZ*nGz*nC*nChi).
            inline label nTot() const
            {
                return nZ_*nGz_*nC_*nChi_;
            }

            //- Opt-1: the flat tabulated mass-fraction node list of a species
            //  (length nTot), used to pre-blend the per-node base thermo.
            inline const List<scalar>& Ynodes(const word& specie) const
            {
                return Y_tables_[specie];
            }

            //- Opt-1: true if the table tabulates this species' mass fraction.
            inline bool hasY(const word& specie) const
            {
                return Y_tables_.found(specie);
            }

        // 3-D legacy entry points
        //  Evaluate at chi = chi_axis_[0]. Required so that the original
        //  FGM combustion model and any existing case that still ships a
        //  3-D fgmProperties keeps working unchanged.

            scalar interpolate (scalar Z, scalar gZ, scalar C) const;
            scalar interpolateT(scalar Z, scalar gZ, scalar C) const;
            scalar interpolateY
            (
                const word& specie, scalar Z, scalar gZ, scalar C
            ) const;
            scalar interpolatePsis(scalar Z, scalar gZ, scalar C) const;

        // Predicates

            inline bool hasSpecie(const word& specie) const
            { return Y_tables_.found(specie); }

            inline bool hasT() const { return !T_table_.empty(); }

            //- True if the table carries the pre-tabulated psisTab field.
            inline bool hasPsis() const { return !psis_table_.empty(); }

            //- True if the table carries the Tier-2 real-gas coefficient
            //  fields (RG_*). The solver then enables the per-cell coefficient
            //  lookup instead of the live O(n^2) mixing.
            inline bool hasRealGasCoeffs() const
            { return !RGcoeff_.empty(); }

            //- True if the table carries a tabulated Lewis number Le_<var>
            //  (Tier-4 differential diffusion). The solver then uses a per-cell
            //  Le(Z,gZ,c[,chi]) for that variable instead of a constant.
            inline bool hasLeField(const word& var) const
            { return Le_tables_.found(var); }

            //- True if the table carries a chi_st axis (nChi > 1).
            inline bool hasChi() const
            { return hasChi_ && !useEnthalpy_ && !useDilution_; }

            //- True if the 4th axis is an enthalpy-defect axis (non-adiabatic
            //  FPV). The solver then transports total enthalpy and forms the
            //  defect coordinate; the table interpolation is unchanged.
            inline bool useEnthalpy() const { return useEnthalpy_; }

            //- True if the 4th axis is a steam-dilution axis W. The solver
            //  then transports a conserved dilution scalar W and passes it
            //  directly as the 4th coordinate; interpolation is unchanged.
            inline bool useDilution() const { return useDilution_; }

            //- Oxidizer-stream total enthalpy at the flamelet boundary [J/kg].
            inline scalar hOx() const { return hOx_; }

            //- Fuel-stream total enthalpy at the flamelet boundary [J/kg].
            inline scalar hFuel() const { return hFuel_; }

            //- True if the table ships the adiabatic-enthalpy offset dhRef.
            inline bool hasDhRef() const { return dhRef_table_.size() > 0; }

            //- True if the table ships Yc_eq(Z) (the c normalisation).
            inline bool hasCnorm() const { return Cnorm_table_.size() > 0; }

            inline bool hasOptTable(const word& nm) const
            { return optTables_.found(nm); }

            inline const List<scalar>& optTable(const word& nm) const
            { return optTables_[nm]; }

            inline scalar pRef() const { return pRef_; }

            //- Yc_eq(Z): equilibrium value of the raw progress variable.
            //  Linear interpolation on the Z axis, clamped at both ends.
            //  Lets the solver transport Yc (whose transport equation has NO
            //  normalisation cross-terms) and normalise only for the lookup.
            scalar interpolateCnorm(const scalar Z) const;

            //- dhRef(Z,gZ,c) [J/kg]: add to the mixing line to obtain the
            //  manifold's own adiabatic enthalpy, so that
            //      dh = h - ((1-Z) hOx + Z hFuel) - dhRef
            //  is zero exactly on the adiabatic manifold. Caller must check
            //  hasDhRef(). The stored field is 4th-axis-replicated, so any
            //  4th coordinate returns the same value; chi_axis_[0] is used.
            scalar interpolateDhRef(scalar Z, scalar gZ, scalar C) const;

        // Access

            inline const List<scalar>& Zaxis()   const { return Z_axis_;   }
            inline const List<scalar>& gZaxis()  const { return gZ_axis_;  }
            inline const List<scalar>& Caxis()   const { return C_axis_;   }
            inline const List<scalar>& chiAxis() const { return chi_axis_; }
            inline const List<scalar>& sourcePV()const { return sourcePV_; }
            inline const wordList& speciesNames() const
            { return speciesNames_; }

            inline label nZ()   const { return nZ_;   }
            inline label nGz()  const { return nGz_;  }
            inline label nC()   const { return nC_;   }
            inline label nChi() const { return nChi_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
