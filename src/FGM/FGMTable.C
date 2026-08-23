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
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include "tabulatedRealGasMixture.H"
#include "IFstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"

// * * * * * * * * * * * * Binary sidecar cache  * * * * * * * * * * * * * * //

namespace
{
    // Bump on any format change: an old-format .bin must read as invalid.
    const char fgmBinMagic[8] = {'F','G','M','B','I','N','0','1'};

    struct FGMBinHeader
    {
        char magic[8];
        int64_t srcSize;
        int64_t srcMTime;
        int32_t nBlocks;
    };

    bool readBinHeader(const Foam::fileName& bin, FGMBinHeader& h)
    {
        std::ifstream is(bin.c_str(), std::ios::binary);
        if (!is) return false;
        is.read(reinterpret_cast<char*>(&h), sizeof(h));
        return bool(is);
    }
}


Foam::fileName Foam::FGMTable::resolvedDictPath
(
    const fvMesh& mesh,
    const word& dictName
)
{
    // In a decomposed run Time::path() is <case>/processorN while the
    // table lives in the undecomposed <case>/constant -- try both.
    for
    (
        const fileName& p :
        {
            mesh.time().path()/mesh.time().constant()/dictName,
            mesh.time().path()/".."/mesh.time().constant()/dictName
        }
    )
    {
        char buf[4096];
        if (::realpath(p.c_str(), buf))
        {
            return fileName(buf);
        }
    }
    return mesh.time().path()/mesh.time().constant()/dictName;
}


bool Foam::FGMTable::probeCache(const fileName& src)
{
    if (!isFile(src) || !isFile(src + ".stub") || !isFile(src + ".bin"))
    {
        return false;
    }

    FGMBinHeader h;
    if (!readBinHeader(src + ".bin", h)) return false;

    return
        memcmp(h.magic, fgmBinMagic, 8) == 0
     && h.srcSize == int64_t(Foam::fileSize(src))
     && h.srcMTime == int64_t(Foam::lastModified(src))
     && h.nBlocks >= 0;
}


void Foam::FGMTable::loadCacheBin()
{
    const fileName bin(cacheSrc_ + ".bin");
    std::ifstream is(bin.c_str(), std::ios::binary);

    FGMBinHeader h;
    is.read(reinterpret_cast<char*>(&h), sizeof(h));

    for (int32_t b = 0; b < h.nBlocks && is; b++)
    {
        int32_t klen = 0;
        is.read(reinterpret_cast<char*>(&klen), sizeof(klen));
        if (!is || klen <= 0 || klen > 256) break;
        std::string k(klen, '\0');
        is.read(&k[0], klen);
        int64_t n = 0;
        is.read(reinterpret_cast<char*>(&n), sizeof(n));
        if (!is || n < 0) break;
        List<scalar> block(static_cast<label>(n));
        is.read
        (
            reinterpret_cast<char*>(block.begin()),
            std::streamsize(n)*std::streamsize(sizeof(scalar))
        );
        if (!is) break;
        cacheBlocks_.insert(word(k), move(block));
    }

    if (!is || cacheBlocks_.size() != label(h.nBlocks))
    {
        // A truncated bin must not half-load: fall back to the parse.
        FatalErrorInFunction
            << "binary table cache " << bin << " is truncated ("
            << cacheBlocks_.size() << " of " << h.nBlocks
            << " blocks) -- delete it and rerun" << exit(FatalError);
    }

    Info<< "    binary table cache HIT: " << bin.name() << " ("
        << h.nBlocks << " blocks)" << endl;
}


Foam::List<Foam::scalar> Foam::FGMTable::readBig(const word& key)
{
    if (cacheActive_)
    {
        auto iter = cacheBlocks_.find(key);
        if (iter == cacheBlocks_.end())
        {
            FatalErrorInFunction
                << "block '" << key << "' missing from the binary table "
                << "cache " << cacheSrc_ + ".bin"
                << " -- delete the .bin/.stub sidecars and rerun"
                << exit(FatalError);
        }
        List<scalar> block(move(iter()));
        cacheBlocks_.erase(iter);
        return block;
    }
    return List<scalar>(lookup(key));
}


bool Foam::FGMTable::foundBig(const word& key) const
{
    return cacheActive_ ? cacheBlocks_.found(key) : found(key);
}


void Foam::FGMTable::writeCache() const
{
    if (!Pstream::master()) return;

    // Every node-count-length array, enumerated from the MEMBERS the
    // parse filled -- by construction exactly the set readBig serves.
    DynamicList<word> keys;
    DynamicList<const List<scalar>*> blocks;
    auto add = [&](const word& k, const List<scalar>& b)
    {
        if (b.size()) { keys.append(k); blocks.append(&b); }
    };
    add("sourcePV", sourcePV_);
    add("T", T_table_);
    add("psisTab", psis_table_);
    add("dhRef", dhRef_table_);
    forAllConstIter(HashTable<List<scalar>>, optTables_, iter)
    {
        add(iter.key(), iter());
    }
    forAllConstIter(HashTable<List<scalar>>, Y_tables_, iter)
    {
        add("Y_" + iter.key(), iter());
    }
    if (RGcoeff_.size())
    {
        const wordList& cn = tabulatedRealGasMixture::coeffNames();
        forAll(RGcoeff_, k) add("RG_" + cn[k], RGcoeff_[k]);
    }
    forAllConstIter(HashTable<List<scalar>>, Le_tables_, iter)
    {
        add("Le_" + iter.key(), iter());
    }

    // --- <src>.bin : header + raw blocks, written to .tmp then renamed
    {
        const fileName tmp(cacheSrc_ + ".bin.tmp");
        std::ofstream os(tmp.c_str(), std::ios::binary);
        if (!os)
        {
            WarningInFunction
                << "cannot write " << tmp << " -- table cache disabled "
                << "(directory not writable?)" << endl;
            return;
        }

        FGMBinHeader h;
        memcpy(h.magic, fgmBinMagic, 8);
        h.srcSize = int64_t(Foam::fileSize(cacheSrc_));
        h.srcMTime = int64_t(Foam::lastModified(cacheSrc_));
        h.nBlocks = int32_t(keys.size());
        os.write(reinterpret_cast<const char*>(&h), sizeof(h));

        forAll(keys, b)
        {
            const int32_t klen = int32_t(keys[b].size());
            os.write(reinterpret_cast<const char*>(&klen), sizeof(klen));
            os.write(keys[b].c_str(), klen);
            const int64_t n = int64_t(blocks[b]->size());
            os.write(reinterpret_cast<const char*>(&n), sizeof(n));
            os.write
            (
                reinterpret_cast<const char*>(blocks[b]->begin()),
                std::streamsize(n)*std::streamsize(sizeof(scalar))
            );
        }
        if (!os) { WarningInFunction << "short write on " << tmp << endl; return; }
        os.close();
        Foam::mv(tmp, cacheSrc_ + ".bin");
    }

    // --- <src>.stub : the dictionary minus the big blocks
    {
        wordHashSet big;
        forAll(keys, b) big.insert(keys[b]);

        dictionary stub;
        forAllConstIter(dictionary, *this, iter)
        {
            if (!big.found(iter().keyword()))
            {
                stub.add(iter());
            }
        }

        const fileName tmp(cacheSrc_ + ".stub.tmp");
        OFstream os(tmp);
        stub.write(os, false);
        if (!os.good())
        {
            WarningInFunction << "short write on " << tmp << endl;
            return;
        }
        Foam::mv(tmp, cacheSrc_ + ".stub");
    }

    Info<< "    binary table cache written: " << cacheSrc_.name()
        << ".bin/.stub (" << keys.size() << " blocks)" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FGMTable::FGMTable
(
    const fvMesh& mesh
)
:
    FGMTable(mesh, "fgmProperties")
{}


Foam::FGMTable::FGMTable
(
    const fvMesh& mesh,
    const word& dictName
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            // With a valid sidecar cache the 388 MB parse is skipped
            // entirely: the base reads nothing, and the body merges the
            // small .stub dictionary before any lookup runs.
            probeCache(resolvedDictPath(mesh, dictName))
              ? IOobject::NO_READ
              : IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    nChi_(1),
    hasChi_(false),
    useEnthalpy_(false),
    useDilution_(false),
    hOx_(0),
    hFuel_(0),
    chi_axis_(1, scalar(0))
{
    Info<< "\nFGM (FPV + beta-PDF) table initialisation [" << dictName
        << "]" << endl;

    // Sidecar cache: probe again (cheap -- two stats and a 24-byte
    // header) and, on a hit, merge the stub so every dictionary lookup
    // below and in the consumers behaves exactly as after a full parse.
    cacheSrc_ = resolvedDictPath(mesh, dictName);
    cacheActive_ = probeCache(cacheSrc_);
    if (cacheActive_)
    {
        IFstream is(cacheSrc_ + ".stub");
        merge(dictionary(is));
        loadCacheBin();
    }

    nZ_ = readLabel(lookup("nZ"));
    nGz_ = readLabel(lookup("nGz"));
    nC_ = readLabel(lookup("nC"));
    Z_axis_ = List<scalar>(lookup("Z"));
    gZ_axis_ = List<scalar>(lookup("gZ"));
    C_axis_ = List<scalar>(lookup("C"));

    // -------- optional chi (scalar-dissipation) axis --------
    if (found("nChi"))
    {
        nChi_ = readLabel(lookup("nChi"));
        if (nChi_ < 1)
        {
            FatalErrorInFunction
                << "nChi must be >= 1 (read " << nChi_ << ")"
                << exit(FatalError);
        }
        if (nChi_ >= 2 || found("chi"))
        {
            chi_axis_ = List<scalar>(lookup("chi"));
            if (chi_axis_.size() != nChi_)
            {
                FatalErrorInFunction
                    << "chi axis size " << chi_axis_.size()
                    << " != nChi = " << nChi_
                    << exit(FatalError);
            }
            hasChi_ = true;
        }
    }

    // -------- optional enthalpy-defect axis (non-adiabatic FPV, method b) --------
    // The 4th axis machinery is generic; here it carries the enthalpy defect
    // dh [J/kg] instead of chi. Detected by 'fourthAxis enthalpy;'. We reuse the
    // chi_axis_/nChi_ slots so interpolateTable() is unchanged, but flag
    // useEnthalpy_ so the solver forms dh = h - ((1-Z)hOx + Z hFuel) per cell.
    if (found("fourthAxis") && word(lookup("fourthAxis")) == "enthalpy")
    {
        nChi_ = readLabel(lookup("nH"));
        chi_axis_ = List<scalar>(lookup("enthalpy"));
        if (chi_axis_.size() != nChi_)
        {
            FatalErrorInFunction
                << "enthalpy axis size " << chi_axis_.size()
                << " != nH = " << nChi_ << exit(FatalError);
        }
        hOx_   = readScalar(lookup("hOx"));
        hFuel_ = readScalar(lookup("hFuel"));
        hasChi_ = true;          // use the 4-D interpolation path
        useEnthalpy_ = true;
        Info<< "    NON-ADIABATIC FPV: 4th axis = enthalpy defect, nH=" << nChi_
            << "  dh=[" << chi_axis_[0] << "," << chi_axis_[nChi_-1] << "] J/kg"
            << nl << "    hOx=" << hOx_ << " hFuel=" << hFuel_ << " J/kg" << endl;

        // Adiabatic-enthalpy offset (optional; see FGMTable.H). Without it the
        // solver keeps the legacy mixing-line reference, which is exact only
        // for a unity-Lewis manifold.
        if (foundBig("dhRef"))
        {
            dhRef_table_ = readBig("dhRef");
            const label nTot = nZ_*nGz_*nC_*nChi_;
            if (dhRef_table_.size() != nTot)
            {
                FatalErrorInFunction
                    << "dhRef size " << dhRef_table_.size() << " != "
                    << nTot << exit(FatalError);
            }
            scalar lo = GREAT, hi = -GREAT;
            forAll(dhRef_table_, i)
            {
                lo = min(lo, dhRef_table_[i]); hi = max(hi, dhRef_table_[i]);
            }
            Info<< "    dhRef (adiabatic offset from the mixing line): ["
                << lo/1e6 << ", " << hi/1e6 << "] MJ/kg -- dh measured from"
                << " the manifold's own adiabatic state" << endl;
        }
        else
        {
            Info<< "    dhRef absent: dh referenced to the mixing line"
                << " (exact only for unity-Lewis manifolds)" << endl;
        }
    }

    // -------- optional steam-dilution axis (H2/O2/H2O power-generation FPV) --
    // 4th axis W = steam mole fraction in the oxidiser stream. Reuses the same
    // chi_axis_/nChi_ 4-D machinery, but the coordinate is a TRANSPORTED
    // conserved dilution scalar W passed directly (no defect subtraction),
    // flagged by 'fourthAxis dilution;' with the axis under 'nW'/'W'.
    if (found("fourthAxis") && word(lookup("fourthAxis")) == "dilution")
    {
        nChi_ = readLabel(lookup("nW"));
        chi_axis_ = List<scalar>(lookup("W"));
        if (chi_axis_.size() != nChi_)
        {
            FatalErrorInFunction
                << "dilution axis size " << chi_axis_.size()
                << " != nW = " << nChi_ << exit(FatalError);
        }
        hasChi_ = true;          // use the 4-D interpolation path
        useDilution_ = true;
        Info<< "    STEAM-DILUTED FPV: 4th axis = oxidiser steam fraction W, "
            << "nW=" << nChi_
            << "  W=[" << chi_axis_[0] << "," << chi_axis_[nChi_-1] << "]"
            << endl;
    }

    Info<< "    axes: nZ=" << nZ_
        << "  nGz=" << nGz_
        << "  nC="  << nC_
        << "  nChi=" << nChi_
        << (hasChi_ ? "  (4-D table)" : "  (3-D legacy table)")
        << endl;

    if (nZ_ < 2 || nGz_ < 2 || nC_ < 2)
    {
        FatalErrorInFunction
            << "Each of Z/gZ/C axes must have at least 2 entries: "
            << "nZ=" << nZ_ << " nGz=" << nGz_ << " nC=" << nC_
            << exit(FatalError);
    }

    if
    (
        Z_axis_.size()   != nZ_
     || gZ_axis_.size()  != nGz_
     || C_axis_.size()   != nC_
    )
    {
        FatalErrorInFunction
            << "FGM axis lengths do not match declared sizes." << nl
            << "Z: "  << Z_axis_.size()  << " (nZ=" << nZ_ << ")" << nl
            << "gZ: " << gZ_axis_.size() << " (nGz=" << nGz_ << ")" << nl
            << "C: "  << C_axis_.size()  << " (nC=" << nC_ << ")"
            << exit(FatalError);
    }

    // -------- main PV source table --------
    // Yc_eq(Z) (optional): the normalisation behind the table's c = Yc/Yc_eq.
    // Needed only when the solver transports the UNNORMALISED Yc, which is the
    // formulation that has no normalisation cross-terms (see fgmFluid
    // transportYc). One value per Z axis node.
    if (found("Cnorm"))
    {
        Cnorm_table_ = List<scalar>(lookup("Cnorm"));
        if (Cnorm_table_.size() != nZ_)
        {
            FatalErrorInFunction
                << "Cnorm size " << Cnorm_table_.size() << " != nZ " << nZ_
                << exit(FatalError);
        }
        scalar lo = GREAT, hi = -GREAT;
        forAll(Cnorm_table_, i)
        {
            lo = min(lo, Cnorm_table_[i]); hi = max(hi, Cnorm_table_[i]);
        }
        Info<< "    Cnorm Yc_eq(Z): [" << lo << ", " << hi << "]" << endl;
    }

    // -------- optional per-node blocks (PEQSI coefficients, W) --------
    // Read if present; production fgmFluid tables carry none of these.
    {
        static const char* optNames[] =
        {
            "PEQSI_cv", "PEQSI_dpdT_v", "PEQSI_dpdv_T",
            "PEQSI_xi", "PEQSI_alpha", "PEQSI_beta",
            "PEQSI_dpdTn", "PEQSI_dpdvn", "PEQSI_betan",
            "W",
            // Chung/Takahashi transport as the tabulator baked it.
            // NOT a substitute for the runtime mu/kappa: measured against
            // the live mixture on the 2-D shear case they agree at the
            // gas end (1.49e-5 vs 1.56e-5 Pa.s) and diverge by two orders
            // at the cold dense end (0.284 vs 0.00133 Pa.s), where Z = 1
            // decane sits below its freezing point and the two Chung
            // implementations extrapolate differently.  Read for that
            // comparison, not for consumption.
            "mu", "kappa",
            // Table-state density.  Diagnostic in origin, but S_Y needs a
            // density response to a composition change and thermo::rho()
            // is a CACHED field -- it does not answer until correct()
            // runs -- so the manifold-direction derivative comes from
            // here instead.
            "PEQSI_rho"
        };
        const label nNode = nZ_*nGz_*nC_*max(nChi_, 1);

        for (const char* nm : optNames)
        {
            if (foundBig(nm))
            {
                List<scalar> tbl(readBig(nm));
                if (tbl.size() != nNode)
                {
                    FatalErrorInFunction
                        << nm << " size " << tbl.size()
                        << " != node count " << nNode << exit(FatalError);
                }
                optTables_.insert(word(nm), move(tbl));
            }
        }
        if (optTables_.size())
        {
            Info<< "    optional blocks: " << optTables_.toc() << endl;
        }
    }

    pRef_ = lookupOrDefault<scalar>("pressure", 0);

    sourcePV_ = readBig("sourcePV");
    const label nTot = nZ_*nGz_*nC_*nChi_;
    Info<< "    sourcePV entries: " << sourcePV_.size()
        << " (expected " << nTot << ")" << nl << endl;
    if (sourcePV_.size() != nTot)
    {
        FatalErrorInFunction
            << "sourcePV size " << sourcePV_.size()
            << " != nZ*nGz*nC*nChi = " << nTot
            << exit(FatalError);
    }

    // -------- optional temperature table --------
    if (foundBig("T"))
    {
        T_table_ = readBig("T");
        if (T_table_.size() != nTot)
        {
            FatalErrorInFunction
                << "T table size " << T_table_.size() << " != " << nTot
                << exit(FatalError);
        }
        Info<< "    tabulated: T" << endl;
    }

    // -------- optional pre-tabulated isentropic compressibility (psis) --------
    // Offline SRK+JANAF re-derivation of psis=psi/(rho*gamma) at the manifold's
    // reference pressure, smoothed at build time (gamma clipped to [1.05,2.5]
    // before the divide) as a principled alternative to the runtime
    // psisCapRatio neighbour-average patch -- see build_psis_table_v2.py and
    // Obsidian "PEP and LAD Spike Suppression" 2026-07-24. Default absent
    // (psisTabulated switch in pressureCorrector.C falls back if not found).
    if (foundBig("psisTab"))
    {
        psis_table_ = readBig("psisTab");
        if (psis_table_.size() != nTot)
        {
            FatalErrorInFunction
                << "psisTab table size " << psis_table_.size() << " != " << nTot
                << exit(FatalError);
        }
        Info<< "    tabulated: psisTab" << endl;
    }

    // -------- optional species composition tables --------
    if (found("species"))
    {
        speciesNames_ = wordList(lookup("species"));
        forAll(speciesNames_, i)
        {
            const word key("Y_" + speciesNames_[i]);
            List<scalar> tbl(readBig(key));
            if (tbl.size() != nTot)
            {
                FatalErrorInFunction
                    << key << " size " << tbl.size() << " != " << nTot
                    << exit(FatalError);
            }
            Y_tables_.insert(speciesNames_[i], tbl);
        }
        Info<< "    tabulated species: " << speciesNames_ << nl << endl;
    }

    // -------- optional Tier-2 real-gas mixture coefficient tables --------
    // Pre-tabulated output of SRKchungTakaMixture::calculateRealGas, keyed
    // RG_<name> in tabulatedRealGasMixture::coeffNames() order. All 13 must be
    // present (and full length) for the solver to enable the lookup; otherwise
    // the table is treated as a legacy (no-RG) table and live mixing is used.
    {
        const wordList& cn = tabulatedRealGasMixture::coeffNames();
        bool hasAll = true;
        forAll(cn, k)
        {
            if (!foundBig("RG_" + cn[k])) { hasAll = false; break; }
        }

        if (hasAll)
        {
            RGcoeff_.setSize(cn.size());
            forAll(cn, k)
            {
                const word key("RG_" + cn[k]);
                List<scalar> tbl(readBig(key));
                if (tbl.size() != nTot)
                {
                    FatalErrorInFunction
                        << key << " size " << tbl.size() << " != " << nTot
                        << exit(FatalError);
                }
                RGcoeff_[k] = tbl;
            }
            Info<< "    tabulated real-gas coefficients (Tier-2): "
                << cn << nl << endl;
        }
    }

    // -------- optional Tier-4 differential-diffusion Lewis-number tables --------
    // Le_<var> flat fields (length nTot) for the control variables. When
    // present the solver applies a per-cell Le(Z,gZ,c[,chi]) to that variable's
    // laminar diffusivity (rho*D = mu/Le), the standard differential-diffusion
    // FGM closure, generalising the constant 'Le' sub-dict. Each is independent:
    // a table may carry Le_Z, Le_C, both, or neither.
    {
        const wordList leVars({word("Z"), word("C")});
        forAll(leVars, k)
        {
            const word key("Le_" + leVars[k]);
            if (foundBig(key))
            {
                List<scalar> tbl(readBig(key));
                if (tbl.size() != nTot)
                {
                    FatalErrorInFunction
                        << key << " size " << tbl.size() << " != " << nTot
                        << exit(FatalError);
                }
                Le_tables_.insert(leVars[k], tbl);
            }
        }
        if (!Le_tables_.empty())
        {
            Info<< "    tabulated differential-diffusion Lewis numbers "
                << "(Tier-4): " << Le_tables_.sortedToc() << nl << endl;
        }
    }

    // 축 버킷. 모든 축이 확정된 뒤(4-D 는 chi/enthalpy/W 를 늦게 채운다) 1회.
    // 균등 축이면 버킷이 있어도 같은 셀을 찍으므로 결과는 비트 동일하다.
    buildBucket(Z_axis_,   Z_buck_);
    buildBucket(gZ_axis_,  gZ_buck_);
    buildBucket(C_axis_,   C_buck_);
    buildBucket(chi_axis_, chi_buck_);
    {
        auto ratio = [](const List<scalar>& a) -> scalar
        {
            if (a.size() < 3) return 1;
            scalar lo = GREAT, hi = 0;
            for (label i = 1; i < a.size(); i++)
            {
                const scalar d = a[i] - a[i-1];
                lo = min(lo, d); hi = max(hi, d);
            }
            return (lo > VSMALL) ? hi/lo : 1;
        };
        const scalar rZ = ratio(Z_axis_);
        if (rZ > 1.05)
        {
            Info<< "    non-uniform Z axis (dZ ratio " << rZ
                << "), bucket index " << Z_buck_.size() << " entries" << nl
                << endl;
        }
    }

    // A full parse just happened: leave the sidecars behind so the next
    // construction -- this case or any case symlinking the same table --
    // skips it.  A cache-hit construction leaves them untouched.
    if (!cacheActive_)
    {
        writeCache();
    }
    cacheBlocks_.clear();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FGMTable::~FGMTable()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::FGMTable::buildBucket(const List<scalar>& axis, List<label>& buck)
{
    buck.clear();
    const label n = axis.size();
    if (n < 2)
    {
        return;
    }
    const scalar span = axis[n - 1] - axis[0];
    if (span <= VSMALL)
    {
        return;
    }
    scalar dmin = GREAT;
    for (label i = 1; i < n; i++)
    {
        dmin = min(dmin, axis[i] - axis[i - 1]);
    }
    if (dmin <= VSMALL)
    {
        return;
    }
    // 축 셀당 버킷 1개 이상. 상한 2^20 (4 MB) 로 병적인 축을 막는다.
    label M = label(std::ceil(span/dmin));
    M = max(n, min(M, label(1 << 20)));
    buck.setSize(M + 1);
    label j = 0;
    for (label k = 0; k <= M; k++)
    {
        const scalar x = axis[0] + span*scalar(k)/scalar(M);
        while (j + 1 < n && axis[j + 1] <= x)
        {
            j++;
        }
        buck[k] = j;
    }
}


void Foam::FGMTable::bracket
(
    const List<scalar>& axis,
    const List<label>& buck,
    scalar v,
    label& i,
    scalar& w
) const
{
    const label n = axis.size();

    // Degenerate axis: collapse to that slice (used when nChi=1).
    if (n <= 1)
    {
        i = 0;
        w = 0;
        return;
    }

    if (v <= axis[0])
    {
        i = 0;
        w = 0;
        return;
    }
    if (v >= axis[n - 1])
    {
        i = n - 2;
        w = 1;
        return;
    }

    // Find j with axis[j-1] < v <= axis[j]. The manifold axes are uniform
    // (Z/gZ/C/dh are linspace), so an initial guess assuming constant
    // spacing lands on the exact cell (or a neighbour), and the two guarded
    // walks below run 0-1 times -> O(1). For a non-uniform axis the guess is
    // just a starting point and the walks recover the correct cell, so this
    // is a drop-in, bit-identical replacement for the old O(n) linear scan.
    const scalar span = axis[n - 1] - axis[0];
    label j = 1;
    const label M = buck.size() - 1;
    if (M >= 1 && span > VSMALL)
    {
        // 버킷: 비균등 축에서도 O(1) 로 시작 셀을 찍는다.
        label k = label((v - axis[0])/span*scalar(M));
        k = max(label(0), min(k, M));
        j = max(label(1), min(buck[k] + 1, n - 1));
    }
    else if (span > VSMALL)
    {
        j = label((v - axis[0])/span*(n - 1)) + 1;
        j = max(label(1), min(j, n - 1));
    }
    while (j < n - 1 && axis[j] < v)
    {
        j++;
    }
    while (j > 1 && axis[j - 1] >= v)
    {
        j--;
    }

    i = j - 1;
    const scalar d = axis[j] - axis[i];
    w = (d > VSMALL) ? (v - axis[i])/d : 0;
}


void Foam::FGMTable::makeStencil
(
    scalar Z,
    scalar gZ,
    scalar C,
    scalar chi,
    FGMStencil& st
) const
{
    label iZ, iG, iC, iK;

    bracket(Z_axis_, Z_buck_,   Z,   iZ, st.wZ);
    bracket(gZ_axis_, gZ_buck_,  gZ,  iG, st.wG);
    bracket(C_axis_, C_buck_,   C,   iC, st.wC);
    bracket(chi_axis_, chi_buck_, chi, iK, st.wK);

    // Degenerate-axis folding: for any axis of length 1 the bracket above
    // returns i=0, w=0; the +1 corner must fold onto the same slice or the
    // stencil walks into the neighbouring slice / past the allocation.
    // (Previously only the chi axis was folded — latent OOB for tables
    // collapsed in Z/gZ/C.  backport RGP-13-GPU 2005bc9)
    const label iZp = (Z_axis_.size()   >= 2) ? (iZ + 1) : iZ;
    const label iGp = (gZ_axis_.size()  >= 2) ? (iG + 1) : iG;
    const label iCp = (C_axis_.size()   >= 2) ? (iC + 1) : iC;
    const label iKp = (chi_axis_.size() >= 2) ? (iK + 1) : iK;

    // Corner order = interpolateTable's c0000..c1111 (Z fastest, then gZ,
    // then C; chi-low slice first).
    label n = 0;
    for (label dK = 0; dK < 2; dK++)
    {
        const label k = dK ? iKp : iK;
        for (label dC = 0; dC < 2; dC++)
        {
            const label c = dC ? iCp : iC;
            for (label dG = 0; dG < 2; dG++)
            {
                const label g = dG ? iGp : iG;
                for (label dZ = 0; dZ < 2; dZ++)
                {
                    const label z = dZ ? iZp : iZ;
                    st.idx[n++] = flatIndex(z, g, c, k);
                }
            }
        }
    }
}


Foam::scalar Foam::FGMTable::interpolate
(
    const List<scalar>& table,
    const FGMStencil& st
) const
{
    const scalar wZ = st.wZ, wG = st.wG, wC = st.wC, wK = st.wK;

    // Gather the 16 corners (same order as interpolateTable's c0000..c1111).
    const scalar c0000 = table[st.idx[0]];
    const scalar c1000 = table[st.idx[1]];
    const scalar c0100 = table[st.idx[2]];
    const scalar c1100 = table[st.idx[3]];
    const scalar c0010 = table[st.idx[4]];
    const scalar c1010 = table[st.idx[5]];
    const scalar c0110 = table[st.idx[6]];
    const scalar c1110 = table[st.idx[7]];
    const scalar c0001 = table[st.idx[8]];
    const scalar c1001 = table[st.idx[9]];
    const scalar c0101 = table[st.idx[10]];
    const scalar c1101 = table[st.idx[11]];
    const scalar c0011 = table[st.idx[12]];
    const scalar c1011 = table[st.idx[13]];
    const scalar c0111 = table[st.idx[14]];
    const scalar c1111 = table[st.idx[15]];

    // Nested blend VERBATIM from the original interpolateTable -- keeps the
    // floating-point evaluation order, hence bit-identical results.

    // Trilinear in (Z, gZ, C) at chi-low slice
    const scalar a00 = c0000*(1 - wZ) + c1000*wZ;
    const scalar a10 = c0100*(1 - wZ) + c1100*wZ;
    const scalar a01 = c0010*(1 - wZ) + c1010*wZ;
    const scalar a11 = c0110*(1 - wZ) + c1110*wZ;
    const scalar a0  = a00*(1 - wG) + a10*wG;
    const scalar a1  = a01*(1 - wG) + a11*wG;
    const scalar A   = a0 *(1 - wC) + a1 *wC;

    // Trilinear in (Z, gZ, C) at chi-high slice
    const scalar b00 = c0001*(1 - wZ) + c1001*wZ;
    const scalar b10 = c0101*(1 - wZ) + c1101*wZ;
    const scalar b01 = c0011*(1 - wZ) + c1011*wZ;
    const scalar b11 = c0111*(1 - wZ) + c1111*wZ;
    const scalar b0  = b00*(1 - wG) + b10*wG;
    const scalar b1  = b01*(1 - wG) + b11*wG;
    const scalar B   = b0 *(1 - wC) + b1 *wC;

    // Linear in chi
    return A*(1 - wK) + B*wK;
}


Foam::scalar Foam::FGMTable::interpolateTable
(
    const List<scalar>& table,
    scalar Z,
    scalar gZ,
    scalar C,
    scalar chi
) const
{
    // Single-field query: build the stencil and evaluate (the shared-stencil
    // fast path is makeStencil once + interpolate per field).
    FGMStencil st;
    makeStencil(Z, gZ, C, chi, st);
    return interpolate(table, st);
}


// -------- 4-D entry points --------

Foam::scalar Foam::FGMTable::interpolate
(
    scalar Z, scalar gZ, scalar C, scalar chi
) const
{
    return interpolateTable(sourcePV_, Z, gZ, C, chi);
}


Foam::scalar Foam::FGMTable::interpolateT
(
    scalar Z, scalar gZ, scalar C, scalar chi
) const
{
    if (T_table_.empty())
    {
        FatalErrorInFunction
            << "Temperature is not tabulated in fgmProperties."
            << exit(FatalError);
    }
    return interpolateTable(T_table_, Z, gZ, C, chi);
}


Foam::scalar Foam::FGMTable::interpolateY
(
    const word& specie,
    scalar Z, scalar gZ, scalar C, scalar chi
) const
{
    if (!Y_tables_.found(specie))
    {
        FatalErrorInFunction
            << "Species '" << specie
            << "' is not tabulated in fgmProperties." << nl
            << "Available: " << speciesNames_
            << exit(FatalError);
    }
    return interpolateTable(Y_tables_[specie], Z, gZ, C, chi);
}


Foam::scalar Foam::FGMTable::interpolatePsis
(
    scalar Z, scalar gZ, scalar C, scalar chi
) const
{
    if (psis_table_.empty())
    {
        FatalErrorInFunction
            << "psisTab is not tabulated in fgmProperties."
            << exit(FatalError);
    }
    return interpolateTable(psis_table_, Z, gZ, C, chi);
}


// -------- 3-D legacy entry points (evaluate at chi_axis_[0]) --------

Foam::scalar Foam::FGMTable::interpolate
(
    scalar Z, scalar gZ, scalar C
) const
{
    return interpolateTable(sourcePV_, Z, gZ, C, chi_axis_[0]);
}


Foam::scalar Foam::FGMTable::interpolateT
(
    scalar Z, scalar gZ, scalar C
) const
{
    if (T_table_.empty())
    {
        FatalErrorInFunction
            << "Temperature is not tabulated in fgmProperties."
            << exit(FatalError);
    }
    return interpolateTable(T_table_, Z, gZ, C, chi_axis_[0]);
}


Foam::scalar Foam::FGMTable::interpolateY
(
    const word& specie,
    scalar Z, scalar gZ, scalar C
) const
{
    if (!Y_tables_.found(specie))
    {
        FatalErrorInFunction
            << "Species '" << specie
            << "' is not tabulated in fgmProperties." << nl
            << "Available: " << speciesNames_
            << exit(FatalError);
    }
    return interpolateTable(Y_tables_[specie], Z, gZ, C, chi_axis_[0]);
}


Foam::scalar Foam::FGMTable::interpolatePsis
(
    scalar Z, scalar gZ, scalar C
) const
{
    if (psis_table_.empty())
    {
        FatalErrorInFunction
            << "psisTab is not tabulated in fgmProperties."
            << exit(FatalError);
    }
    return interpolateTable(psis_table_, Z, gZ, C, chi_axis_[0]);
}


// -------- Tier-2 real-gas coefficient interpolation --------

void Foam::FGMTable::interpolateRealGasCoeffs
(
    scalar Z, scalar gZ, scalar C, scalar chi,
    List<scalar>& coeffs
) const
{
    if (RGcoeff_.empty())
    {
        return;
    }

    coeffs.setSize(RGcoeff_.size());
    forAll(RGcoeff_, k)
    {
        coeffs[k] = interpolateTable(RGcoeff_[k], Z, gZ, C, chi);
    }
}


void Foam::FGMTable::cInterval(scalar C, scalar& cLo, scalar& cHi) const
{
    if (C_axis_.size() <= 1)
    {
        cLo = cHi = (C_axis_.size() == 1 ? C_axis_[0] : C);
        return;
    }

    label iC = 0;
    scalar wC = 0;
    bracket(C_axis_, C_buck_, C, iC, wC);

    cLo = C_axis_[iC];
    cHi = C_axis_[iC + 1];
}


// -------- Opt-1 base-blend node interpolation stencil --------

void Foam::FGMTable::interpStencil
(
    scalar Z, scalar gZ, scalar C, scalar chi,
    label nodes[16], scalar weights[16]
) const
{
    label iZ, iG, iC, iK;
    scalar wZ, wG, wC, wK;

    bracket(Z_axis_, Z_buck_,   Z,   iZ, wZ);
    bracket(gZ_axis_, gZ_buck_,  gZ,  iG, wG);
    bracket(C_axis_, C_buck_,   C,   iC, wC);
    bracket(chi_axis_, chi_buck_, chi, iK, wK);

    // Fold the chi-high neighbour to the same slice for a 3-D table (mirrors
    // interpolateTable), so the chi-high corners are valid indices with the
    // weight (1-wK)/wK distribution -- wK=0 here, so they contribute nothing.
    const label iKp = (chi_axis_.size() >= 2) ? (iK + 1) : iK;

    label m = 0;
    for (label dZ = 0; dZ <= 1; dZ++)
    {
        const scalar fZ = dZ ? wZ : (scalar(1) - wZ);
        for (label dG = 0; dG <= 1; dG++)
        {
            const scalar fG = dG ? wG : (scalar(1) - wG);
            for (label dC = 0; dC <= 1; dC++)
            {
                const scalar fC = dC ? wC : (scalar(1) - wC);
                for (label dK = 0; dK <= 1; dK++)
                {
                    const scalar fK = dK ? wK : (scalar(1) - wK);
                    nodes[m] =
                        flatIndex(iZ + dZ, iG + dG, iC + dC, dK ? iKp : iK);
                    weights[m] = fZ*fG*fC*fK;
                    m++;
                }
            }
        }
    }
}


// -------- Tier-4 differential-diffusion Lewis-number interpolation --------

Foam::scalar Foam::FGMTable::interpolateLe
(
    const word& var,
    scalar Z, scalar gZ, scalar C, scalar chi
) const
{
    if (!Le_tables_.found(var))
    {
        FatalErrorInFunction
            << "Lewis number Le_" << var
            << " is not tabulated in fgmProperties." << nl
            << "Available: " << Le_tables_.sortedToc()
            << exit(FatalError);
    }
    return interpolateTable(Le_tables_[var], Z, gZ, C, chi);
}


Foam::scalar Foam::FGMTable::interpolateDhRef
(
    scalar Z, scalar gZ, scalar C
) const
{
    if (dhRef_table_.empty())
    {
        FatalErrorInFunction
            << "dhRef is not tabulated in fgmProperties." << exit(FatalError);
    }
    // 4th-axis-replicated: any coordinate gives the same value.
    return interpolateTable(dhRef_table_, Z, gZ, C, chi_axis_[0]);
}


Foam::scalar Foam::FGMTable::interpolateCnorm(const scalar Z) const
{
    if (Cnorm_table_.empty())
    {
        FatalErrorInFunction
            << "Cnorm (Yc_eq(Z)) is not tabulated in fgmProperties; it is "
            << "required when transportYc is enabled." << exit(FatalError);
    }
    const scalar Zc = max(min(Z, Z_axis_[nZ_ - 1]), Z_axis_[0]);
    label i = 0;
    while (i < nZ_ - 2 && Z_axis_[i + 1] < Zc)
    {
        i++;
    }
    const scalar dZ = Z_axis_[i + 1] - Z_axis_[i];
    const scalar w = (dZ > VSMALL) ? (Zc - Z_axis_[i])/dZ : 0;
    return (1 - w)*Cnorm_table_[i] + w*Cnorm_table_[i + 1];
}


// ************************************************************************* //
