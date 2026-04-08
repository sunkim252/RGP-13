# FGM Combustion Model — Porting Guide & Code Review

## 1. Overview

FGM (Flamelet-Generated Manifold) combustion model ported from `references/DNS_FGM` (OpenFOAM 5~7, DNS-only standalone solver) to an OpenFOAM-13 `combustionModel` plugin.

**Key change:** Instead of a standalone solver (`fgmPremixedFoam`), FGM is now a plugin library loaded by the existing `multicomponentFluid` solver via runtime selection.

## 2. Architecture Comparison

| Aspect | Original (DNS_FGM) | OF-13 Port (RGP-13) |
|--------|-------------------|---------------------|
| Solver | `fgmPremixedFoam` (standalone) | `multicomponentFluid` (existing) |
| Thermo | `fgmThermo` (custom class) | Standard `hePsiThermo` + `multicomponentMixture` |
| PV transport | `PVEqn.H` (dedicated equation) | Species equation (PV as Yi) |
| Combustion | Embedded in solver | `combustionModel` plugin (`FGM`) |
| Table lookup | `lookupFGM` (part of solver) | `FGMTable` (part of library) |
| Turbulence | None (DNS only) | laminar / RAS / LES via `momentumTransport` |

## 3. New Files

### `src/FGM/FGMTable.H` / `FGMTable.C`
- Reads `constant/fgmProperties` dictionary
- Stores 4 lookup tables: PV axis, sourcePV, T, rho (each `List<scalar>`)
- `interpolate()`: 1D linear interpolation with PV clamped to [0,1]
- Ported from `lookupFGM.C` with edge case handling preserved

### `src/FGM/FGM.H` / `FGM.C`
- Inherits `combustionModel` (OF-13 base class)
- Registered as `TypeName("FGM")` in runtime selection table
- Constructor reads `progressVariable` name from `FGMCoeffs` (default: `"PV"`)
- Finds PV index in `thermo.species()` list

#### Implemented methods:
| Method | Implementation |
|--------|---------------|
| `correct()` | Loops over cells, interpolates `sourcePV` from FGM table using current PV field |
| `R(label speciei)` | Returns `sourcePV_` for PV species, zero for others |
| `R(volScalarField& Y)` | Builds `fvScalarMatrix` with explicit source for PV species |
| `Qdot()` | Returns zero (Step 1 — heat release via species enthalpy) |
| `read()` | Returns true (table is static) |

## 4. Porting Decisions

### 4.1 PV as a species
The original DNS_FGM treats PV as a standalone scalar field with its own transport equation. In OF-13, PV is registered as a species in the `multicomponentMixture` species list. This allows:
- Reuse of the existing species transport equation loop
- Automatic turbulent diffusion via `thermophysicalTransport`
- No solver modification needed

### 4.2 Qdot = 0
In Step 1, `Qdot()` returns zero. The energy equation still solves for enthalpy, but temperature evolution is driven by species transport (PV changes composition → thermo properties change). Full `Qdot` integration requires computing heat of combustion from the FGM table, planned for Step 2.

### 4.3 Table-driven density not enforced
The original DNS_FGM sets `rho` directly from the table. In OF-13, density comes from the EOS (`perfectGas` or `SRKGas`). The FGM table's `rho` column is read but not currently used. For stricter FGM coupling, a custom thermo model could override `rho()` with table values.

## 5. Case Setup

### `constant/combustionProperties`
```
combustionModel FGM;
FGMCoeffs { progressVariable PV; }
```

### `constant/physicalProperties`
```
thermoType
{
    type            hePsiThermo;
    mixture         multicomponentMixture;
    transport       sutherland;
    thermo          janaf;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleEnthalpy;
}
species (PV air);
// ... species definitions with JANAF coefficients
```

### `constant/momentumTransport`
```
// DNS
simulationType laminar;

// RAS (also works)
simulationType RAS;
RAS { model kEpsilon; turbulence on; }

// LES (also works)
simulationType LES;
LES { model Smagorinsky; turbulence on; delta cubeRootVol; }
```

### `system/controlDict`
```
solver multicomponentFluid;
libs ("libRGP13realFluid.so");
```

### `system/fvSchemes`
```
div(phi,Yi_h) Gauss limitedLinear 1;   // Required for multicomponentFluid
```

### `system/fvSolution`
```
Yi { solver PBiCGStab; preconditioner DILU; tolerance 1e-6; relTol 0.1; }
YiFinal { $Yi; relTol 0; }
```

## 6. Build & Test

```bash
openfoam13
cd /home/openfoam/RGP-13/RGP-13-realFluid
./Allwmake

cd /home/openfoam/RGP-13/test/bunsenFlame_FGM
blockMesh
foamRun
```

### Test result (3 timesteps, dt=1e-8):
```
Time = 3e-07s
DILUPBiCGStab: Solving for PV, Initial residual = 0.122, Final residual = 1.2e-11, No Iterations 1
continuity errors: sum local = 2.7e-12, global = -1.1e-12
ExecutionTime = 3.66 s
```

## 7. Code Review Checklist

### Correctness
- [ ] `FGMTable::interpolate()` — boundary handling at PV=0 and PV=1
- [ ] `FGM::correct()` — cell loop matches table indexing
- [ ] `FGM::R(label)` — correct species index comparison
- [ ] `FGM::R(volScalarField&)` — source term sign and dimensions
- [ ] Species name matching: `pvName_` vs `thermo.species()` lookup

### Dimensions
- [ ] `sourcePV_`: `[kg/m³/s]` = `dimDensity/dimTime`
- [ ] `Qdot_`: `[W/m³]` = `dimEnergy/dimVolume/dimTime`
- [ ] `R(label)` returns `volScalarField::Internal` (not `volScalarField`)

### Robustness
- [ ] Division by zero guard in `interpolate()` (smallValue = 1e-5)
- [ ] PV clamping to [0, 1] in interpolation
- [ ] Table size validation (minimum 2 entries)
- [ ] Species "PV" not found in thermo → `FatalError` with clear message

### Future improvements
- [ ] Step 2: 2D table (Z, c) for non-premixed / partially premixed
- [ ] Step 3: PV variance transport for RAS
- [ ] Full Qdot from table (hBurnt - hUnburnt) × sourcePV
- [ ] Table-driven density override option
- [ ] Integration with SRK/Chung real-fluid thermo
