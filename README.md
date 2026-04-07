# RGP-13: Real Fluid Models for OpenFOAM-13

Supercritical combustion modeling framework for OpenFOAM-13, implementing real-gas thermophysical models based on Oefelein's approach.

## Models

| Component | Model | Description |
|-----------|-------|-------------|
| **Equation of State** | SRK (Soave-Redlich-Kwong) | Cubic EOS with temperature-dependent attraction |
| **Viscosity & Thermal Conductivity** | Chung et al. (1988) | High-pressure transport properties |
| **Mass Diffusivity** | Takahashi (1975) + Fuller | Pressure-corrected binary diffusion |
| **Mixing Rules** | Quadratic combining rules | Pre-computed interaction matrices |

## Installation

```bash
# Inside your OpenFOAM-13 environment:
./install.sh /path/to/OpenFOAM-13

# Or if FOAM_SRC is already set:
./install.sh

# Preview changes without applying:
./install.sh --dry-run

# Remove installed files:
./install.sh --uninstall
```

## Usage

Set `thermoType` in `constant/physicalProperties`:

```
thermoType
{
    type            hePsiThermo;
    mixture         SRKchungTakaMixture;
    transport       chung;
    thermo          janaf;
    equationOfState SRKGas;
    specie          rfSpecie;
    energy          sensibleEnthalpy;
}
```

Each species requires an `rfProperties` sub-dictionary with critical properties:

```
O2
{
    specie      { molWeight 31.998; }
    rfProperties
    {
        Tc      154.58;     // Critical temperature [K]
        Pc      5.043e6;    // Critical pressure [Pa]
        Vc      73.37;      // Critical volume [cm3/mol]
        omega   0.0222;     // Acentric factor
        miui    0.0;        // Dipole moment [Debye]
        kappai  0.0;        // Association factor
        sigmvi  16.6;       // Fuller diffusion volume
    }
    thermodynamics { /* JANAF coefficients */ }
    transport {}
}
```

## File Structure

```
src/
├── rfSpecie/              Base specie class with critical properties
├── SRKGas/                Soave-Redlich-Kwong equation of state
├── chungTransport/        Chung viscosity/conductivity + Takahashi diffusion
├── SRKchungTakaMixture/   Quadratic mixing rules for multicomponent
└── include/               Template instantiation macros
```

## References

- Oefelein, J.C. (2006). "Mixing and combustion of cryogenic oxygen-hydrogen shear-coaxial jet flames at supercritical pressure"
- Chung, T.H. et al. (1988). "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties"
- Takahashi, S. (1975). "Preparation of a Generalized Chart for the Diffusion Coefficients of Gases at High Pressures"
- Soave, G. (1972). "Equilibrium constants from a modified Redlich-Kwong equation of state"

## Acknowledgments

Based on [realFluidFoam-8](https://github.com/) by Combustion & Propulsion Lab, UNIST (Prof. Chun Sang Yoo).
Ported to OpenFOAM-13 with bug fixes, API adaptation, and performance optimizations.
