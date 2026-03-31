# TriPhase V16 - Verification and Derivative Scripts

All verification and derivative scripts for the TriPhase V16 framework release.

- **Paper verification scripts** (Papers I-IV): reproduce every calculation in the peer review papers
- **Master derivative chain scripts** (864 scripts): independently verify every row of the V16 Master Matrix across 18 mathematical frameworks

All scripts derive from three observed vacuum properties: epsilon_0, mu_0, Z_0.

C.R. Fuccillo, Magnetic Innovative Solutions LLC

---

## Paper Verification Scripts

| Paper | Title | Directory | Scripts |
|-------|-------|-----------|---------|
| I | CMB Acoustic Peak Structure from Vacuum Electromagnetic Properties | `Paper_I_CMB/` | 10 |
| II | Dark Energy Equation of State from Vacuum Mode Structure | `Paper_II_w0/` | 3 |
| III | Local Laboratory Evidence for Vacuum Mode Structure | `Paper_III_Local/` | -- |
| IV | Hubble Tension Resolution via Three-Phase Vacuum Structure | `Paper_IV_Hubble/` | 3 |

### Quick Start (Papers)

```bash
# Run master verification (all 8 checks)
python verify_calculations.py

# Run individual verifications
python Paper_I_CMB/verify_cmb_l1_from_epsilon0.py
python Paper_II_w0/verify_w0_from_mode_structure.py
python Paper_IV_Hubble/robustness_analysis.py
```

### Key Results (Papers)

| Quantity | Derived | Observed | Error |
|----------|---------|----------|-------|
| l1 (CMB peak 1) | 220.02 | 220.0 +/- 0.5 | 0.01% |
| l2 (CMB peak 2) | 532.78 | 533.2 +/- 5.2 | 0.08% |
| l3 (CMB peak 3) | 816.91 | 816.9 +/- 2.8 | 0.001% |
| l3/l2 ratio | 1.500001 | 1.5000 | 0.00008% |
| w0 (dark energy) | -0.83796 | -0.838 +/- 0.055 | 0.005% |

---

## Master Derivative Chain Scripts

864 standalone Python scripts (48 targets x 18 mathematical frameworks) corresponding to the V16 Master Derivative Chain (S3). Each script independently derives a physical quantity from epsilon_0, mu_0, Z_0 and compares to the CODATA/PDG reference value.

### Derivative Targets (48)

| Target | Folder |
|--------|--------|
| Fine-structure constant inverse | `Derivatives/alpha_inverse/` |
| Bottom quark mass | `Derivatives/bottom_quark_mass/` |
| Speed of light squared | `Derivatives/c_squared/` |
| Charm quark mass | `Derivatives/charm_quark_mass/` |
| Coulomb pressure | `Derivatives/coulomb_pressure/` |
| Critical density | `Derivatives/critical_density/` |
| Dark energy pressure | `Derivatives/dark_energy_pressure/` |
| Dark energy scale | `Derivatives/dark_energy_scale/` |
| Dark energy w0 | `Derivatives/dark_energy_w0/` |
| Down quark mass | `Derivatives/down_quark_mass/` |
| Einstein field equation | `Derivatives/einstein_field_equation/` |
| Electromagnetic pressure | `Derivatives/electromagnetic_pressure/` |
| Electron g-2 | `Derivatives/electron_g2/` |
| Electron mass | `Derivatives/electron_mass/` |
| Elementary charge | `Derivatives/elementary_charge/` |
| Energy per mode | `Derivatives/energy_per_mode/` |
| Gravitational constant | `Derivatives/gravitational_constant/` |
| Gravity-pressure slope | `Derivatives/gravity_pressure_slope/` |
| Reduced Planck constant | `Derivatives/hbar/` |
| Higgs mass | `Derivatives/Higgs_mass/` |
| Horizon 18-step | `Derivatives/horizon_18step/` |
| Hubble constant | `Derivatives/hubble_constant/` |
| Hydrostatic pressure | `Derivatives/hydrostatic_pressure/` |
| 3.5 keV line | `Derivatives/keV_3p5_line/` |
| Lyman-alpha wavelength | `Derivatives/lyman_alpha/` |
| Magneton ratio | `Derivatives/magneton_ratio/` |
| Matter density | `Derivatives/matter_density/` |
| MOND acceleration | `Derivatives/MOND_acceleration/` |
| Muon mass | `Derivatives/muon_mass/` |
| Neutrino mass | `Derivatives/neutrino_mass/` |
| Neutron mass | `Derivatives/neutron_mass/` |
| Photon deceleration | `Derivatives/photon_deceleration/` |
| Proton/electron mass ratio | `Derivatives/proton_electron_ratio/` |
| Proton mass | `Derivatives/proton_mass/` |
| Proton radius | `Derivatives/proton_radius/` |
| Rydberg constant | `Derivatives/rydberg_constant/` |
| Speed of light | `Derivatives/speed_of_light/` |
| Strange quark mass | `Derivatives/strange_quark_mass/` |
| Tau mass | `Derivatives/tau_mass/` |
| Thermal pressure | `Derivatives/thermal_pressure/` |
| Top quark mass | `Derivatives/top_quark_mass/` |
| Triangular T(17) | `Derivatives/triangular_T17/` |
| Up quark mass | `Derivatives/up_quark_mass/` |
| Vacuum rigidity | `Derivatives/vacuum_rigidity/` |
| Vector frame | `Derivatives/vector_frame/` |
| Velocity spacing | `Derivatives/velocity_spacing/` |
| W boson mass | `Derivatives/W_boson_mass/` |
| Z boson mass | `Derivatives/Z_boson_mass/` |

### 18 Mathematical Frameworks

Each target folder contains 18 scripts, one per framework:

1. Anchor Primitive
2. Dimensional Analysis
3. Energy Density
4. Error Correcting Code
5. Fourier Harmonic
6. Group Theory
7. Information Theory
8. Lattice Gauge
9. Number Theory
10. Periodic
11. QFT
12. Renormalization Group
13. Statistical Mechanics
14. Symmetry Breaking
15. Thermodynamics
16. Topological
17. WaveMechanics Coupling
18. WaveMechanics Primitive

### Quick Start (Derivatives)

```bash
# Run any individual derivative
python Derivatives/gravitational_constant/gravitational_constant_Anchor_Primitive.py
python Derivatives/alpha_inverse/alpha_inverse_WaveMechanics_Primitive.py
python Derivatives/dark_energy_w0/dark_energy_w0_Thermodynamics.py
```

Each script prints the derived value, reference value, and percent error.

---

## Repository Structure

```
triphase-v16/
├── verify_calculations.py          # Master verification (8 checks, all papers)
├── Paper_I_CMB/                    # CMB acoustic peak verification (10 scripts)
├── Paper_II_w0/                    # Dark energy equation of state (3 scripts)
├── Paper_III_Local/                # Local laboratory evidence
├── Paper_IV_Hubble/                # Hubble tension analysis (3 scripts)
├── Data/                           # Raw observational data (Planck 2018)
├── Common/                         # Shared utilities
├── Derivatives/                    # Master derivative chain (864 scripts)
│   ├── alpha_inverse/              # 18 framework scripts
│   ├── gravitational_constant/     # 18 framework scripts
│   ├── ...                         # (48 target folders total)
│   └── Z_boson_mass/               # 18 framework scripts
└── Docs/                           # Documentation
```

## Requirements

- Python 3.8+
- NumPy
- SciPy
- Plotly (for interactive charts -- opens in browser)

## Charts

The three figure scripts (`fig_*.py`) generate interactive HTML charts using Plotly
that open automatically in your default web browser. No GUI framework required.

## License

MIT License

## Author

Christian R. Fuccillo
MIS Magnetic Innovative Solutions LLC
