# Paper V - SI Derivation of Distance and Gravitational Coupling

Verification scripts for:
"SI Derivation of Cosmic Distance and Gravitational Coupling from Vacuum Electromagnetic Properties"
C.R. Fuccillo, Magnetic Innovative Solutions LLC

## Scripts

| Script | Description | Requirements |
|--------|-------------|-------------|
| `DISTANCE_CONVERGENCE_GRAPH.py` | 4-panel figure comparing standard (d = cz/H0) vs logarithmic (d = (c/H0)ln(1+z)) distance formulas. Shows standard diverges at high z while TriPhase converges. Includes megamaser data overlay. | numpy, matplotlib |
| `generate_individual_graphs.py` | Generates 3 individual publication figures: (1) smoking-gun divergence plot to z=10, (2) maser distance comparison at low z, (3) correction factor bar chart by redshift | numpy, matplotlib |
| `G_units_derivation.py` | SI unit dimensional analysis showing G = epsilon_0 * c^6 / I_P^2, explaining how F/m units become m^3/(kg*s^2) and why the numerical factor 7.5 = 15/2 from mode counting (0.5% error) | numpy |
| `geometric_vs_redshift_comparison.py` | Comprehensive comparison of geometric distances (megamasers, eclipsing binaries, lensing) to redshift-derived distances under standard and TriPhase cosmology | numpy, scipy |

## Key Results

- **Logarithmic distance formula**: d = (c/H0) * ln(1+z) gives geometric = wavelength distance at ALL redshifts (0% divergence)
- **Standard formula divergence**: d = cz/H0 diverges from geometric distance by 9% at z=0.1, 51% at z=0.5, 107% at z=1.0
- **G from vacuum properties**: G = epsilon_0 * c^6 / I_P^2 = (15/2) * epsilon_0 in SI, connecting gravitational and electromagnetic coupling through the Planck current
- **Mode counting origin**: The factor 15/2 = 7.5 traces to mode counting (5 active modes * 3 phases / 2 quadratures), error 0.5%

## Running

```
python DISTANCE_CONVERGENCE_GRAPH.py
python generate_individual_graphs.py
python G_units_derivation.py
python geometric_vs_redshift_comparison.py
```

Each script displays step-by-step calculations and pauses for review (Press Enter to exit).
Graph scripts save figures to the current directory and display them interactively.
