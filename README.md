# TriPhase V16 - Peer Review Verification Scripts

Companion verification scripts for the TriPhase paper series (Papers I-IV).
Every calculation in all four papers can be reproduced by running these scripts.

C.R. Fuccillo, Magnetic Innovative Solutions LLC

## Papers

| Paper | Title | Directory | Scripts |
|-------|-------|-----------|---------|
| I | CMB Acoustic Peak Structure from Vacuum Electromagnetic Properties | `Paper_I_CMB/` | 10 |
| II | Dark Energy Equation of State from Vacuum Mode Structure | `Paper_II_w0/` | 3 |
| III | Local Laboratory Evidence for Vacuum Mode Structure | `Paper_III_Local/` | -- |
| IV | Hubble Tension Resolution via Three-Phase Vacuum Structure | `Paper_IV_Hubble/` | 3 |

## Quick Start

```bash
# Run master verification (all 8 checks)
python verify_calculations.py

# Run individual verifications
python Paper_I_CMB/verify_cmb_l1_from_epsilon0.py
python Paper_II_w0/verify_w0_from_mode_structure.py
python Paper_IV_Hubble/robustness_analysis.py
```

## Structure

```
triphase-v16/
├── verify_calculations.py          # Master verification (8 checks, all papers)
├── Paper_I_CMB/                    # CMB acoustic peak verification
│   ├── verify_cmb_l1_from_epsilon0.py      # l1 = 220.02 from epsilon_0
│   ├── verify_cmb_l2_from_alpha_mu0.py     # l2 = 532.78 from alpha^-1, mu_0
│   ├── verify_cmb_l3_from_Z0.py            # l3 = 816.91 from Z_0
│   ├── verify_cmb_peak_ratio.py            # l3/l2 = 3/2 (phases/quadratures)
│   ├── verify_111_projection.py            # [111] crystallographic projection
│   ├── verify_bessel_first_zero.py         # x1 = 2.4048 verification
│   ├── verify_bessel_mode_ratios.py        # Bessel zero ratios vs mode integers
│   ├── x1_bessel_verification.py           # Extended Bessel zero analysis
│   ├── analyze_cmb_peaks.py                # Planck TT power spectrum peak fitting
│   └── fig_bessel_ratios.py                # Interactive Bessel ratio chart (Plotly)
├── Paper_II_w0/                    # Dark energy equation of state
│   ├── verify_w0_from_mode_structure.py    # w0 = -181/216 from mode counting
│   ├── w0_balance_equation_verification.py # Balance equation verification
│   └── fig_w0_comparison.py                # w0 measurements comparison (Plotly)
├── Paper_III_Local/                # Local laboratory evidence
├── Paper_IV_Hubble/                # Hubble tension analysis
│   ├── fig_hubble_tension.py               # H0 measurement timeline (Plotly)
│   ├── robustness_analysis.py              # Systematic robustness checks
│   └── monte_carlo_slope_amplification.py  # Monte Carlo slope analysis
├── Data/                           # Raw observational data
│   ├── COM_PowerSpect_CMB-TT-binned_R3.01.txt  # Planck 2018 TT (binned)
│   └── COM_PowerSpect_CMB-TT-full_R3.01.txt    # Planck 2018 TT (full)
└── Common/                         # Shared utilities
```

## Key Results

| Quantity | Derived | Observed | Error |
|----------|---------|----------|-------|
| l1 (CMB peak 1) | 220.02 | 220.0 +/- 0.5 | 0.01% |
| l2 (CMB peak 2) | 532.78 | 533.2 +/- 5.2 | 0.08% |
| l3 (CMB peak 3) | 816.91 | 816.9 +/- 2.8 | 0.001% |
| l3/l2 ratio | 1.500001 | 1.5000 | 0.00008% |
| w0 (dark energy) | -0.83796 | -0.838 +/- 0.055 | 0.005% |

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
