# Paper I: CMB Acoustic Peak Structure from Vacuum Electromagnetic Properties

Verification scripts for Paper I. Derives three CMB acoustic peak positions
from three vacuum electromagnetic properties:

| Peak | Vacuum Property | Derived | Observed (Planck 2018) |
|------|----------------|---------|----------------------|
| l1 | epsilon_0 | 220.02 | 220.0 +/- 0.5 |
| l2 | alpha^-1, mu_0 | 532.78 | 533.2 +/- 5.2 |
| l3 | Z_0 | 816.91 | 816.9 +/- 2.8 |

## Scripts

- `verify_cmb_l1_from_epsilon0.py` -- Derives l1 from vacuum permittivity
- `verify_cmb_l2_from_alpha_mu0.py` -- Derives l2 from fine structure constant and permeability
- `verify_cmb_l3_from_Z0.py` -- Derives l3 from vacuum impedance
- `verify_cmb_peak_ratio.py` -- Verifies l3/l2 = 3/2 = phases/quadratures
- `verify_111_projection.py` -- [111] crystallographic projection (90 deg -> 120 deg)
- `verify_bessel_first_zero.py` -- Bessel J0 first zero x1 = 2.4048
- `verify_bessel_mode_ratios.py` -- Bessel zero ratios match mode integer ratios
- `x1_bessel_verification.py` -- Extended Bessel zero verification
- `analyze_cmb_peaks.py` -- Gaussian peak fitting on Planck TT power spectrum
- `fig_bessel_ratios.py` -- Interactive Plotly chart of Bessel ratio comparisons
