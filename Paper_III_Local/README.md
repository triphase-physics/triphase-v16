# Paper III - Local Laboratory Evidence

Verification scripts for:
"Local Laboratory Evidence for Three-Phase Vacuum Structure"
C.R. Fuccillo, Magnetic Innovative Solutions LLC

## Scripts

| Script | Description | Requirements |
|--------|-------------|-------------|
| `x1_bessel_verification.py` | Derives Bessel first zero x1 = 2.4048 from vacuum properties (epsilon_0, mu_0, Z_0) via mode structure. Error: 0.12% | numpy |
| `fig_bessel_ratios.py` | Computes 9 Bessel J0 zero ratios, compares to mode integer predictions (all within 0.55%), generates interactive two-panel figure | numpy, scipy, plotly |

## Key Results

- **Bessel first zero from vacuum**: x1 = [epsilon_0 * mu_0 * 10^18 / (6+1)^2] * cbrt(pi * Z_0) = 2.4019 (exact: 2.4048, error 0.12%)
- **9 Bessel zero ratios match mode integers** to better than 0.55%:
  - x3/x1 = 18/5 (0.04%), x8/x2 = 22/5 (0.26%), x5/x1 = 137/22 (0.30%)
- **Walking droplet quantization** (Bush/MIT): same Bessel zeros govern mm-scale droplet orbits and vacuum mode structure across 21 orders of magnitude in energy

## Running

```
python x1_bessel_verification.py
python fig_bessel_ratios.py
```

Each script displays step-by-step calculations and pauses for review (Press Enter to exit).
