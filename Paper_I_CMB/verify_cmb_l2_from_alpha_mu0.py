"""
Paper I Calculation: CMB Second Acoustic Peak from alpha^-1 and mu_0
====================================================================
Companion script for:
  "CMB Acoustic Peak Structure from Vacuum Electromagnetic Properties"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Derives the second CMB acoustic peak position l2 = 532.84 from the fine
structure constant inverse alpha^-1 and vacuum permeability mu_0.

Requires: numpy, scipy
"""
import sys
import traceback


def main():
    import numpy as np
    from scipy.special import jn_zeros

    # ============================================================
    # MECHANISM (What this script does and why)
    # ============================================================
    print("=" * 72)
    print("  CMB SECOND ACOUSTIC PEAK FROM alpha^-1 AND mu_0 -- Paper I")
    print("  Magnetic Innovative Solutions LLC")
    print("=" * 72)
    print()
    print("WHAT THIS CALCULATES:")
    print("-" * 72)
    print("""
The second CMB acoustic peak corresponds to the first rarefaction
harmonic of the baryon-photon fluid.  Planck 2018 Gaussian-fitted
position: l2 = 533.2 +/- 5.2.

This derivation uses two vacuum properties:
  alpha^-1 = 137.035999177  (fine structure constant inverse, CODATA 2022)
  mu_0 = 1.25663706212e-6 H/m  (vacuum permeability, CODATA 2022)

and the first Bessel J_0 zero:
  x_1 = 2.4048255577...

The formula (two stages):
  Stage 1 -- Quantization:
    l2_base = pi * alpha^-1 * mu_0 * 10^6 * [alpha^-1 / (alpha^-1 + 2*x1/19)]
            = 540.00

  Stage 2 -- Mode coupling (rarefaction correction):
    l2 = l2_base * [1 - x1/180]
       = 540.00 * 0.98664 = 532.84

Where:
  pi            -- geometric factor from angular projection
  alpha^-1      -- fine structure constant inverse (137.036)
  mu_0 * 10^6   -- permeability scaled to CMB units
  19 = 17 + 2   -- mode hierarchy integer + quadratures
  x1/180        -- rarefaction correction (Bessel zero / mode product)
  180 = 5 x 36  -- active modes x 6^2
""")

    # ============================================================
    # STEP 1: Input constants
    # ============================================================
    print("=" * 72)
    print("  STEP 1: Input Constants")
    print("=" * 72)
    print()

    alpha_inv = 137.035999177
    mu_0 = 1.25663706212e-6  # H/m
    x1 = jn_zeros(0, 1)[0]

    print(f"  Fine structure constant inverse (CODATA 2022):")
    print(f"    alpha^-1 = {alpha_inv:.9f}")
    print()
    print(f"  Vacuum permeability (CODATA 2022):")
    print(f"    mu_0 = {mu_0:.11e} H/m")
    print()
    print(f"  First zero of Bessel J_0 (from scipy):")
    print(f"    x_1 = {x1:.10f}")
    print()
    print(f"  Mode integers:")
    print(f"    17       (mode hierarchy integer)")
    print(f"    2        (quadratures per phase)")
    print(f"    19 = 17 + 2")
    print()

    l2_observed = 533.2
    l2_obs_err = 5.2
    print(f"  Observed value (Planck 2018, Gaussian fit):")
    print(f"    l2 = {l2_observed:.1f} +/- {l2_obs_err:.1f}")
    print()

    # ============================================================
    # STEP 2: Step-by-step calculation
    # ============================================================
    print("=" * 72)
    print("  STEP 2: Step-by-Step Calculation")
    print("=" * 72)
    print()
    print(f"  Full formula:")
    print(f"    l2 = pi * alpha^-1 * mu_0 * 10^6 * [alpha^-1 / (alpha^-1 + 2*x1/19)] * [1 - x1/180]")
    print()

    # Part A: Raw product
    raw = np.pi * alpha_inv * mu_0 * 1e6
    print(f"  Part A: Raw product")
    print(f"    pi * alpha^-1 * mu_0 * 10^6")
    print(f"    = {np.pi:.10f} * {alpha_inv:.9f} * {mu_0:.11e} * 10^6")
    print(f"    = {np.pi:.10f} * {alpha_inv:.9f} * {mu_0 * 1e6:.11f}")
    print(f"    = {raw:.10f}")
    print()

    # Part B: Quantization correction (Bessel)
    two_x1_over_19 = 2 * x1 / 19
    denom = alpha_inv + two_x1_over_19
    quant_corr = alpha_inv / denom

    print(f"  Part B: Quantization correction (Bessel)")
    print(f"    2 * x1 / 19 = 2 * {x1:.10f} / 19")
    print(f"                 = {2*x1:.10f} / 19")
    print(f"                 = {two_x1_over_19:.10f}")
    print()
    print(f"    alpha^-1 + 2*x1/19 = {alpha_inv:.9f} + {two_x1_over_19:.10f}")
    print(f"                        = {denom:.10f}")
    print()
    print(f"    Quantization = {alpha_inv:.9f} / {denom:.10f}")
    print(f"                 = {quant_corr:.10f}")
    print()

    l2_base = raw * quant_corr
    print(f"  Stage 1 result (base):")
    print(f"    l2_base = {raw:.10f} * {quant_corr:.10f}")
    print(f"            = {l2_base:.6f}")
    print()

    # Part C: Mode coupling correction (rarefaction)
    rarefaction = 1 - x1 / 180
    print(f"  Part C: Mode coupling correction (rarefaction)")
    print(f"    x1 / 180 = {x1:.10f} / 180")
    print(f"             = {x1/180:.10f}")
    print(f"    [1 - x1/180] = 1 - {x1/180:.10f}")
    print(f"                 = {rarefaction:.10f}")
    print()
    print(f"    Physical origin: 180 = 5 x 36 = active_modes x 6^2")
    print(f"    The second peak is a rarefaction (underdensity), so the")
    print(f"    correction shifts it DOWN from the base quantization value.")
    print()

    # Part D: Final result
    l2_derived = l2_base * rarefaction
    print(f"  Part D: Final result")
    print(f"    l2 = {l2_base:.6f} * {rarefaction:.10f}")
    print(f"       = {l2_derived:.6f}")
    print()

    # ============================================================
    # STEP 3: Compare to observation
    # ============================================================
    print("=" * 72)
    print("  STEP 3: Comparison to Planck 2018")
    print("=" * 72)
    print()

    error_abs = abs(l2_derived - l2_observed)
    error_pct = error_abs / l2_observed * 100
    sigma_away = error_abs / l2_obs_err

    print(f"  Derived:  l2 = {l2_derived:.4f}")
    print(f"  Observed: l2 = {l2_observed:.1f} +/- {l2_obs_err:.1f}  (Planck 2018, Gaussian fit)")
    print()
    print(f"  Absolute error: |{l2_derived:.4f} - {l2_observed:.1f}| = {error_abs:.4f}")
    print(f"  Percent error:  {error_abs:.4f} / {l2_observed:.1f} * 100 = {error_pct:.4f}%")
    print(f"  Sigma:          {error_abs:.4f} / {l2_obs_err:.1f} = {sigma_away:.2f} sigma")
    print()

    passed = error_pct < 1.0
    status = "PASSED" if passed else "FAILED"
    print(f"  Result: {status} (threshold: < 1% error)")
    print()

    print("=" * 72)
    print(f"  COMPLETE -- l2 = {l2_derived:.2f} from alpha^-1 and mu_0, error = {error_pct:.4f}%")
    print("=" * 72)
    print()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
