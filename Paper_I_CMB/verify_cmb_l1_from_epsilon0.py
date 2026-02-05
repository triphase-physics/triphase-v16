"""
Paper I Calculation: CMB First Acoustic Peak from Vacuum Permittivity
=====================================================================
Companion script for:
  "CMB Acoustic Peak Structure from Vacuum Electromagnetic Properties"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Derives the first CMB acoustic peak position l1 = 220 from vacuum
permittivity epsilon_0 and the first Bessel J0 zero.

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
    print("  CMB FIRST ACOUSTIC PEAK FROM VACUUM PERMITTIVITY -- Paper I")
    print("  Magnetic Innovative Solutions LLC")
    print("=" * 72)
    print()
    print("WHAT THIS CALCULATES:")
    print("-" * 72)
    print("""
The cosmic microwave background (CMB) has acoustic peaks -- bright
spots at specific angular scales caused by sound waves in the early
universe.  The first peak at multipole l = 220 is one of the most
precisely measured numbers in cosmology (Planck 2018: 220.0 +/- 0.5).

Standard cosmology (LCDM) derives l1 from 6 fitted parameters.
This derivation uses a single measured property of empty space:

  epsilon_0 = 8.8541878128 x 10^-12 F/m  (vacuum permittivity)

combined with the first zero of the Bessel function J_0:

  x_1 = 2.4048255577...  (a mathematical constant like pi)

The formula:
  l1 = 5^2 * epsilon_0 * 10^12 * [22 / (22 + x1/18)]

Where:
  5^2 = 25    -- mode integer squared (coupling modes)
  10^12       -- scale factor (vacuum scale to CMB scale)
  22          -- total mode-coupling integer (Paper I, Table 1)
  18 = 3 x 6 -- phases x total modes
  x1/18       -- Bessel correction from cylindrical mode geometry

This gives l1 from vacuum properties with zero free parameters.
""")

    # ============================================================
    # STEP 1: Input constants
    # ============================================================
    print("=" * 72)
    print("  STEP 1: Input Constants")
    print("=" * 72)
    print()

    # Vacuum permittivity (CODATA 2022)
    epsilon_0 = 8.8541878128e-12  # F/m
    print(f"  Vacuum permittivity (CODATA 2022):")
    print(f"    epsilon_0 = {epsilon_0:.10e} F/m")
    print()

    # First Bessel J0 zero
    x1 = jn_zeros(0, 1)[0]
    print(f"  First zero of Bessel J_0 (from scipy):")
    print(f"    x_1 = {x1:.10f}")
    print()

    # Mode integers
    print(f"  Mode integers from three-phase structure:")
    print(f"    N_phases = 3")
    print(f"    N_total  = 6   (3 phases x 2 quadratures)")
    print(f"    18 = 3 x 6     (phases x total modes)")
    print(f"    22              (mode-coupling integer, Paper I Table 1)")
    print(f"    25 = 5^2        (coupling mode integer squared)")
    print()

    # Observed value
    l1_observed = 220.0
    l1_obs_err = 0.5
    print(f"  Observed value (Planck 2018):")
    print(f"    l1 = {l1_observed:.1f} +/- {l1_obs_err:.1f}")
    print()

    # ============================================================
    # STEP 2: Build the formula piece by piece
    # ============================================================
    print("=" * 72)
    print("  STEP 2: Step-by-Step Calculation")
    print("=" * 72)
    print()
    print(f"  Formula: l1 = 5^2 * epsilon_0 * 10^12 * [22 / (22 + x1/18)]")
    print()

    # Part A: Raw product
    raw = 25 * epsilon_0 * 1e12
    print(f"  Part A: Raw product")
    print(f"    5^2 * epsilon_0 * 10^12")
    print(f"    = 25 * {epsilon_0:.10e} * 10^12")
    print(f"    = 25 * {epsilon_0 * 1e12:.10f}")
    print(f"    = {raw:.10f}")
    print()

    # Part B: Bessel correction numerator/denominator
    x1_over_18 = x1 / 18
    denom = 22 + x1_over_18
    correction = 22 / denom
    print(f"  Part B: Bessel correction factor")
    print(f"    x1 / 18 = {x1:.10f} / 18")
    print(f"            = {x1_over_18:.10f}")
    print()
    print(f"    22 + x1/18 = 22 + {x1_over_18:.10f}")
    print(f"               = {denom:.10f}")
    print()
    print(f"    Correction = 22 / {denom:.10f}")
    print(f"               = {correction:.10f}")
    print()

    # Part C: Final result
    l1_derived = raw * correction
    print(f"  Part C: Final result")
    print(f"    l1 = {raw:.10f} * {correction:.10f}")
    print(f"       = {l1_derived:.6f}")
    print()

    # ============================================================
    # STEP 3: Compare to observation
    # ============================================================
    print("=" * 72)
    print("  STEP 3: Comparison to Planck 2018")
    print("=" * 72)
    print()

    error_abs = abs(l1_derived - l1_observed)
    error_pct = error_abs / l1_observed * 100
    sigma_away = error_abs / l1_obs_err

    print(f"  Derived:  l1 = {l1_derived:.4f}")
    print(f"  Observed: l1 = {l1_observed:.1f} +/- {l1_obs_err:.1f}")
    print()
    print(f"  Absolute error: |{l1_derived:.4f} - {l1_observed:.1f}| = {error_abs:.4f}")
    print(f"  Percent error:  {error_abs:.4f} / {l1_observed:.1f} * 100 = {error_pct:.4f}%")
    print(f"  Sigma:          {error_abs:.4f} / {l1_obs_err:.1f} = {sigma_away:.2f} sigma")
    print()

    passed = error_pct < 1.0
    status = "PASSED" if passed else "FAILED"
    print(f"  Result: {status} (threshold: < 1% error)")
    print()

    print("=" * 72)
    print(f"  COMPLETE -- l1 = {l1_derived:.4f} from epsilon_0, error = {error_pct:.4f}%")
    print("=" * 72)
    print()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
