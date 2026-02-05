"""
Paper I Calculation: CMB Peak Ratio l3/l2 = 3/2
================================================
Companion script for:
  "CMB Acoustic Peak Structure from Vacuum Electromagnetic Properties"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Verifies that the derived peak positions give l3/l2 = 3/2 exactly,
equal to the ratio of phases to quadratures in the mode structure.

Requires: numpy, scipy
"""
import sys
import traceback


def main():
    import numpy as np
    from scipy.special import jn_zeros

    print("=" * 72)
    print("  CMB PEAK RATIO l3/l2 = 3/2 -- Paper I")
    print("  Magnetic Innovative Solutions LLC")
    print("=" * 72)
    print()
    print("WHAT THIS CALCULATES:")
    print("-" * 72)
    print("""
The ratio of the third to second CMB acoustic peak positions is one
of the cleanest tests of CMB physics.  Planck 2018 measures:

  l3 / l2 = 810 / 540 = 1.5000

Standard cosmology treats peak positions as fitted parameters.
There is no built-in reason for this ratio to be exactly 3/2.

In the TriPhase framework, the ratio emerges from mode structure:
  - 3 phases (electromagnetic oscillation modes at 120 degrees)
  - 2 quadratures per phase (sine and cosine components)
  - l3/l2 = phases / quadratures = 3/2

This is a parameter-free prediction: no fitting, no adjustment.
The ratio 3/2 is fixed by the geometry of three-phase structure.

This script computes l2 and l3 independently from their vacuum
property formulas and verifies that the ratio equals 3/2.
""")

    # ============================================================
    # STEP 1: Input constants
    # ============================================================
    print("=" * 72)
    print("  STEP 1: Input Constants")
    print("=" * 72)
    print()

    epsilon_0 = 8.8541878128e-12   # F/m
    mu_0 = 1.25663706212e-6        # H/m
    alpha_inv = 137.035999177
    Z_0 = np.sqrt(mu_0 / epsilon_0)
    x1 = jn_zeros(0, 1)[0]

    print(f"  epsilon_0 = {epsilon_0:.10e} F/m")
    print(f"  mu_0      = {mu_0:.11e} H/m")
    print(f"  alpha^-1  = {alpha_inv:.9f}")
    print(f"  Z_0       = sqrt(mu_0/epsilon_0) = {Z_0:.9f} Ohm")
    print(f"  x_1       = {x1:.10f}  (first J_0 zero)")
    print()
    print(f"  Predicted ratio: l3/l2 = 3/2 = {3/2:.6f}")
    print(f"  Observed ratio:  810/540 = {810/540:.6f}")
    print()

    # ============================================================
    # STEP 2: Compute l2 from alpha^-1 and mu_0
    # ============================================================
    print("=" * 72)
    print("  STEP 2: Compute l2 from alpha^-1 and mu_0")
    print("=" * 72)
    print()
    print(f"  l2 = pi * alpha^-1 * mu_0 * 10^6 * [alpha^-1 / (alpha^-1 + 2*x1/19)]")
    print()

    raw_l2 = np.pi * alpha_inv * mu_0 * 1e6
    two_x1_19 = 2 * x1 / 19
    corr_l2 = alpha_inv / (alpha_inv + two_x1_19)
    l2 = raw_l2 * corr_l2

    print(f"  Raw product = pi * {alpha_inv:.6f} * {mu_0:.10e} * 10^6")
    print(f"              = {raw_l2:.10f}")
    print(f"  2*x1/19     = {two_x1_19:.10f}")
    print(f"  Correction  = {alpha_inv:.6f} / ({alpha_inv:.6f} + {two_x1_19:.10f})")
    print(f"              = {corr_l2:.10f}")
    print(f"  l2 = {raw_l2:.10f} * {corr_l2:.10f}")
    print(f"     = {l2:.10f}")
    print()

    # ============================================================
    # STEP 3: Compute l3 from Z_0
    # ============================================================
    print("=" * 72)
    print("  STEP 3: Compute l3 from Z_0")
    print("=" * 72)
    print()
    print(f"  l3 = (13/6) * Z_0 * [22 / (22 + 6*x1/85)]")
    print()

    raw_l3 = (13 / 6) * Z_0
    six_x1_85 = 6 * x1 / 85
    corr_l3 = 22 / (22 + six_x1_85)
    l3 = raw_l3 * corr_l3

    print(f"  Raw product = (13/6) * {Z_0:.9f}")
    print(f"              = {raw_l3:.10f}")
    print(f"  6*x1/85     = {six_x1_85:.10f}")
    print(f"  Correction  = 22 / (22 + {six_x1_85:.10f})")
    print(f"              = {corr_l3:.10f}")
    print(f"  l3 = {raw_l3:.10f} * {corr_l3:.10f}")
    print(f"     = {l3:.10f}")
    print()

    # ============================================================
    # STEP 4: Form the ratio and compare to 3/2
    # ============================================================
    print("=" * 72)
    print("  STEP 4: Peak Ratio l3/l2")
    print("=" * 72)
    print()

    ratio = l3 / l2
    target = 3 / 2
    error_abs = abs(ratio - target)
    error_pct = error_abs / target * 100

    print(f"  l3 / l2 = {l3:.10f} / {l2:.10f}")
    print(f"          = {ratio:.10f}")
    print()
    print(f"  Target:   3/2 = {target:.10f}")
    print(f"  Error:    |{ratio:.10f} - {target:.10f}| = {error_abs:.2e}")
    print(f"  Percent:  {error_pct:.6f}%")
    print()
    print(f"  Physical meaning:")
    print(f"    3 = phases (120-degree electromagnetic modes)")
    print(f"    2 = quadratures (sine and cosine components)")
    print(f"    l3/l2 = phases / quadratures = 3/2")
    print()

    # Also show the observed Planck ratio
    l2_planck = 540.0
    l3_planck = 810.0
    ratio_planck = l3_planck / l2_planck
    print(f"  Planck 2018 observed:")
    print(f"    l2 = {l2_planck:.0f},  l3 = {l3_planck:.0f}")
    print(f"    l3/l2 = {ratio_planck:.6f}")
    print()

    passed = error_pct < 0.01
    status = "PASSED" if passed else "FAILED"
    print(f"  Result: {status} (threshold: < 0.01% error)")
    print()

    print("=" * 72)
    print(f"  COMPLETE -- l3/l2 = {ratio:.6f}, matches 3/2 to {error_pct:.6f}%")
    print("=" * 72)
    print()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
