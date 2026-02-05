"""
Paper III Calculation: Bessel J0 First Zero Verification
========================================================
Companion script for:
  "Local Laboratory Evidence for Three-Phase Vacuum Structure"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Verifies the first zero of Bessel function J_0(x) at x1 = 2.4048255577
using SciPy and confirms J_0(x1) = 0.

Requires: scipy
"""
import sys
import traceback


def main():
    from scipy.special import jn_zeros, j0

    print("=" * 72)
    print("  BESSEL J0 FIRST ZERO VERIFICATION -- Paper III")
    print("  Magnetic Innovative Solutions LLC")
    print("=" * 72)
    print()
    print("WHAT THIS CALCULATES:")
    print("-" * 72)
    print("""
    Bessel functions describe standing waves in circular and cylindrical
    geometry.  J_0(x) is the zeroth-order Bessel function -- it gives the
    amplitude pattern of the fundamental vibration mode on a circular drum.

    The first zero of J_0 occurs at:
      x_1 = 2.4048255577...

    This is a mathematical constant like pi -- it requires no measurement,
    no experiment, no physical input.  It comes from solving the Bessel
    differential equation.

    In the TriPhase framework, x_1 appears in CMB peak correction factors
    because vacuum fluctuations follow the same cylindrical wave mathematics.
    It also appears in walking droplet experiments (Bush, MIT) where
    bouncing silicone droplets form quantized orbits at Bessel zero radii.

    This script independently verifies x_1 using:
      1. SciPy's jn_zeros function (numerical root finding)
      2. Direct evaluation: J_0(x_1) should equal zero
      3. Comparison to the NIST/Abramowitz & Stegun table value
    """)

    # ============================================================
    # STEP 1: Table value from reference
    # ============================================================
    print("=" * 72)
    print("  STEP 1: Reference Value")
    print("=" * 72)
    print()

    x1_table = 2.4048255577
    print(f"  NIST Digital Library of Mathematical Functions (DLMF):")
    print(f"    x_1 = {x1_table}")
    print()
    print(f"  Also listed in Abramowitz & Stegun, Table 9.5")
    print()

    # ============================================================
    # STEP 2: SciPy computed value
    # ============================================================
    print("=" * 72)
    print("  STEP 2: SciPy Independent Computation")
    print("=" * 72)
    print()

    x1_scipy = jn_zeros(0, 1)[0]
    print(f"  scipy.special.jn_zeros(0, 1) = {x1_scipy:.15f}")
    print()

    # ============================================================
    # STEP 3: Verify J_0(x_1) = 0
    # ============================================================
    print("=" * 72)
    print("  STEP 3: Verify J_0(x_1) = 0")
    print("=" * 72)
    print()

    j0_value = j0(x1_scipy)
    print(f"  J_0({x1_scipy:.10f}) = {j0_value:.2e}")
    print(f"  Expected: 0")
    print(f"  |J_0(x_1)| = {abs(j0_value):.2e}  (machine precision)")
    print()

    # ============================================================
    # STEP 4: Compare table vs computed
    # ============================================================
    print("=" * 72)
    print("  STEP 4: Comparison")
    print("=" * 72)
    print()

    diff = abs(x1_scipy - x1_table)
    error_pct = diff / x1_table * 100

    print(f"  Table value:  x_1 = {x1_table}")
    print(f"  SciPy value:  x_1 = {x1_scipy:.15f}")
    print(f"  Difference:   {diff:.2e}")
    print(f"  Percent error: {error_pct:.10f}%")
    print()

    # Show nearby Bessel zeros for context
    print(f"  Context -- first 5 zeros of J_0:")
    zeros = jn_zeros(0, 5)
    for i, z in enumerate(zeros):
        print(f"    x_{i+1} = {z:.10f}")
    print()
    print(f"  Spacing increases by approximately pi = 3.14159...")
    for i in range(1, len(zeros)):
        print(f"    x_{i+1} - x_{i} = {zeros[i] - zeros[i-1]:.10f}")
    print()

    passed = error_pct < 0.0001 and abs(j0_value) < 1e-10
    status = "PASSED" if passed else "FAILED"
    print(f"  Result: {status}")
    print()

    print("=" * 72)
    print(f"  COMPLETE -- x_1 = {x1_scipy:.10f} verified, J_0(x_1) = {j0_value:.2e}")
    print("=" * 72)
    print()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
