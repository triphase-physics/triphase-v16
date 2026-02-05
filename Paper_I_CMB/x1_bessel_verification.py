"""
Paper III Calculation: Bessel First Zero from Vacuum Properties
===============================================================
Companion script for:
  "Local Laboratory Evidence for Three-Phase Vacuum Structure"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Derives the first Bessel zero x_1 = 2.4048 from vacuum electromagnetic
properties using the balance equation with mode structure.

Requires: numpy
"""
import sys
import traceback


def main():
    import numpy as np

    print("=" * 72)
    print("  BESSEL FIRST ZERO FROM VACUUM PROPERTIES -- Paper III")
    print("  Magnetic Innovative Solutions LLC")
    print("=" * 72)
    print()
    print("WHAT THIS CALCULATES:")
    print("-" * 72)
    print("""
    The first zero of Bessel function J_0(x) is a mathematical constant:
      x_1 = 2.4048255577...

    It describes where the fundamental circular standing wave first
    crosses zero amplitude.  It appears in drum vibrations, electromagnetic
    waveguides, diffraction patterns, and walking droplet orbits.

    This script derives x_1 from vacuum electromagnetic properties:

      x_1 = [epsilon_0 * mu_0 * 10^18 / (6+1)^2] * cbrt(pi * Z_0)

    Where:
      epsilon_0 * mu_0 = 1/c^2  (vacuum product)
      10^18            = scale factor (SI to mode units)
      (6+1)^2 = 49    = mode interaction space:
                         6 oscillation modes + 1 ground mode, squared
                         for 2D radial symmetry (Bessel = 2D)
      cbrt(pi * Z_0)  = cylindrical geometry factor:
                         pi from circular symmetry
                         Z_0 from vacuum impedance
                         cube root for per-dimension scaling

    If the vacuum mode structure is real, it should encode the Bessel
    zero -- and it does, to within a few percent.
    """)

    # ============================================================
    # STEP 1: Vacuum properties
    # ============================================================
    print("=" * 72)
    print("  STEP 1: Vacuum Properties (CODATA 2022)")
    print("=" * 72)
    print()

    epsilon_0 = 8.8541878128e-12  # F/m
    mu_0 = 1.25663706212e-6       # H/m
    Z_0 = np.sqrt(mu_0 / epsilon_0)
    c = 1 / np.sqrt(epsilon_0 * mu_0)

    print(f"  epsilon_0 = {epsilon_0:.10e} F/m")
    print(f"  mu_0      = {mu_0:.11e} H/m")
    print(f"  Z_0       = sqrt(mu_0/epsilon_0) = {Z_0:.9f} Ohm")
    print(f"  c         = 1/sqrt(epsilon_0*mu_0) = {c:.3f} m/s")
    print()

    # ============================================================
    # STEP 2: Mode structure
    # ============================================================
    print("=" * 72)
    print("  STEP 2: Mode Structure")
    print("=" * 72)
    print()

    total_modes = 6   # 3 phases x 2 quadratures
    ground_mode = 1
    mode_factor = (total_modes + ground_mode) ** 2

    print(f"  3 phases x 2 quadratures = {total_modes} oscillation modes")
    print(f"  + {ground_mode} ground mode = {total_modes + ground_mode} total")
    print(f"  Squared for 2D radial symmetry: ({total_modes}+{ground_mode})^2 = {mode_factor}")
    print()
    print(f"  Why squared (not cubed)?")
    print(f"  Bessel functions describe 2D circular symmetry.")
    print(f"  We work in the radial plane, not 3D volume.")
    print()

    # ============================================================
    # STEP 3: Vacuum product and scaling
    # ============================================================
    print("=" * 72)
    print("  STEP 3: Vacuum Product and Scale Factor")
    print("=" * 72)
    print()

    product = epsilon_0 * mu_0
    product_scaled = product * 1e18

    print(f"  epsilon_0 * mu_0 = {product:.6e}")
    print(f"  = 1/c^2 = 1/({c:.3f})^2 = {1/c**2:.6e}  [check]")
    print()
    print(f"  Scale factor 10^18 converts SI to mode units:")
    print(f"    epsilon_0 ~ 10^-12, mu_0 ~ 10^-6")
    print(f"    Product ~ 10^-18")
    print(f"    10^18 normalizes to dimensionless mode space")
    print()
    print(f"  epsilon_0 * mu_0 * 10^18 = {product_scaled:.10f}")
    print()

    # ============================================================
    # STEP 4: Cylindrical geometry factor
    # ============================================================
    print("=" * 72)
    print("  STEP 4: Cylindrical Geometry Factor")
    print("=" * 72)
    print()

    piZ0 = np.pi * Z_0
    cbrt_piZ0 = piZ0 ** (1/3)

    print(f"  pi * Z_0 = {np.pi:.10f} * {Z_0:.9f}")
    print(f"           = {piZ0:.10f}")
    print()
    print(f"  cbrt(pi * Z_0) = ({piZ0:.6f})^(1/3)")
    print(f"                  = {cbrt_piZ0:.10f}")
    print()
    print(f"  Why cube root?  Same as w0 derivation:")
    print(f"  Extracting per-dimension contribution from volume quantity.")
    print()

    # ============================================================
    # STEP 5: Assemble the result
    # ============================================================
    print("=" * 72)
    print("  STEP 5: Calculate x_1")
    print("=" * 72)
    print()
    print(f"  x_1 = [epsilon_0 * mu_0 * 10^18 / (6+1)^2] * cbrt(pi * Z_0)")
    print()

    coefficient = product_scaled / mode_factor
    x1_calc = coefficient * cbrt_piZ0

    print(f"  Coefficient = {product_scaled:.10f} / {mode_factor}")
    print(f"              = {coefficient:.10f}")
    print()
    print(f"  x_1 = {coefficient:.10f} * {cbrt_piZ0:.10f}")
    print(f"      = {x1_calc:.10f}")
    print()

    # ============================================================
    # STEP 6: Compare to exact value
    # ============================================================
    print("=" * 72)
    print("  STEP 6: Comparison to Exact Bessel Zero")
    print("=" * 72)
    print()

    x1_exact = 2.4048255577

    error_abs = abs(x1_calc - x1_exact)
    error_pct = error_abs / x1_exact * 100

    print(f"  Derived:  x_1 = {x1_calc:.6f}")
    print(f"  Exact:    x_1 = {x1_exact}")
    print()
    print(f"  Absolute error: {error_abs:.6f}")
    print(f"  Percent error:  {error_pct:.2f}%")
    print()

    if error_pct < 5.0:
        print(f"  The vacuum properties encode the Bessel zero")
        print(f"  to within {error_pct:.1f}% -- consistent with the")
        print(f"  mode structure interpretation.")
    else:
        print(f"  The error is {error_pct:.1f}% -- larger than expected.")
        print(f"  The approximation captures the correct order and scale.")
    print()

    print("=" * 72)
    print(f"  COMPLETE -- x_1 = {x1_calc:.4f} from vacuum properties (exact: {x1_exact})")
    print("=" * 72)
    print()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
