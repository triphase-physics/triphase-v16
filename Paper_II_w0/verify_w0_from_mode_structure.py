"""
Paper II Calculation: Dark Energy w0 = -5/6 from Mode Structure
================================================================
Companion script for:
  "A Geometric Derivation of the Dark Energy Equation of State"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Derives the dark energy equation of state w0 = -5/6 = -0.8333
from three-phase vacuum mode geometry.  No physical constants needed.

Requires: (no dependencies -- pure arithmetic)
"""
import sys
import traceback


def main():
    print("=" * 72)
    print("  DARK ENERGY w0 FROM MODE STRUCTURE -- Paper II")
    print("  Magnetic Innovative Solutions LLC")
    print("=" * 72)
    print()
    print("WHAT THIS CALCULATES:")
    print("-" * 72)
    print("""
    Dark energy drives the accelerating expansion of the universe.
    Its equation of state parameter w0 relates pressure to energy density:
      P = w0 * rho * c^2

    Standard LCDM assumes w0 = -1 exactly (cosmological constant).
    DESI DR2 (2025) observes:  w0 = -0.838 +/- 0.055

    This derivation gets w0 from pure geometry -- no physical constants,
    no fitted parameters, no adjustable numbers.  Only mode counting:

      1. Three phases give 6 modes (3 phases x 2 quadratures)
      2. Vacuum energy splits: 1/6 kinetic (fluctuations), 5/6 background
      3. Background sector drives dark energy: w0 = -(5/6)

    The result -5/6 = -0.8333... matches DESI to 0.6% (0.09 sigma).

    NOTE: An alternate derivation via self-coupling suppression gives
    -181/216 = -0.83796, which is within 0.6% of -5/6. Both values fall
    within the DESI measurement uncertainty.
    """)

    # ============================================================
    # STEP 1: Mode structure
    # ============================================================
    print("=" * 72)
    print("  STEP 1: Mode Structure")
    print("=" * 72)
    print()

    phases = 3
    quadratures = 2
    total_modes = phases * quadratures

    print(f"  Three-phase electromagnetic structure:")
    print(f"    Phases:      {phases}   (120-degree separated oscillation modes)")
    print(f"    Quadratures: {quadratures}   (sine and cosine per phase)")
    print(f"    Total modes: {phases} x {quadratures} = {total_modes}")
    print()

    # ============================================================
    # STEP 2: Energy partition
    # ============================================================
    print("=" * 72)
    print("  STEP 2: Vacuum Energy Partition")
    print("=" * 72)
    print()

    kinetic_frac = 1 / total_modes
    background_frac = 1 - kinetic_frac

    print(f"  Of {total_modes} total modes:")
    print(f"    Kinetic sector (fluctuations):   1/{total_modes} = {kinetic_frac:.6f}")
    print(f"    Background sector (dark energy): {total_modes-1}/{total_modes} = {background_frac:.6f}")
    print()
    print(f"  The background sector consists of {total_modes-1} of {total_modes} modes.")
    print(f"  This gives the dark energy fraction: {total_modes-1}/{total_modes} = 5/6")
    print()

    # ============================================================
    # STEP 3: Equation of state
    # ============================================================
    print("=" * 72)
    print("  STEP 3: Dark Energy Equation of State")
    print("=" * 72)
    print()

    w0 = -background_frac
    w0_exact = -5 / 6

    print(f"  The background pressure opposes expansion (negative pressure).")
    print(f"  Equation of state: w0 = -(background fraction)")
    print()
    print(f"    w0 = -{total_modes-1}/{total_modes}")
    print(f"       = -5/6")
    print(f"       = {w0:.10f}")
    print()
    print(f"  Exact fraction check:")
    print(f"    -5/6     = {w0_exact:.10f}")
    print(f"    Computed = {w0:.10f}")
    print(f"    Match:     {abs(w0 - w0_exact):.2e}")
    print()

    # ============================================================
    # STEP 4: Compare to DESI observation
    # ============================================================
    print("=" * 72)
    print("  STEP 4: Comparison to DESI DR2 (2025)")
    print("=" * 72)
    print()

    desi_w0 = -0.838
    desi_err = 0.055
    deviation = abs(w0 - desi_w0)
    sigma = deviation / desi_err
    error_pct = deviation / abs(desi_w0) * 100

    print(f"  DESI DR2 + CMB + Pantheon+ (2025):")
    print(f"    w0 = {desi_w0} +/- {desi_err}")
    print()
    print(f"  This derivation:")
    print(f"    w0 = -5/6 = {w0:.4f}")
    print()
    print(f"  Deviation:  |{w0:.4f} - ({desi_w0})| = {deviation:.4f}")
    print(f"  Error:      {error_pct:.2f}%")
    print(f"  In sigma:   {deviation:.4f} / {desi_err} = {sigma:.2f} sigma")
    print()
    print(f"  LCDM prediction:  w0 = -1")
    print(f"  LCDM deviation:   |{desi_w0} - (-1)| = {abs(desi_w0 - (-1)):.3f}")
    print(f"  LCDM in sigma:    {abs(desi_w0 - (-1)):.3f} / {desi_err} = {abs(desi_w0 - (-1))/desi_err:.1f} sigma")
    print()

    # ============================================================
    # STEP 5: Alternate derivation via self-coupling suppression
    # ============================================================
    print("=" * 72)
    print("  STEP 5: Alternate Derivation (Self-Coupling Suppression)")
    print("=" * 72)
    print()

    print(f"  An alternate approach accounts for self-coupling suppression:")
    print()

    suppression = kinetic_frac ** 2
    print(f"  Suppression = (1/{total_modes})^2 = 1/{total_modes**2} = {suppression:.10f}")

    k = 1 - suppression
    print(f"  Effective scaling k = 1 - 1/{total_modes**2} = {k:.10f} = 35/36")

    n = k / 2
    print(f"  Geometric mean n = k/2 = {n:.10f} = 35/72")

    w0_alt = n / 3 - 1
    w0_alt_frac = -181 / 216
    print(f"  w0 = n/3 - 1 = {w0_alt:.10f} = -181/216 = {w0_alt_frac:.10f}")
    print()

    diff_alt = abs(w0_alt - w0)
    diff_alt_pct = diff_alt / abs(w0) * 100
    print(f"  Difference from -5/6:")
    print(f"    |-181/216 - (-5/6)| = {diff_alt:.6f}")
    print(f"    Error: {diff_alt_pct:.2f}%")
    print()
    print(f"  Both -5/6 and -181/216 are within DESI measurement uncertainty.")
    print()

    passed = abs(w0 - w0_exact) < 1e-10 and sigma < 1.0
    status = "PASSED" if passed else "FAILED"
    print(f"  Result: {status}")
    print(f"  (Fraction exact to machine precision, within 1 sigma of DESI)")
    print()

    print("=" * 72)
    print(f"  COMPLETE -- w0 = -5/6 = {w0:.4f}, DESI match at {sigma:.2f} sigma")
    print("=" * 72)
    print()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
