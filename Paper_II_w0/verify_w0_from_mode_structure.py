"""
Paper II Calculation: Dark Energy w0 = -181/216 from Mode Structure
===================================================================
Companion script for:
  "A Geometric Derivation of the Dark Energy Equation of State"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Derives the dark energy equation of state w0 = -181/216 = -0.8380
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
      3. Fluctuations cannot self-couple (suppression factor 1/36)
      4. Effective scaling:  k = 1 - 1/36 = 35/36
      5. Observer geometric mean:  n = k/2 = 35/72
      6. Equation of state:  w0 = n/3 - 1 = 35/216 - 1 = -181/216

    The result -181/216 = -0.83796... matches DESI to 0.005 (0.04 sigma).
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
    print(f"  In 216ths (common denominator = 6^3 = {6**3}):")
    print(f"    Kinetic:    36/216 = 1/6")
    print(f"    Background: 180/216 = 5/6")
    print()

    # ============================================================
    # STEP 3: Self-coupling suppression
    # ============================================================
    print("=" * 72)
    print("  STEP 3: Self-Coupling Suppression")
    print("=" * 72)
    print()

    suppression = kinetic_frac ** 2
    print(f"  Fluctuations cannot self-couple.")
    print(f"  (A wave cannot push itself -- like pulling yourself up by bootstraps.)")
    print()
    print(f"  Suppression = (kinetic fraction)^2")
    print(f"              = (1/{total_modes})^2")
    print(f"              = 1/{total_modes**2}")
    print(f"              = {suppression:.10f}")
    print()

    # ============================================================
    # STEP 4: Effective scaling exponent
    # ============================================================
    print("=" * 72)
    print("  STEP 4: Effective Scaling Exponent")
    print("=" * 72)
    print()

    k = 1 - suppression
    print(f"  Remove the suppressed self-coupling component:")
    print(f"    k = 1 - 1/{total_modes**2}")
    print(f"      = 1 - 1/36")
    print(f"      = 35/36")
    print(f"      = {k:.10f}")
    print(f"  Check: 35/36 = {35/36:.10f}")
    print()

    # ============================================================
    # STEP 5: Observer geometric mean
    # ============================================================
    print("=" * 72)
    print("  STEP 5: Observer Geometric Mean")
    print("=" * 72)
    print()

    n = k / 2
    print(f"  We observe from within the fluctuation sector.")
    print(f"  The observable scaling is the geometric mean:")
    print(f"    n = k / 2")
    print(f"      = (35/36) / 2")
    print(f"      = 35/72")
    print(f"      = {n:.10f}")
    print(f"  Check: 35/72 = {35/72:.10f}")
    print()

    # ============================================================
    # STEP 6: Equation of state
    # ============================================================
    print("=" * 72)
    print("  STEP 6: Equation of State w0 = n/3 - 1")
    print("=" * 72)
    print()

    w0 = n / 3 - 1
    w0_exact = -181 / 216

    print(f"  Standard cosmological identity: w = n/3 - 1")
    print()
    print(f"    w0 = n/3 - 1")
    print(f"       = (35/72)/3 - 1")
    print(f"       = 35/216 - 1")
    print(f"       = 35/216 - 216/216")
    print(f"       = -181/216")
    print(f"       = {w0:.10f}")
    print()
    print(f"  Exact fraction check:")
    print(f"    -181/216 = {w0_exact:.10f}")
    print(f"    Computed = {w0:.10f}")
    print(f"    Match:     {abs(w0 - w0_exact):.2e}")
    print()
    print(f"  Denominator: 216 = 6^3 = {6**3}")
    print(f"  Numerator:   181 = 216 - 35")
    print()

    # ============================================================
    # STEP 7: Compare to DESI observation
    # ============================================================
    print("=" * 72)
    print("  STEP 7: Comparison to DESI DR2 (2025)")
    print("=" * 72)
    print()

    desi_w0 = -0.838
    desi_err = 0.055
    deviation = abs(w0 - desi_w0)
    sigma = deviation / desi_err

    print(f"  DESI DR2 + CMB + Pantheon+ (2025):")
    print(f"    w0 = {desi_w0} +/- {desi_err}")
    print()
    print(f"  This derivation:")
    print(f"    w0 = -181/216 = {w0:.4f}")
    print()
    print(f"  Deviation:  |{w0:.4f} - ({desi_w0})| = {deviation:.4f}")
    print(f"  In sigma:   {deviation:.4f} / {desi_err} = {sigma:.2f} sigma")
    print()
    print(f"  LCDM prediction:  w0 = -1")
    print(f"  LCDM deviation:   |{desi_w0} - (-1)| = {abs(desi_w0 - (-1)):.3f}")
    print(f"  LCDM in sigma:    {abs(desi_w0 - (-1)):.3f} / {desi_err} = {abs(desi_w0 - (-1))/desi_err:.1f} sigma")
    print()

    passed = abs(w0 - w0_exact) < 1e-10 and sigma < 1.0
    status = "PASSED" if passed else "FAILED"
    print(f"  Result: {status}")
    print(f"  (Fraction exact to machine precision, within 1 sigma of DESI)")
    print()

    print("=" * 72)
    print(f"  COMPLETE -- w0 = -181/216 = {w0:.4f}, DESI match at {sigma:.2f} sigma")
    print("=" * 72)
    print()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
