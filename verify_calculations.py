"""
Paper I Master Verification: All TriPhase CMB Peak Calculations
===============================================================
Companion script for:
  "CMB Acoustic Peak Structure from Vacuum Electromagnetic Properties"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Runs 8 independent verifications of the TriPhase framework:
  1. [111] crystallographic projection (90 -> 120 deg)
  2. Bessel J0 first zero (x1 = 2.4048)
  3. Bessel zero ratios vs mode integers
  4. CMB peak l1 = 220 from epsilon_0
  5. CMB peak l2 = 532.84 from alpha^-1 and mu_0
  6. CMB peak l3 = 816.90 from Z_0
  7. Peak ratio l3/l2 = 3/2
  8. Dark energy w0 = -181/216 from mode structure

Requires: numpy, scipy
"""
import sys
import traceback


def main():
    import numpy as np
    from scipy.special import jn_zeros, j0

    print("=" * 72)
    print("  MASTER VERIFICATION SUITE -- Paper I")
    print("  Magnetic Innovative Solutions LLC")
    print("=" * 72)
    print()
    print("WHAT THIS CALCULATES:")
    print("-" * 72)
    print("""
    This script runs every calculation in the TriPhase framework and
    verifies each one independently.  It serves as a one-stop check:
    run it, and either all 8 verifications pass or they do not.

    The framework claims that CMB acoustic peak positions can be derived
    from three measured vacuum electromagnetic properties:

      epsilon_0 = 8.854e-12 F/m   (vacuum permittivity)
      mu_0      = 1.257e-6  H/m   (vacuum permeability)
      Z_0       = 376.73 Ohm      (vacuum impedance = sqrt(mu_0/epsilon_0))

    plus one mathematical constant:

      x1 = 2.4048  (first zero of Bessel function J_0)

    and one measured coupling constant:

      alpha^-1 = 137.036  (fine structure constant inverse)

    Each verification is independent.  No result depends on a previous
    verification passing.  All intermediate calculations are printed
    so you can follow along with a calculator.
    """)

    # ============================================================
    # PHYSICAL CONSTANTS (CODATA 2022)
    # ============================================================
    epsilon_0 = 8.8541878128e-12   # F/m
    mu_0 = 1.25663706212e-6        # H/m
    c = 299792458                   # m/s
    Z_0 = np.sqrt(mu_0 / epsilon_0)
    alpha_inv = 137.035999177
    pi = np.pi

    # Bessel J0 zeros (mathematical constants)
    bessel_zeros = jn_zeros(0, 10)
    x1 = bessel_zeros[0]
    x2 = bessel_zeros[1]
    x3 = bessel_zeros[2]
    x5 = bessel_zeros[4]
    x8 = bessel_zeros[7]

    print("=" * 72)
    print("  INPUT CONSTANTS")
    print("=" * 72)
    print()
    print(f"  epsilon_0  = {epsilon_0:.10e} F/m")
    print(f"  mu_0       = {mu_0:.11e} H/m")
    print(f"  Z_0        = {Z_0:.6f} Ohm")
    print(f"  alpha^-1   = {alpha_inv:.9f}")
    print(f"  x1 (J0)    = {x1:.10f}")
    print(f"  c          = {c} m/s")
    print()

    results = []


    # ============================================================
    # VERIFICATION 1: [111] Crystallographic Projection
    # ============================================================
    def verify_111_projection():
        print("=" * 72)
        print("  VERIFICATION 1: [111] Projection (90 deg -> 120 deg)")
        print("=" * 72)
        print()
        print(f"  Three perpendicular axes (x, y, z) viewed along the cube")
        print(f"  body diagonal [111] project to 120-degree separation.")
        print(f"  This is the geometric origin of three-phase structure.")
        print()

        x_hat = np.array([1, 0, 0])
        y_hat = np.array([0, 1, 0])
        z_hat = np.array([0, 0, 1])
        n_hat = np.array([1, 1, 1]) / np.sqrt(3)

        def project(v, n):
            return v - np.dot(v, n) * n

        x_proj = project(x_hat, n_hat)
        y_proj = project(y_hat, n_hat)
        z_proj = project(z_hat, n_hat)

        def angle_between(v1, v2):
            cos_t = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            return np.degrees(np.arccos(np.clip(cos_t, -1, 1)))

        a_xy = angle_between(x_proj, y_proj)
        a_yz = angle_between(y_proj, z_proj)
        a_zx = angle_between(z_proj, x_proj)

        print(f"  n_hat = [1,1,1]/sqrt(3) = [{n_hat[0]:.10f}, ...]")
        print(f"  Angle(x_proj, y_proj) = {a_xy:.6f} deg")
        print(f"  Angle(y_proj, z_proj) = {a_yz:.6f} deg")
        print(f"  Angle(z_proj, x_proj) = {a_zx:.6f} deg")
        print(f"  Sum = {a_xy + a_yz + a_zx:.6f} deg")
        print()

        max_err = max(abs(a_xy - 120), abs(a_yz - 120), abs(a_zx - 120))
        passed = max_err < 0.0001
        print(f"  Max deviation from 120: {max_err:.2e} deg")
        print(f"  Result: {'PASSED' if passed else 'FAILED'}")
        print()
        return passed


    # ============================================================
    # VERIFICATION 2: Bessel J0 First Zero
    # ============================================================
    def verify_bessel_zero():
        print("=" * 72)
        print("  VERIFICATION 2: Bessel J0 First Zero (x1 = 2.4048)")
        print("=" * 72)
        print()
        print(f"  Bessel functions describe standing waves in circular geometry.")
        print(f"  x1 is the first radius where J_0 crosses zero.")
        print()

        x1_table = 2.4048255577
        x1_scipy = jn_zeros(0, 1)[0]
        j0_val = j0(x1_scipy)

        print(f"  Table value (NIST):  x1 = {x1_table}")
        print(f"  SciPy computed:      x1 = {x1_scipy:.15f}")
        print(f"  J_0(x1) = {j0_val:.2e}  (should be 0)")
        print()

        diff = abs(x1_scipy - x1_table)
        error = diff / x1_table * 100
        print(f"  Difference: {diff:.2e}")
        print(f"  Error: {error:.10f}%")

        passed = error < 0.0001 and abs(j0_val) < 1e-10
        print(f"  Result: {'PASSED' if passed else 'FAILED'}")
        print()
        return passed


    # ============================================================
    # VERIFICATION 3: Bessel Zero Ratios vs Mode Integers
    # ============================================================
    def verify_bessel_ratios():
        print("=" * 72)
        print("  VERIFICATION 3: Bessel Zero Ratios vs Mode Integers")
        print("=" * 72)
        print()
        print(f"  Ratios of J_0 zeros match mode-counting integer ratios")
        print(f"  to better than 0.5% accuracy.")
        print()

        ratios = [
            ("x3/x1", x3, x1, 18, 5, "phases*total / active"),
            ("x8/x2", x8, x2, 22, 5, "coupling / active"),
            ("x5/x1", x5, x1, 137, 22, "alpha^-1 / coupling"),
        ]

        print(f"  {'Ratio':<8s}  {'Bessel':>10s}  {'Integer':>10s}  {'Error':>8s}  Origin")
        print(f"  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*24}")

        all_ok = True
        max_error = 0
        for label, num_z, den_z, i_num, i_den, origin in ratios:
            b_ratio = num_z / den_z
            i_ratio = i_num / i_den
            err = abs(b_ratio - i_ratio) / i_ratio * 100
            max_error = max(max_error, err)
            if err >= 1.0:
                all_ok = False
            print(f"  {label:<8s}  {b_ratio:>10.6f}  {i_num:>4d}/{i_den:<4d}  {err:>7.4f}%  {origin}")

        print()
        print(f"  Max error: {max_error:.4f}%")
        print(f"  Result: {'PASSED' if all_ok else 'FAILED'}")
        print()
        return all_ok


    # ============================================================
    # VERIFICATION 4: CMB Peak l1 = 220 from epsilon_0
    # ============================================================
    def verify_l1():
        print("=" * 72)
        print("  VERIFICATION 4: l1 = 220 from epsilon_0")
        print("=" * 72)
        print()
        print(f"  Formula: l1 = 5^2 * epsilon_0 * 10^12 * [22/(22 + x1/18)]")
        print()

        raw = 25 * epsilon_0 * 1e12
        print(f"  Step 1: Raw = 25 * {epsilon_0:.10e} * 10^12 = {raw:.6f}")

        x1_18 = x1 / 18
        print(f"  Step 2: x1/18 = {x1:.10f} / 18 = {x1_18:.6f}")

        denom = 22 + x1_18
        print(f"  Step 3: 22 + x1/18 = {denom:.6f}")

        corr = 22 / denom
        print(f"  Step 4: Correction = 22 / {denom:.6f} = {corr:.6f}")

        l1 = raw * corr
        print(f"  Step 5: l1 = {raw:.6f} * {corr:.6f} = {l1:.4f}")

        obs = 220
        obs_err = 1
        error = abs(l1 - obs) / obs * 100
        sigma = abs(l1 - obs) / obs_err
        print()
        print(f"  Observed (Planck 2018): l1 = {obs} +/- {obs_err}")
        print(f"  Derived: l1 = {l1:.4f}")
        print(f"  Error: {error:.4f}% ({sigma:.2f} sigma)")

        passed = error < 1.0
        print(f"  Result: {'PASSED' if passed else 'FAILED'}")
        print()
        return passed, l1


    # ============================================================
    # VERIFICATION 5: CMB Peak l2 = 532.84 from alpha^-1 and mu_0
    # ============================================================
    def verify_l2():
        print("=" * 72)
        print("  VERIFICATION 5: l2 = 532.84 from alpha^-1 and mu_0")
        print("=" * 72)
        print()
        print(f"  Formula: l2 = pi * alpha^-1 * mu_0 * 10^6 * [alpha^-1/(alpha^-1 + 2*x1/19)] * [1 - x1/180]")
        print()

        raw = pi * alpha_inv * mu_0 * 1e6
        print(f"  Step 1: Raw = pi * {alpha_inv:.6f} * {mu_0:.10e} * 10^6 = {raw:.6f}")

        two_x1_19 = 2 * x1 / 19
        print(f"  Step 2: 2*x1/19 = 2 * {x1:.10f} / 19 = {two_x1_19:.6f}")

        denom = alpha_inv + two_x1_19
        print(f"  Step 3: alpha^-1 + 2*x1/19 = {denom:.6f}")

        corr = alpha_inv / denom
        print(f"  Step 4: Quantization = {alpha_inv:.6f} / {denom:.6f} = {corr:.6f}")

        l2_base = raw * corr
        print(f"  Step 5: l2_base = {raw:.6f} * {corr:.6f} = {l2_base:.4f}")

        rarefaction = 1 - x1 / 180
        print(f"  Step 6: Rarefaction = [1 - x1/180] = 1 - {x1/180:.6f} = {rarefaction:.6f}")

        l2 = l2_base * rarefaction
        print(f"  Step 7: l2 = {l2_base:.4f} * {rarefaction:.6f} = {l2:.4f}")

        obs = 533.2
        obs_err = 5.2
        error = abs(l2 - obs) / obs * 100
        sigma = abs(l2 - obs) / obs_err
        print()
        print(f"  Observed (Planck 2018, Gaussian fit): l2 = {obs} +/- {obs_err}")
        print(f"  Derived: l2 = {l2:.4f}")
        print(f"  Error: {error:.4f}% ({sigma:.2f} sigma)")

        passed = error < 1.0
        print(f"  Result: {'PASSED' if passed else 'FAILED'}")
        print()
        return passed, l2, l2_base


    # ============================================================
    # VERIFICATION 6: CMB Peak l3 = 816.90 from Z_0
    # ============================================================
    def verify_l3():
        print("=" * 72)
        print("  VERIFICATION 6: l3 = 816.90 from Z_0")
        print("=" * 72)
        print()
        print(f"  Formula: l3 = (13/6) * Z_0 * [22/(22 + 6*x1/85)] * [1 + x1/282]")
        print()

        raw = (13 / 6) * Z_0
        print(f"  Step 1: Raw = (13/6) * {Z_0:.6f} = {raw:.6f}")

        six_x1_85 = 6 * x1 / 85
        print(f"  Step 2: 6*x1/85 = 6 * {x1:.10f} / 85 = {six_x1_85:.6f}")

        denom = 22 + six_x1_85
        print(f"  Step 3: 22 + 6*x1/85 = {denom:.6f}")

        corr = 22 / denom
        print(f"  Step 4: Quantization = 22 / {denom:.6f} = {corr:.6f}")

        l3_base = raw * corr
        print(f"  Step 5: l3_base = {raw:.6f} * {corr:.6f} = {l3_base:.4f}")

        compression = 1 + x1 / 282
        print(f"  Step 6: Compression = [1 + x1/282] = 1 + {x1/282:.6f} = {compression:.6f}")

        l3 = l3_base * compression
        print(f"  Step 7: l3 = {l3_base:.4f} * {compression:.6f} = {l3:.4f}")

        obs = 816.9
        obs_err = 2.8
        error = abs(l3 - obs) / obs * 100
        sigma = abs(l3 - obs) / obs_err
        print()
        print(f"  Observed (Planck 2018, Gaussian fit): l3 = {obs} +/- {obs_err}")
        print(f"  Derived: l3 = {l3:.4f}")
        print(f"  Error: {error:.4f}% ({sigma:.2f} sigma)")

        passed = error < 1.0
        print(f"  Result: {'PASSED' if passed else 'FAILED'}")
        print()
        return passed, l3, l3_base


    # ============================================================
    # VERIFICATION 7: Peak Ratio l3/l2 = 3/2
    # ============================================================
    def verify_ratio(l2_base, l3_base):
        print("=" * 72)
        print("  VERIFICATION 7: Peak Ratio l3_base/l2_base = 3/2")
        print("=" * 72)
        print()
        print(f"  TriPhase predicts l3_base/l2_base = 3/2 = phases/quadratures.")
        print(f"  The ratio applies to base quantization positions (before mode")
        print(f"  coupling corrections shift peaks to observed positions).")
        print()

        ratio = l3_base / l2_base
        target = 3 / 2
        error = abs(ratio - target) / target * 100

        print(f"  l3_base = {l3_base:.4f}")
        print(f"  l2_base = {l2_base:.4f}")
        print(f"  l3_base/l2_base = {ratio:.6f}")
        print(f"  Target = 3/2 = {target:.6f}")
        print(f"  Error: {error:.6f}%")

        passed = error < 0.01
        print(f"  Result: {'PASSED' if passed else 'FAILED'}")
        print()
        return passed


    # ============================================================
    # VERIFICATION 8: w0 = -181/216 from Mode Structure
    # ============================================================
    def verify_w0():
        print("=" * 72)
        print("  VERIFICATION 8: w0 = -181/216 from Mode Structure")
        print("=" * 72)
        print()
        print(f"  Derives dark energy equation of state from vacuum geometry.")
        print(f"  No free parameters -- pure mode counting.")
        print()

        phases = 3
        quadratures = 2
        N = phases * quadratures
        print(f"  Step 1: N = {phases} phases x {quadratures} quadratures = {N} modes")

        kinetic = 1 / N
        print(f"  Step 2: Kinetic partition = 1/N = 1/{N} = {kinetic:.6f}")
        print(f"          In 216ths: {N**2}/{N**3} = 36/216")

        supp = kinetic ** 2
        print(f"  Step 3: Self-coupling suppression = (1/{N})^2 = 1/{N**2} = {supp:.10f}")

        k = 1 - supp
        print(f"  Step 4: Effective scaling k = 1 - 1/{N**2} = {N**2 - 1}/{N**2} = {k:.6f}")

        n = k / 2
        print(f"  Step 5: Geometric mean n = k/2 = {k:.6f}/2 = {n:.6f}")

        w0 = n / 3 - 1
        w0_frac = -181 / 216
        print(f"  Step 6: w0 = n/3 - 1 = {n:.6f}/3 - 1 = {w0:.10f}")
        print(f"          -181/216 = {w0_frac:.10f}")
        print(f"          Difference: {abs(w0 - w0_frac):.2e}")

        desi = -0.838
        desi_err = 0.055
        sigma = abs(w0 - desi) / desi_err
        print()
        print(f"  DESI DR2 (2025): w0 = {desi} +/- {desi_err}")
        print(f"  Derived: w0 = {w0:.6f}")
        print(f"  Deviation: {sigma:.2f} sigma")

        passed = abs(w0 - w0_frac) < 1e-10 and sigma < 1.0
        print(f"  Result: {'PASSED' if passed else 'FAILED'}")
        print()
        return passed


    # ============================================================
    # RUN ALL VERIFICATIONS
    # ============================================================
    print()
    print("=" * 72)
    print("  RUNNING 8 VERIFICATIONS")
    print("=" * 72)
    print()

    results.append(("1. [111] Projection (90->120 deg)", verify_111_projection()))
    results.append(("2. Bessel Zero x1 = 2.4048", verify_bessel_zero()))
    results.append(("3. Bessel Ratios vs Mode Integers", verify_bessel_ratios()))

    p4, l1_val = verify_l1()
    results.append(("4. CMB Peak l1 = 220", p4))

    p5, l2_val, l2_base = verify_l2()
    results.append(("5. CMB Peak l2 = 532.84", p5))

    p6, l3_val, l3_base = verify_l3()
    results.append(("6. CMB Peak l3 = 816.90", p6))

    results.append(("7. Peak Ratio l3_base/l2_base = 3/2", verify_ratio(l2_base, l3_base)))
    results.append(("8. w0 = -181/216", verify_w0()))

    # ============================================================
    # SUMMARY
    # ============================================================
    print("=" * 72)
    print("  VERIFICATION SUMMARY")
    print("=" * 72)
    print()

    all_passed = True
    for name, passed in results:
        symbol = "[OK]  " if passed else "[FAIL]"
        print(f"  {symbol} {name}: {'PASSED' if passed else 'FAILED'}")
        if not passed:
            all_passed = False

    print()
    n_pass = sum(1 for _, p in results if p)
    n_total = len(results)

    if all_passed:
        print(f"  ALL {n_total} VERIFICATIONS PASSED")
        print(f"  CMB peaks derived from vacuum electromagnetic constants.")
        print(f"  w0 = -181/216 derived from mode structure geometry.")
    else:
        print(f"  {n_pass}/{n_total} VERIFICATIONS PASSED -- REVIEW REQUIRED")

    print()
    print("=" * 72)
    print(f"  COMPLETE -- {n_pass}/{n_total} verifications passed")
    print("=" * 72)
    print()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
