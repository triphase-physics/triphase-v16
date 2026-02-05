"""
Paper II Calculation: w0 Balance Equation from Vacuum Properties
================================================================
Companion script for:
  "A Geometric Derivation of the Dark Energy Equation of State"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Derives w0 from the balance equation using vacuum permittivity and
permeability: w0 = -2 * cbrt(mu_0) / (5^3 * cbrt(epsilon_0)).

Requires: numpy
"""
import sys
import traceback


def main():
    import numpy as np

    print("=" * 72)
    print("  w0 BALANCE EQUATION VERIFICATION -- Paper II")
    print("  Magnetic Innovative Solutions LLC")
    print("=" * 72)
    print()
    print("WHAT THIS CALCULATES:")
    print("-" * 72)
    print("""
    Paper II derives the dark energy equation of state w0 = -181/216
    from pure mode-counting geometry (see verify_w0_from_mode_structure.py).

    This script provides an ALTERNATIVE derivation using the measured
    vacuum electromagnetic properties epsilon_0 and mu_0 directly.

    The balance equation:
      w0 = -2 * cbrt(mu_0) / (5^3 * cbrt(epsilon_0))

    Physical interpretation:
      - Left side:  w0 * 125 * cbrt(epsilon_0) = pressure EOS x capacitive
      - Right side: -2 * cbrt(mu_0) = inductive contribution

    This says: the dark energy equation of state is the ratio of
    magnetic (inductive) to electric (capacitive) vacuum energy,
    scaled by the coupling mode count 5^3 = 125.

    Cube roots appear because pressure = energy/volume and volume
    scales as length^3 -- the cube root extracts the per-dimension
    contribution.
    """)

    # ============================================================
    # STEP 1: Vacuum electromagnetic properties
    # ============================================================
    print("=" * 72)
    print("  STEP 1: Vacuum Electromagnetic Properties (CODATA 2022)")
    print("=" * 72)
    print()

    epsilon_0 = 8.8541878128e-12  # F/m
    mu_0 = 1.25663706212e-6       # H/m
    Z_0 = np.sqrt(mu_0 / epsilon_0)
    c = 1 / np.sqrt(epsilon_0 * mu_0)

    print(f"  epsilon_0 = {epsilon_0:.10e} F/m  (vacuum permittivity)")
    print(f"  mu_0      = {mu_0:.11e} H/m  (vacuum permeability)")
    print(f"  Z_0       = sqrt(mu_0/epsilon_0) = {Z_0:.9f} Ohm")
    print(f"  c         = 1/sqrt(epsilon_0*mu_0) = {c:.3f} m/s")
    print()

    # ============================================================
    # STEP 2: Cube roots (pressure = energy/volume)
    # ============================================================
    print("=" * 72)
    print("  STEP 2: Cube Roots of Vacuum Properties")
    print("=" * 72)
    print()

    cbrt_eps = epsilon_0 ** (1/3)
    cbrt_mu = mu_0 ** (1/3)

    print(f"  Pressure = Energy / Volume,  Volume = Length^3")
    print(f"  Cube root extracts the per-dimension contribution.")
    print()
    print(f"  cbrt(epsilon_0) = ({epsilon_0:.10e})^(1/3)")
    print(f"                   = {cbrt_eps:.10e}")
    print()
    print(f"  cbrt(mu_0)      = ({mu_0:.11e})^(1/3)")
    print(f"                   = {cbrt_mu:.10e}")
    print()

    # ============================================================
    # STEP 3: Mode structure
    # ============================================================
    print("=" * 72)
    print("  STEP 3: Mode Structure")
    print("=" * 72)
    print()

    phases = 3
    quadratures = 2
    total_modes = phases * quadratures
    coupling_modes = total_modes - 1

    print(f"  3 phases x 2 quadratures = {total_modes} total modes")
    print(f"  {total_modes} - 1 ground mode = {coupling_modes} coupling modes")
    print(f"  5^3 = {coupling_modes**3}  (coupling modes cubed)")
    print()

    # ============================================================
    # STEP 4: The balance equation
    # ============================================================
    print("=" * 72)
    print("  STEP 4: Balance Equation Calculation")
    print("=" * 72)
    print()
    print(f"  Formula: w0 = -2 * cbrt(mu_0) / (5^3 * cbrt(epsilon_0))")
    print()

    numerator = -2 * cbrt_mu
    denominator = (coupling_modes ** 3) * cbrt_eps
    w0_calc = numerator / denominator

    print(f"  Numerator:   -2 * cbrt(mu_0)")
    print(f"             = -2 * {cbrt_mu:.10e}")
    print(f"             = {numerator:.10e}")
    print()
    print(f"  Denominator: 5^3 * cbrt(epsilon_0)")
    print(f"             = {coupling_modes**3} * {cbrt_eps:.10e}")
    print(f"             = {denominator:.10e}")
    print()
    print(f"  w0 = {numerator:.10e} / {denominator:.10e}")
    print(f"     = {w0_calc:.10f}")
    print()

    # ============================================================
    # STEP 5: Compare to DESI and mode-counting result
    # ============================================================
    print("=" * 72)
    print("  STEP 5: Comparison")
    print("=" * 72)
    print()

    w0_mode = -181/216
    w0_desi = -0.838
    w0_desi_err = 0.055

    diff_mode = abs(w0_calc - w0_mode)
    diff_desi = abs(w0_calc - w0_desi)
    sigma_desi = diff_desi / w0_desi_err

    print(f"  Balance equation: w0 = {w0_calc:.6f}")
    print(f"  Mode counting:   w0 = -181/216 = {w0_mode:.6f}")
    print(f"  DESI DR2 (2025): w0 = {w0_desi} +/- {w0_desi_err}")
    print()
    print(f"  vs mode counting: |{w0_calc:.6f} - ({w0_mode:.6f})| = {diff_mode:.6f}")
    print(f"  vs DESI:          |{w0_calc:.6f} - ({w0_desi})| = {diff_desi:.6f} ({sigma_desi:.2f} sigma)")
    print()

    print(f"  Physical interpretation:")
    print(f"    Left:  w0 * 125 * cbrt(epsilon_0) = pressure x capacitive")
    print(f"    Right: -2 * cbrt(mu_0) = inductive contribution")
    print(f"    Balance: capacitive storage = inductive release")
    print()

    passed = sigma_desi < 1.0
    status = "PASSED" if passed else "FAILED"
    print(f"  Result: {status} (within 1 sigma of DESI)")
    print()

    print("=" * 72)
    print(f"  COMPLETE -- w0 = {w0_calc:.6f} from balance equation")
    print("=" * 72)
    print()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
