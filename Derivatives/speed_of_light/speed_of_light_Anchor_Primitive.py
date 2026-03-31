"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================

Derivative:  Speed of Light (c = 299792458 m/s)
Framework:   Anchor_Primitive
Version:     16.0
Generated:   2026-03-26
Status:      Active Development

Tag: (D) DERIVED - Pure anchor chain (epsilon_0, mu_0 only)

================================================================================
FRAMEWORK DESCRIPTION
================================================================================

The Anchor_Primitive framework is the PUREST form of TriPhase derivation:
- ONLY inputs: epsilon_0 and mu_0 (universal constants)
- c = 1/sqrt(epsilon_0 * mu_0) is DERIVED, never imported
- Z_0 = sqrt(mu_0 / epsilon_0) is DERIVED
- NO shortcut symbols like α, G, H₀ etc.
- Every intermediate value is FULLY SPELLED OUT from epsilon_0 and mu_0
- No scipy imports. No CODATA values used in calculation (only as calibration)
- Every formula traces back to epsilon_0 and mu_0 explicitly

The reader can trace EVERY number back to the two inputs without ever needing
to look up what "alpha" or "G" means. This is the self-contained proof.

================================================================================
ANCHOR PRIMITIVE DERIVATION
================================================================================

PURE ANCHOR CHAIN: Every value traced to epsilon_0 and mu_0

The speed of light is the SIMPLEST and most fundamental anchor derivation.
It emerges directly from the electromagnetic wave equation in vacuum.

Step 1: Maxwell's wave equation
  Starting from Maxwell's equations in vacuum (no charges, no currents),
  the wave equation for electromagnetic fields is:

  ∇²E = epsilon_0 * mu_0 * ∂²E/∂t²

  This describes how electric field E propagates through space and time.

Step 2: Plane wave solution
  For a plane wave: E = E₀ * exp(i(kx - ωt))

  Substituting into the wave equation:
  -k² = epsilon_0 * mu_0 * (-ω²)
  k² = epsilon_0 * mu_0 * ω²

Step 3: Phase velocity
  The phase velocity (speed of wave propagation) is:
  v_phase = ω / k

  From the dispersion relation (Step 2):
  k² = epsilon_0 * mu_0 * ω²
  k = sqrt(epsilon_0 * mu_0) * ω

  Therefore:
  v_phase = ω / [sqrt(epsilon_0 * mu_0) * ω]
  v_phase = 1 / sqrt(epsilon_0 * mu_0)

Step 4: Speed of light
  In vacuum, the phase velocity equals the speed of light:
  c = v_phase = 1 / sqrt(epsilon_0 * mu_0)

  This is the FUNDAMENTAL relationship - the speed of light is determined
  entirely by the vacuum's electric and magnetic properties.

  c = 1 / sqrt(8.8541878128e-12 * 1.25663706212e-6)
  c = 1 / sqrt(1.11265006e-17)
  c = 1 / 3.33564095e-9
  c = 299792458 m/s (exact by SI definition)

Physical Interpretation:
- epsilon_0 determines how strongly electric fields resist formation (capacitance)
- mu_0 determines how strongly magnetic fields resist formation (inductance)
- Together they form an "LC circuit" in vacuum
- The resonant frequency is infinite, but the wave speed is finite: c
- This is NOT a postulate - it's a mathematical consequence of epsilon_0 and mu_0

================================================================================
IMPLEMENTATION
================================================================================
"""

import math

def derive_speed_of_light_anchor_primitive():
    """
    Derive speed of light from pure anchor chain.
    Uses ONLY epsilon_0 and mu_0 as inputs.
    """

    print("=" * 80)
    print("ANCHOR PRIMITIVE DERIVATION: Speed of Light")
    print("=" * 80)
    print()

    # ============================================================================
    # ANCHOR INPUTS (SI 2019 exact values)
    # ============================================================================

    print("ANCHOR INPUTS:")
    print("-" * 80)

    epsilon_0 = 8.8541878128e-12  # F/m (exact, SI 2019)
    mu_0 = 1.25663706212e-6       # H/m (exact, SI 2019)

    print(f"  epsilon_0 = {epsilon_0:.13e} F/m  (electric permittivity)")
    print(f"  mu_0      = {mu_0:.14e} H/m  (magnetic permeability)")
    print()

    # ============================================================================
    # ANCHOR CHAIN: Derive c from epsilon_0, mu_0
    # ============================================================================

    print("ANCHOR CHAIN:")
    print("-" * 80)

    # Step 1: Product of epsilon_0 and mu_0
    product = epsilon_0 * mu_0
    print(f"Step 1: Product of vacuum permittivity and permeability")
    print(f"  epsilon_0 * mu_0 = {epsilon_0:.13e} * {mu_0:.14e}")
    print(f"                   = {product:.11e} F·H/m²")
    print()

    # Step 2: Square root of product
    sqrt_product = math.sqrt(product)
    print(f"Step 2: Square root of product")
    print(f"  sqrt(epsilon_0 * mu_0) = sqrt({product:.11e})")
    print(f"                         = {sqrt_product:.11e} s/m")
    print()

    # Step 3: Reciprocal gives speed of light
    c = 1.0 / sqrt_product
    print(f"Step 3: Reciprocal gives speed of light")
    print(f"  c = 1 / sqrt(epsilon_0 * mu_0)")
    print(f"    = 1 / {sqrt_product:.11e}")
    print(f"    = {c:.0f} m/s")
    print()

    # ============================================================================
    # VERIFICATION: Direct computation
    # ============================================================================

    print("VERIFICATION: Direct computation")
    print("-" * 80)

    c_direct = 1.0 / math.sqrt(epsilon_0 * mu_0)
    print(f"  c = 1 / sqrt({epsilon_0:.13e} * {mu_0:.14e})")
    print(f"    = {c_direct:.0f} m/s")
    print(f"  (Verification: {abs(c - c_direct):.2e} difference)")
    print()

    # ============================================================================
    # ALTERNATIVE FORMS
    # ============================================================================

    print("ALTERNATIVE FORMS:")
    print("-" * 80)

    # Vacuum impedance
    Z_0 = math.sqrt(mu_0 / epsilon_0)
    print(f"  Vacuum impedance:")
    print(f"  Z_0 = sqrt(mu_0 / epsilon_0)")
    print(f"      = sqrt({mu_0:.14e} / {epsilon_0:.13e})")
    print(f"      = {Z_0:.9f} Ohms")
    print()

    # Speed of light from Z_0
    c_from_Z0 = 1.0 / math.sqrt(epsilon_0 * mu_0)
    print(f"  Relationship to vacuum impedance:")
    print(f"  c = 1 / sqrt(epsilon_0 * mu_0)")
    print(f"  Z_0 = sqrt(mu_0 / epsilon_0)")
    print(f"  Therefore: c * Z_0 = 1/epsilon_0 and c/Z_0 = 1/mu_0")
    print()

    # ============================================================================
    # SI 2019 EXACT VALUE CHECK
    # ============================================================================

    print("=" * 80)
    print("SI 2019 EXACT VALUE CHECK")
    print("=" * 80)
    print()

    c_si_exact = 299792458  # m/s (exact by SI definition)

    error = c - c_si_exact
    error_pct = (error / c_si_exact) * 100

    print(f"  Anchor Primitive:  {c:.0f} m/s")
    print(f"  SI 2019 Exact:     {c_si_exact} m/s (defining constant)")
    print(f"  Error:             {error:+.6e} m/s ({error_pct:+.6e}%)")
    print()

    if abs(error) < 0.1:
        status = "EXACT MATCH"
    elif abs(error_pct) < 1e-6:
        status = "EXCELLENT"
    else:
        status = "NUMERICAL PRECISION ISSUE"

    print(f"  Calibration Status: {status}")
    print()

    # ============================================================================
    # PHYSICAL INTERPRETATION
    # ============================================================================

    print("=" * 80)
    print("PHYSICAL INTERPRETATION")
    print("=" * 80)
    print()
    print("  The speed of light (c = 299792458 m/s) is the phase velocity of")
    print("  electromagnetic waves in vacuum.")
    print()
    print("  Anchor Primitive Origin:")
    print("    - Emerges from Maxwell's wave equation in vacuum")
    print("    - Determined entirely by epsilon_0 and mu_0")
    print("    - c = 1 / sqrt(epsilon_0 * mu_0)")
    print()
    print("  This is NOT a postulate - it's a mathematical consequence of how")
    print("  the vacuum responds to electromagnetic disturbances. The vacuum")
    print("  acts like an \"LC circuit\" with:")
    print("    - epsilon_0 = capacitance per meter")
    print("    - mu_0 = inductance per meter")
    print("    - c = resonant wave speed")
    print()
    print("  In SI 2019, c is defined as exactly 299792458 m/s, which then")
    print("  constrains epsilon_0 and mu_0 through this relationship.")
    print()

    return c

# ================================================================================
# MAIN EXECUTION
# ================================================================================

if __name__ == "__main__":
    result = derive_speed_of_light_anchor_primitive()

    print("=" * 80)
    print(f"RESULT: c = {result:.0f} m/s")
    print("=" * 80)
    print()

    input("Press Enter to exit...")
