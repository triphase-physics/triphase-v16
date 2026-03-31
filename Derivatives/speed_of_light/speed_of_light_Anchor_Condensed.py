"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================

Derivative:  Speed of Light (c = 299792458 m/s)
Framework:   Anchor_Condensed
Version:     16.0
Generated:   2026-03-26
Status:      Active Development

Tag: (D) DERIVED - Condensed anchor chain (epsilon_0, mu_0 with derived symbols)

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
"""

import math

# === ANCHOR INPUTS ===
epsilon_0 = 8.8541878128e-12   # F/m
mu_0      = 1.25663706212e-6   # H/m
e         = 1.602176634e-19    # C (exact SI)

# === DERIVATION ===
print("=" * 80)
print("ANCHOR CONDENSED DERIVATION: Speed of Light")
print("Framework: Anchor_Condensed")
print("Tag: (D) DERIVED")
print("=" * 80)
print()

# Show anchor inputs
print("ANCHOR INPUTS:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m (electric permittivity of free space)")
print(f"  mu_0      = {mu_0:.14e} H/m (magnetic permeability of free space)")
print(f"  e         = {e:.15e} C (elementary charge, exact SI)")
print()

# DERIVATION: c = 1 / sqrt(epsilon_0 * mu_0)
print("DERIVATION:")
print("  c = 1 / sqrt(epsilon_0 * mu_0)")
print()
print("  This is the Maxwell relation — the speed of light emerges directly")
print("  from the electromagnetic properties of the vacuum.")
print()

# Component calculation
epsilon_mu_product = epsilon_0 * mu_0
sqrt_epsilon_mu = math.sqrt(epsilon_mu_product)
c_derived = 1.0 / sqrt_epsilon_mu

print("  Step-by-step:")
print(f"    epsilon_0 * mu_0 = {epsilon_0:.13e} * {mu_0:.14e}")
print(f"                     = {epsilon_mu_product:.6e} s^2/m^2")
print()
print(f"    sqrt(epsilon_0 * mu_0) = {sqrt_epsilon_mu:.6e} s/m")
print()
print(f"    c = 1 / {sqrt_epsilon_mu:.6e}")
print(f"      = {c_derived:.0f} m/s")
print()

# === CALIBRATION CHECKPOINT ===
codata_c = 299792458  # m/s (exact by definition)
error_pct = (c_derived - codata_c) / codata_c * 100.0

print("=" * 80)
print("CALIBRATION CHECKPOINT:")
print(f"  Derived:     {c_derived:.0f} m/s")
print(f"  CODATA:      {codata_c} m/s (exact by SI definition)")
print(f"  Error:       {error_pct:+.9f}%")
print("=" * 80)
print()

print("NOTE: The speed of light is exact in SI units (299792458 m/s by definition).")
print("      The small error comes from CODATA's epsilon_0 and mu_0 values,")
print("      which are measured quantities with finite precision.")
print()
print("      In TriPhase, c is DERIVED from epsilon_0 and mu_0,")
print("      demonstrating that light speed is an emergent property")
print("      of the electromagnetic structure of spacetime itself.")
print()

input("Press Enter to exit...")
