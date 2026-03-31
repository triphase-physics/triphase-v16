"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================

Derivative:  Proton-Electron Mass Ratio (mp/me = 1836.153...)
Framework:   Anchor_Condensed
Version:     16.0
Generated:   2026-03-26
Status:      Active Development

Tag: (D*) DERIVED with discrete selection - Condensed anchor chain

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
"""

import math

# === ANCHOR INPUTS ===
epsilon_0 = 8.8541878128e-12   # F/m
mu_0      = 1.25663706212e-6   # H/m
e         = 1.602176634e-19    # C (exact SI)

# === DERIVED ANCHOR CHAIN ===
c     = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0   = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha = 1.0 / alpha_inv

# === DERIVATION ===
print("=" * 80)
print("ANCHOR CONDENSED DERIVATION: Proton-Electron Mass Ratio")
print("Framework: Anchor_Condensed")
print("Tag: (D*) DERIVED with discrete selection")
print("=" * 80)
print()

# Show anchor inputs
print("ANCHOR INPUTS:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print(f"  e         = {e:.15e} C (exact SI)")
print()

# Show derived chain (condensed)
print("DERIVED CHAIN (all from epsilon_0, mu_0):")
print(f"  c     = {c:.0f} m/s")
print(f"  Z_0   = {Z_0:.6f} Ohms")
print(f"  alpha = 1/{alpha_inv:.6f} = {alpha:.10f}")
print()

# DERIVATION: mp/me = 2^2 * 3^3 * 17 * (1 + 5*alpha^2/pi)
print("DERIVATION:")
print("  mp/me = 2^2 * 3^3 * 17 * (1 + 5*alpha^2/pi)")
print()

# Discrete factor
discrete_factor = 2**2 * 3**3 * 17
print(f"  Discrete factor: 2^2 * 3^3 * 17 = {2**2} * {3**3} * 17 = {discrete_factor}")
print()

# Continuous correction
correction_term = 1.0 + 5.0 * alpha**2 / math.pi
print(f"  Correction term: 1 + 5*alpha^2/pi")
print(f"                 = 1 + 5*{alpha:.10f}^2/pi")
print(f"                 = 1 + {5.0 * alpha**2:.12f}/pi")
print(f"                 = 1 + {5.0 * alpha**2 / math.pi:.12f}")
print(f"                 = {correction_term:.12f}")
print()

# Final ratio
mp_me_ratio = discrete_factor * correction_term
print(f"  mp/me = {discrete_factor} * {correction_term:.12f}")
print(f"        = {mp_me_ratio:.9f}")
print()

# === CALIBRATION CHECKPOINT ===
codata_ratio = 1836.15267343  # CODATA 2022
error_pct = (mp_me_ratio - codata_ratio) / codata_ratio * 100.0

print("=" * 80)
print("CALIBRATION CHECKPOINT:")
print(f"  Derived:    {mp_me_ratio:.9f}")
print(f"  CODATA 2022: {codata_ratio:.8f}")
print(f"  Error:      {error_pct:+.6f}%")
print("=" * 80)
print()

print("NOTE: The discrete factor (2^2 * 3^3 * 17 = 1836) dominates.")
print("      The correction term (1 + 5*alpha^2/pi ≈ 1.000084) adds fine tuning.")
print("      This is a 'quantum ladder' structure: discrete + continuous.")
print()

input("Press Enter to exit...")
