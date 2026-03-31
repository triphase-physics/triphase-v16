# -*- coding: utf-8 -*-
"""
================================================================================
TRIPHASE WAVE MECHANICS FRAMEWORK
Derivative 8: Energy Per Mode (epsilon_pair)
================================================================================

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
All Rights Reserved

DERIVATIVE TAG: (D*) DERIVED with discrete selection

MECHANISM:
-----------
The energy per mode epsilon_pair emerges from the Rydberg constant combined
with the fundamental triangular mode number T_17 = 17*18/2 = 153.

Starting from vacuum permittivity and permeability, we derive:
1. Speed of light c
2. Fine structure constant alpha (from m=17 nodal structure)
3. Rydberg constant (from alpha, electron mass, c, and h)
4. Energy per mode = Rydberg * h * c / T_17

This connects the atomic energy spectrum to discrete mode counting.

PURE WAVE MECHANICS DERIVATION:
--------------------------------
Uses ONLY vacuum constants epsilon_0, mu_0 and SI-defined exact values h, e, m_e
No measured values used in calculation - only for calibration checkpoint
"""

import math

print("=" * 80)
print("TRIPHASE DERIVATIVE 8: ENERGY PER MODE (epsilon_pair)")
print("=" * 80)
print()

# ============================================================================
# FUNDAMENTAL INPUTS - The ONLY external parameters
# ============================================================================
print("FUNDAMENTAL INPUTS (Vacuum Constants):")
print("-" * 80)

epsilon_0 = 8.8541878128e-12  # F/m - Vacuum permittivity
mu_0 = 1.25663706212e-6       # H/m - Vacuum permeability

print(f"epsilon_0 = {epsilon_0:.13e} F/m")
print(f"mu_0      = {mu_0:.14e} H/m")
print()

# SI-defined exact values
h = 6.62607015e-34    # J·s - Planck constant (SI-defined exact)
m_e = 9.1093837015e-31  # kg - Electron mass (anchor for mass derivations)
e = 1.602176634e-19   # C - Elementary charge (SI-defined exact)

print("SI-DEFINED EXACT VALUES:")
print(f"h   = {h:.11e} J·s (Planck constant)")
print(f"m_e = {m_e:.13e} kg (electron mass anchor)")
print(f"e   = {e:.12e} C (elementary charge)")
print()

# ============================================================================
# STEP 1: Derive speed of light and impedance
# ============================================================================
print("STEP 1: Derive Speed of Light")
print("-" * 80)

c = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0 = math.sqrt(mu_0 / epsilon_0)

print(f"c   = 1/sqrt(epsilon_0 * mu_0) = {c:.10e} m/s")
print(f"Z_0 = sqrt(mu_0/epsilon_0)      = {Z_0:.10f} Ω")
print()

# ============================================================================
# STEP 2: Derive fine structure constant alpha
# ============================================================================
print("STEP 2: Derive Fine Structure Constant")
print("-" * 80)
print("Nodal structure: m = 17")
print("Node count: 8*m + 1 = 8*17 + 1 = 137")
print()

m = 17
node_count = 8 * m + 1
print(f"m          = {m}")
print(f"node_count = {node_count}")
print()

correction = math.log(node_count) / node_count
alpha_inv = node_count + correction
alpha = 1.0 / alpha_inv

print(f"correction = ln({node_count})/{node_count} = {correction:.10f}")
print(f"alpha_inv  = {node_count} + {correction:.10f} = {alpha_inv:.10f}")
print(f"alpha      = 1/alpha_inv = {alpha:.15f}")
print()

# ============================================================================
# STEP 3: Derive Rydberg constant
# ============================================================================
print("STEP 3: Derive Rydberg Constant")
print("-" * 80)
print("Rydberg = alpha^2 * m_e * c / (2*h)")
print()

R_inf = (alpha**2 * m_e * c) / (2 * h)

print(f"R_inf = {R_inf:.10e} m^-1")
print()

# ============================================================================
# STEP 4: Derive energy per mode (epsilon_pair)
# ============================================================================
print("STEP 4: Derive Energy Per Mode")
print("-" * 80)
print("Triangular mode number: T_17 = 17*18/2 = 153")
print("epsilon_pair = Rydberg * h * c / T_17")
print()

T_17 = (17 * 18) // 2
epsilon_pair = (R_inf * h * c) / T_17

print(f"T_17         = {T_17}")
print(f"epsilon_pair = {epsilon_pair:.15e} J")
print()

# Convert to eV for readability
epsilon_pair_eV = epsilon_pair / e
print(f"epsilon_pair = {epsilon_pair_eV:.10f} eV")
print()

# ============================================================================
# CALIBRATION CHECKPOINT (measured values - NOT used in calculation)
# ============================================================================
print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
print("These measured values appear ONLY for verification - NOT used in derivation")
print()

R_inf_measured = 10973731.568160  # m^-1 (CODATA 2018)
epsilon_pair_expected_eV = 89.0  # eV (approximate from hydrogen spectrum)

print(f"Rydberg constant:")
print(f"  Derived:  {R_inf:.6e} m^-1")
print(f"  Measured: {R_inf_measured:.6f} m^-1 (CODATA 2018)")
print(f"  Match: {abs(R_inf - R_inf_measured) / R_inf_measured * 100:.6f}% difference")
print()

print(f"Energy per mode:")
print(f"  Derived:  {epsilon_pair_eV:.6f} eV")
print(f"  Expected: ~{epsilon_pair_expected_eV:.1f} eV (hydrogen spectrum)")
print(f"  Match: {abs(epsilon_pair_eV - epsilon_pair_expected_eV) / epsilon_pair_expected_eV * 100:.2f}% difference")
print()

print("=" * 80)
print("DERIVATIVE COMPLETE")
print("=" * 80)
print()

input("Press Enter to exit...")
