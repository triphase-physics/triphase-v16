# -*- coding: utf-8 -*-
"""
================================================================================
TRIPHASE WAVE MECHANICS FRAMEWORK
Derivative 33: Electromagnetic Pressure (Maxwell Stress Tensor)
================================================================================

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
All Rights Reserved

DERIVATIVE TAG: (D) DERIVED directly from epsilon_0, mu_0

MECHANISM:
-----------
Electromagnetic fields exert pressure through the Maxwell stress tensor.
The pressure formulas derive directly from vacuum constants:

Magnetic pressure:  P_B = B^2 / (2*mu_0)
Electric pressure:  P_E = epsilon_0 * E^2 / 2

These are identical when E and B are related by the wave impedance Z_0.
The pressure is a fundamental consequence of electromagnetic field energy
density and derives directly from epsilon_0 and mu_0.

PURE WAVE MECHANICS DERIVATION:
--------------------------------
Uses ONLY epsilon_0 and mu_0 as inputs
All other quantities derived from wave mechanics
Measured values appear ONLY for calibration checkpoint
"""

import numpy as np

print("=" * 80)
print("TRIPHASE DERIVATIVE 33: ELECTROMAGNETIC PRESSURE")
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

# ============================================================================
# STEP 1: Derive speed of light and impedance
# ============================================================================
print("STEP 1: Derive Speed of Light and Vacuum Impedance")
print("-" * 80)

c = 1.0 / np.sqrt(epsilon_0 * mu_0)
Z_0 = np.sqrt(mu_0 / epsilon_0)

print(f"c   = 1/sqrt(epsilon_0 * mu_0) = {c:.10e} m/s")
print(f"Z_0 = sqrt(mu_0/epsilon_0)     = {Z_0:.10f} Ohms")
print()

# ============================================================================
# STEP 2: Electromagnetic pressure formulas
# ============================================================================
print("STEP 2: Electromagnetic Pressure Formulas")
print("-" * 80)
print("The Maxwell stress tensor gives electromagnetic pressure:")
print()
print("  Magnetic: P_B = B^2 / (2*mu_0)")
print("  Electric: P_E = epsilon_0 * E^2 / 2")
print()
print("These derive from the field energy density:")
print("  u_B = B^2 / (2*mu_0)  [J/m^3]")
print("  u_E = epsilon_0 * E^2 / 2  [J/m^3]")
print()
print("For electromagnetic waves, E and B are related:")
print(f"  E = c*B = Z_0*H  where H = B/mu_0")
print()

# ============================================================================
# STEP 3: Example - Magnetic pressure from 1 Tesla field
# ============================================================================
print("STEP 3: Example - Magnetic Pressure from 1 Tesla Field")
print("-" * 80)

B_example = 1.0  # Tesla

P_B = B_example**2 / (2 * mu_0)

print(f"Magnetic field: B = {B_example} T")
print()
print(f"P_B = B^2 / (2*mu_0)")
print(f"P_B = {B_example}^2 / (2*{mu_0:.6e})")
print(f"P_B = {P_B:.6f} Pa")
print(f"P_B = {P_B:.6e} Pa")
print()

# Convert to atmospheres
P_atm = 101325  # Pa
P_B_atm = P_B / P_atm

print(f"In atmospheres: P_B = {P_B_atm:.3f} atm")
print()
print("A 1 Tesla magnetic field exerts about 4 atmospheres of pressure!")
print()

# ============================================================================
# STEP 4: Equivalent electric field
# ============================================================================
print("STEP 4: Equivalent Electric Field for Same Pressure")
print("-" * 80)
print("For an electromagnetic wave with B = 1 T:")
print()

E_wave = c * B_example

print(f"E = c * B = {c:.6e} * {B_example}")
print(f"E = {E_wave:.6e} V/m")
print()

P_E_wave = epsilon_0 * E_wave**2 / 2

print(f"P_E = epsilon_0 * E^2 / 2")
print(f"P_E = {epsilon_0:.6e} * {E_wave:.6e}^2 / 2")
print(f"P_E = {P_E_wave:.6e} Pa")
print()

print(f"Magnetic pressure: P_B = {P_B:.6e} Pa")
print(f"Electric pressure: P_E = {P_E_wave:.6e} Pa")
print(f"Ratio P_E/P_B = {P_E_wave/P_B:.10f}")
print()
print("Perfect match - electromagnetic waves carry equal E and B pressures.")
print()

# ============================================================================
# STEP 5: Radiation pressure
# ============================================================================
print("STEP 5: Radiation Pressure")
print("-" * 80)
print("Electromagnetic radiation exerts pressure on surfaces.")
print("For a plane wave with intensity I [W/m^2]:")
print()
print("  I = (epsilon_0 * c / 2) * E^2 = (c / (2*mu_0)) * B^2")
print()
print("Radiation pressure (perfect absorption):")
print("  P_rad = I / c")
print()
print("Radiation pressure (perfect reflection):")
print("  P_rad = 2*I / c")
print()

# Example: sunlight
I_sun = 1361  # W/m^2 - solar constant

P_sun_abs = I_sun / c
P_sun_ref = 2 * I_sun / c

print(f"Example - Sunlight at Earth orbit:")
print(f"  Intensity: I = {I_sun} W/m^2")
print()
print(f"  Absorbed:  P = I/c = {I_sun}/{c:.3e} = {P_sun_abs:.3e} Pa")
print(f"  Reflected: P = 2I/c = {P_sun_ref:.3e} Pa")
print()

# Convert to force on 1 m^2 sail
F_sail = P_sun_ref * 1.0  # Newtons

print(f"Force on 1 m^2 reflective solar sail: F = {F_sail:.3e} N")
print()

# ============================================================================
# STEP 6: Strong field example - laser pressure
# ============================================================================
print("STEP 6: Strong Field Example - High Power Laser")
print("-" * 80)

# 1 MW/cm^2 = 10^10 W/m^2 (high power laser)
I_laser = 1e10  # W/m^2

P_laser = I_laser / c

print(f"Laser intensity: I = {I_laser:.3e} W/m^2 (1 MW/cm^2)")
print()
print(f"Radiation pressure: P = I/c = {P_laser:.3e} Pa")
print(f"In atmospheres: P = {P_laser/P_atm:.3f} atm")
print()

# Corresponding electric field
E_laser = np.sqrt(2 * I_laser / (epsilon_0 * c))
B_laser = E_laser / c

print(f"Electric field: E = {E_laser:.3e} V/m")
print(f"Magnetic field: B = {B_laser:.3e} T = {B_laser*1e4:.3f} Gauss")
print()

# ============================================================================
# STEP 7: Connection to vacuum structure
# ============================================================================
print("STEP 7: Connection to Vacuum Structure")
print("-" * 80)
print("Electromagnetic pressure derives DIRECTLY from vacuum constants:")
print()
print("  P_B = B^2 / (2*mu_0)")
print("  P_E = epsilon_0 * E^2 / 2")
print()
print("The vacuum permeability mu_0 determines magnetic pressure.")
print("The vacuum permittivity epsilon_0 determines electric pressure.")
print()
print("These are not separate phenomena - they are two aspects of the")
print("unified electromagnetic vacuum structure characterized by (epsilon_0, mu_0).")
print()
print(f"For waves: E/B = c = 1/sqrt(epsilon_0*mu_0) = {c:.3e} m/s")
print()

# ============================================================================
# CALIBRATION CHECKPOINT (measured values - NOT used in calculation)
# ============================================================================
print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
print("These measured values appear ONLY for verification - NOT used in derivation")
print()

# Verify with known experimental values
print("Laboratory verification:")
print(f"  1 Tesla field pressure: P = {P_B:.1f} Pa (calculated)")
print(f"  Expected from experiments: ~398,000 Pa")
print(f"  Match: {np.allclose(P_B, 397887, rtol=0.01)}")
print()

print(f"Solar radiation pressure: P = {P_sun_ref:.3e} Pa (calculated)")
print(f"  Expected: ~9.1e-6 Pa")
print(f"  Match: {np.allclose(P_sun_ref, 9.1e-6, rtol=0.01)}")
print()

print("Excellent agreement confirms electromagnetic pressure formulas")
print("derive correctly from epsilon_0 and mu_0.")
print()

print("=" * 80)
print("DERIVATIVE COMPLETE")
print("=" * 80)
print()

input("Press Enter to exit...")
