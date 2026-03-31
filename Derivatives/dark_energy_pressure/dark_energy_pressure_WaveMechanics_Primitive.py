# -*- coding: utf-8 -*-
"""
================================================================================
TRIPHASE WAVE MECHANICS FRAMEWORK
Derivative 40: Dark Energy Pressure
================================================================================

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
All Rights Reserved

DERIVATIVE TAG: (D) DERIVED - w_0 from mode counting, rho_DE observational

MECHANISM:
-----------
Dark energy pressure follows the equation of state:

P_DE = w_0 * rho_DE * c^2

where:
  w_0 = -(17/18)^2 = -0.8919753... (DERIVED from mode counting)
  rho_DE = dark energy density (OBSERVATIONAL)
  c = speed of light (DERIVED from epsilon_0, mu_0)

The equation of state parameter w_0 derives from fundamental mode structure:
  - 17 = fine structure mode (gives alpha via 8*17+1 = 137)
  - 18 = 2*3*3 = combined mode structure

Negative pressure (w_0 < -1/3) drives cosmic acceleration.

PURE WAVE MECHANICS DERIVATION:
--------------------------------
Uses ONLY epsilon_0 and mu_0 as inputs
w_0 derives from mode counting
rho_DE is observational (measured)
Measured values appear ONLY for calibration checkpoint
"""

import numpy as np

print("=" * 80)
print("TRIPHASE DERIVATIVE 40: DARK ENERGY PRESSURE")
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
# STEP 1: Derive speed of light
# ============================================================================
print("STEP 1: Derive Speed of Light")
print("-" * 80)

c = 1.0 / np.sqrt(epsilon_0 * mu_0)

print(f"c = 1/sqrt(epsilon_0 * mu_0) = {c:.10e} m/s")
print()

# ============================================================================
# STEP 2: Derive equation of state parameter w_0
# ============================================================================
print("STEP 2: Derive Dark Energy Equation of State w_0")
print("-" * 80)
print("From TriPhase mode counting:")
print()

m_alpha = 17    # Fine structure mode
m_combined = 18  # Combined mode: 2*3*3

print(f"Fine structure mode:  m_alpha = {m_alpha}")
print(f"  (Gives alpha via nodal structure 8*{m_alpha}+1 = {8*m_alpha+1})")
print()
print(f"Combined mode: m_combined = {m_combined} = 2×3×3")
print(f"  Factor 2: electromagnetic duality (E/B)")
print(f"  Factor 3×3: spatial mode structure")
print()

ratio = m_alpha / m_combined
w_0 = -(ratio**2)

print(f"w_0 = -(m_alpha / m_combined)^2")
print(f"w_0 = -({m_alpha}/{m_combined})^2")
print(f"w_0 = -({ratio:.15f})^2")
print(f"w_0 = {w_0:.15f}")
print()

# ============================================================================
# STEP 3: Dark energy equation of state
# ============================================================================
print("STEP 3: Dark Energy Equation of State")
print("-" * 80)
print("The equation relating pressure to energy density:")
print()
print("  P = w * rho * c^2")
print()
print("where:")
print("  w = equation of state parameter")
print("  rho = energy density [kg/m^3]")
print("  c = speed of light [m/s]")
print()

print("Different cosmic components:")
print("  w = 0      : Matter (dust)")
print("  w = 1/3    : Radiation")
print("  w = -1     : Cosmological constant (Lambda)")
print(f"  w = {w_0:.3f} : TriPhase dark energy")
print()

print(f"Our w_0 = {w_0:.3f} < -1/3 implies:")
print("  - Negative pressure (repulsive gravity)")
print("  - Drives cosmic acceleration")
print("  - Dynamical, not constant vacuum energy")
print()

# ============================================================================
# STEP 4: Dark energy density (observational)
# ============================================================================
print("STEP 4: Dark Energy Density (Observational)")
print("-" * 80)
print("From cosmological observations (Planck 2018, DESI 2024):")
print()

# Dark energy density from observations
Omega_DE = 0.685  # Dark energy density parameter
rho_c = 9.47e-27  # kg/m^3 - critical density
rho_DE = Omega_DE * rho_c

print(f"Dark energy fraction: Omega_DE = {Omega_DE}")
print(f"Critical density:     rho_c = {rho_c:.3e} kg/m^3")
print()
print(f"Dark energy density: rho_DE = Omega_DE * rho_c")
print(f"rho_DE = {Omega_DE} * {rho_c:.3e}")
print(f"rho_DE = {rho_DE:.3e} kg/m^3")
print()

# ============================================================================
# STEP 5: Calculate dark energy pressure
# ============================================================================
print("STEP 5: Calculate Dark Energy Pressure")
print("-" * 80)

P_DE = w_0 * rho_DE * c**2

print(f"P_DE = w_0 * rho_DE * c^2")
print(f"P_DE = {w_0:.6f} * {rho_DE:.3e} * {c:.3e}^2")
print(f"P_DE = {P_DE:.6e} Pa")
print()

# Scientific notation
exponent = int(np.floor(np.log10(abs(P_DE))))
mantissa = P_DE / 10**exponent

print(f"P_DE = {mantissa:.10f} × 10^{exponent} Pa")
print()

print("The NEGATIVE pressure drives cosmic acceleration.")
print()

# ============================================================================
# STEP 6: Dark energy dominance
# ============================================================================
print("STEP 6: Dark Energy Dominance in Universe")
print("-" * 80)

# Matter density
Omega_m = 0.315  # Matter density parameter
rho_m = Omega_m * rho_c
P_m = 0  # Matter has negligible pressure (w_m ~ 0)

print(f"Matter density parameter: Omega_m = {Omega_m}")
print(f"Matter density:           rho_m = {rho_m:.3e} kg/m^3")
print(f"Matter pressure:          P_m ~ 0 Pa (dust-like)")
print()

print("Dark energy:")
print(f"  rho_DE = {rho_DE:.3e} kg/m^3")
print(f"  P_DE = {P_DE:.3e} Pa")
print()

ratio_rho = rho_DE / rho_m
print(f"Density ratio: rho_DE / rho_m = {ratio_rho:.3f}")
print()
print("Dark energy dominates the energy budget of the universe (~68%).")
print("Its negative pressure drives accelerating expansion.")
print()

# ============================================================================
# STEP 7: Acceleration criterion
# ============================================================================
print("STEP 7: Cosmic Acceleration Criterion")
print("-" * 80)
print("The universe accelerates when total equation of state w_tot < -1/3")
print()

w_crit = -1.0/3.0

print(f"Critical value: w_crit = -1/3 = {w_crit:.6f}")
print(f"Dark energy:    w_0 = {w_0:.6f}")
print()

if w_0 < w_crit:
    print(f"✓ w_0 < w_crit: Universe accelerates")
    margin = abs(w_0 - w_crit)
    print(f"  Margin: {margin:.6f} (strongly accelerating)")
else:
    print(f"✗ w_0 > w_crit: Universe decelerates")

print()

# Effective w_tot for universe
w_tot = (Omega_m * 0 + Omega_DE * w_0) / (Omega_m + Omega_DE)

print(f"Effective total w_tot = (Omega_m*w_m + Omega_DE*w_0)/(Omega_m + Omega_DE)")
print(f"w_tot = ({Omega_m}*0 + {Omega_DE}*{w_0:.3f})/({Omega_m}+{Omega_DE})")
print(f"w_tot = {w_tot:.6f}")
print()

if w_tot < w_crit:
    print(f"✓ w_tot = {w_tot:.3f} < -1/3: Current universe accelerates")
else:
    print(f"✗ w_tot > -1/3: Current universe decelerates")

print()

# ============================================================================
# STEP 8: Energy density vs pressure magnitude
# ============================================================================
print("STEP 8: Energy Density vs Pressure Magnitude")
print("-" * 80)

u_DE = rho_DE * c**2  # Energy density in J/m^3

print("Dark energy energy density:")
print(f"  u_DE = rho_DE * c^2 = {u_DE:.6e} J/m^3")
print()

print("Dark energy pressure:")
print(f"  P_DE = w_0 * rho_DE * c^2 = {P_DE:.6e} Pa")
print()

print("Relation:")
print(f"  P_DE = w_0 * u_DE")
print(f"  P_DE/u_DE = w_0 = {w_0:.6f}")
print()

print(f"Pressure magnitude: |P_DE| = {abs(P_DE):.3e} Pa")
print(f"Energy density:      u_DE = {u_DE:.3e} J/m^3")
print(f"Ratio: |P_DE|/u_DE = {abs(P_DE)/u_DE:.3f}")
print()

# ============================================================================
# CALIBRATION CHECKPOINT (measured values - NOT used in calculation)
# ============================================================================
print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
print("These measured values appear ONLY for verification - NOT used in derivation")
print()

# Planck 2018 w_0 (assuming wCDM model)
w_0_Planck = -1.03
w_0_err = 0.03

print("Planck 2018 (wCDM model):")
print(f"  w_0 = {w_0_Planck:.2f} ± {w_0_err:.2f}")
print()

print("TriPhase prediction:")
print(f"  w_0 = {w_0:.3f}")
print()

deviation = abs(w_0 - w_0_Planck)
sigma = deviation / w_0_err

print(f"Comparison:")
print(f"  Difference: {deviation:.3f}")
print(f"  Sigma: {sigma:.2f}")
print()

if sigma < 2.0:
    print(f"  ✓ Within 2-sigma of Planck measurement")
elif sigma < 3.0:
    print(f"  Within 3-sigma of Planck measurement")
else:
    print(f"  Beyond 3-sigma - suggests different dark energy model")

print()

print("Note: TriPhase w_0 = -0.892 predicts quintessence-like dark energy,")
print("distinct from cosmological constant (w = -1). This is testable with")
print("future surveys (DESI, Euclid, Roman).")
print()

print("=" * 80)
print("DERIVATIVE COMPLETE")
print("=" * 80)
print()

input("Press Enter to exit...")
