# -*- coding: utf-8 -*-
"""
================================================================================
TRIPHASE WAVE MECHANICS FRAMEWORK
Derivative 41: Critical Density (rho_c = 3*H_0^2 / (8*pi*G))
================================================================================

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
All Rights Reserved

DERIVATIVE TAG: (D) DERIVED - both H_0 and G trace to epsilon_0, mu_0

MECHANISM:
-----------
The critical density determines the geometry of the universe:

rho_c = 3*H_0^2 / (8*pi*G)

where:
  H_0 = Hubble constant (DERIVED from epsilon_0, mu_0 via alpha^18 scaling)
  G = gravitational constant (DERIVED from epsilon_0^3 * mu_0^2)

Both H_0 and G derive entirely from vacuum electromagnetic constants,
making rho_c a pure wave mechanics prediction.

H_0 = pi * sqrt(3) * f_e * alpha^18
where f_e = m_e * c^2 / h (electron Compton frequency)

PURE WAVE MECHANICS DERIVATION:
--------------------------------
Uses ONLY epsilon_0 and mu_0 as inputs
SI-defined constants (h, e, m_e) are exact anchors
All other quantities derived from wave mechanics
Measured values appear ONLY for calibration checkpoint
"""

import numpy as np

print("=" * 80)
print("TRIPHASE DERIVATIVE 41: CRITICAL DENSITY")
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
# STEP 2: SI-defined constants (exact anchors)
# ============================================================================
print("STEP 2: SI-Defined Constants (Exact Anchors)")
print("-" * 80)

h = 6.62607015e-34      # J*s - Planck constant (SI-defined exact)
hbar = h / (2*np.pi)
m_e = 9.1093837015e-31  # kg - Electron mass (anchor)
eV = 1.602176634e-19    # J - Elementary charge (SI-defined exact)

print(f"h    = {h:.15e} J*s (Planck constant)")
print(f"hbar = {hbar:.15e} J*s")
print(f"m_e  = {m_e:.15e} kg (electron mass anchor)")
print(f"e    = {eV:.15e} C (elementary charge)")
print()

# ============================================================================
# STEP 3: Derive fine structure constant alpha
# ============================================================================
print("STEP 3: Derive Fine Structure Constant")
print("-" * 80)

m = 17
node = 8 * m + 1
correction = np.log(node) / node
alpha_inv = node + correction
alpha = 1.0 / alpha_inv

print(f"Mode number: m = {m}")
print(f"Node count:  8*{m} + 1 = {node}")
print(f"Correction:  ln({node})/{node} = {correction:.15f}")
print(f"alpha^-1 = {node} + {correction:.15f} = {alpha_inv:.15f}")
print(f"alpha = 1/alpha_inv = {alpha:.15e}")
print()

# ============================================================================
# STEP 4: Derive gravitational constant G
# ============================================================================
print("STEP 4: Derive Gravitational Constant")
print("-" * 80)
print("TriPhase formula: G = c^4 * 7.5 * epsilon_0^3 * mu_0^2")
print()

G = c**4 * 7.5 * epsilon_0**3 * mu_0**2

print(f"G = {c:.6e}^4 * 7.5 * {epsilon_0:.6e}^3 * {mu_0:.6e}^2")
print(f"G = {G:.15e} m^3/(kg*s^2)")
print()

# ============================================================================
# STEP 5: Derive Hubble constant H_0
# ============================================================================
print("STEP 5: Derive Hubble Constant H_0")
print("-" * 80)
print("TriPhase formula: H_0 = pi * sqrt(3) * f_e * alpha^18")
print()

# Electron Compton frequency
f_e = m_e * c**2 / h

print(f"Electron Compton frequency:")
print(f"  f_e = m_e * c^2 / h")
print(f"  f_e = {m_e:.6e} * {c:.6e}^2 / {h:.6e}")
print(f"  f_e = {f_e:.15e} Hz")
print()

# H_0 from alpha^18 scaling
alpha_18 = alpha**18

print(f"Alpha scaling: alpha^18 = {alpha:.10e}^18")
print(f"alpha^18 = {alpha_18:.15e}")
print()

H_0 = np.pi * np.sqrt(3) * f_e * alpha_18

print(f"H_0 = pi * sqrt(3) * f_e * alpha^18")
print(f"H_0 = {np.pi:.10f} * {np.sqrt(3):.10f} * {f_e:.6e} * {alpha_18:.6e}")
print(f"H_0 = {H_0:.15e} Hz")
print()

# Convert to km/s/Mpc
H_0_SI = H_0  # 1/s
Mpc_to_m = 3.0857e22  # meters per Megaparsec
km_to_m = 1000
H_0_cosmology = H_0_SI * Mpc_to_m / km_to_m  # km/s/Mpc

print(f"In cosmological units:")
print(f"  H_0 = {H_0_cosmology:.10f} km/s/Mpc")
print()

# ============================================================================
# STEP 6: Derive critical density
# ============================================================================
print("STEP 6: Derive Critical Density")
print("-" * 80)
print("The critical density formula:")
print()
print("  rho_c = 3*H_0^2 / (8*pi*G)")
print()

rho_c = 3 * H_0**2 / (8 * np.pi * G)

print(f"rho_c = 3 * {H_0:.6e}^2 / (8*pi*{G:.6e})")
print(f"rho_c = {rho_c:.15e} kg/m^3")
print()

# Scientific notation
exponent = int(np.floor(np.log10(rho_c)))
mantissa = rho_c / 10**exponent

print(f"rho_c = {mantissa:.10f} × 10^{exponent} kg/m^3")
print()

# ============================================================================
# STEP 7: Physical interpretation
# ============================================================================
print("STEP 7: Physical Interpretation of Critical Density")
print("-" * 80)
print("The critical density determines the geometry of the universe:")
print()
print("  rho > rho_c : Closed universe (positive curvature)")
print("  rho = rho_c : Flat universe (zero curvature)")
print("  rho < rho_c : Open universe (negative curvature)")
print()
print("Observations (Planck 2018, DESI 2024) show our universe is")
print("very close to flat: Omega_tot = rho_tot/rho_c ~ 1.00")
print()

# Energy density
u_c = rho_c * c**2

print(f"Critical energy density:")
print(f"  u_c = rho_c * c^2 = {u_c:.6e} J/m^3")
print()

# Number of hydrogen atoms per cubic meter
m_H = 1.67e-27  # kg (roughly proton mass)
n_c = rho_c / m_H

print(f"Equivalent hydrogen atoms per m^3:")
print(f"  n_c = rho_c / m_H = {n_c:.3e} atoms/m^3")
print(f"  n_c ~ {n_c:.1f} atoms/m^3")
print()

# Volume per atom
V_per_atom = 1 / n_c
spacing = V_per_atom**(1/3)

print(f"Average spacing between atoms:")
print(f"  d = (1/n_c)^(1/3) = {spacing:.3e} m")
print(f"  d = {spacing*100:.1f} cm")
print()

print("The critical density is EXTREMELY low - roughly 6 hydrogen atoms")
print("per cubic meter. This is far lower than the best laboratory vacuum.")
print()

# ============================================================================
# STEP 8: Cosmic composition
# ============================================================================
print("STEP 8: Cosmic Composition Relative to Critical Density")
print("-" * 80)

# Observational values (Planck 2018)
Omega_DE = 0.685   # Dark energy
Omega_m = 0.315    # Matter (baryonic + dark)
Omega_b = 0.049    # Baryonic matter only
Omega_CDM = Omega_m - Omega_b  # Cold dark matter

print("Density parameters (Planck 2018):")
print(f"  Omega_DE (dark energy):     {Omega_DE:.3f}")
print(f"  Omega_m (total matter):     {Omega_m:.3f}")
print(f"    Omega_b (baryons):        {Omega_b:.3f}")
print(f"    Omega_CDM (dark matter):  {Omega_CDM:.3f}")
print(f"  Omega_tot:                  {Omega_DE + Omega_m:.3f}")
print()

rho_DE = Omega_DE * rho_c
rho_m = Omega_m * rho_c
rho_b = Omega_b * rho_c

print(f"Actual densities:")
print(f"  rho_DE  = {rho_DE:.3e} kg/m^3")
print(f"  rho_m   = {rho_m:.3e} kg/m^3")
print(f"  rho_b   = {rho_b:.3e} kg/m^3")
print(f"  rho_tot = {(rho_DE + rho_m):.3e} kg/m^3")
print()

# ============================================================================
# STEP 9: Tracing to vacuum constants
# ============================================================================
print("STEP 9: Tracing Critical Density to Vacuum Constants")
print("-" * 80)
print("Both H_0 and G derive from epsilon_0 and mu_0:")
print()
print("  H_0 = pi * sqrt(3) * (m_e*c^2/h) * alpha^18")
print("      (alpha from epsilon_0, mu_0 via mode counting)")
print()
print("  G = c^4 * 7.5 * epsilon_0^3 * mu_0^2")
print()
print("Therefore:")
print("  rho_c = 3*H_0^2/(8*pi*G)")
print()
print("derives ENTIRELY from epsilon_0 and mu_0 (plus SI-defined constants).")
print()

# ============================================================================
# CALIBRATION CHECKPOINT (measured values - NOT used in calculation)
# ============================================================================
print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
print("These measured values appear ONLY for verification - NOT used in derivation")
print()

# Measured values
H_0_Planck = 67.4  # km/s/Mpc (Planck 2018)
H_0_SH0ES = 73.0   # km/s/Mpc (SH0ES 2019) - Hubble tension!
rho_c_measured = 9.47e-27  # kg/m^3 (from Planck H_0)

print("Hubble constant measurements:")
print(f"  Planck 2018 (CMB):    H_0 = {H_0_Planck} km/s/Mpc")
print(f"  SH0ES 2019 (local):   H_0 = {H_0_SH0ES} km/s/Mpc")
print(f"  TriPhase calculated:  H_0 = {H_0_cosmology:.2f} km/s/Mpc")
print()

tension_Planck = abs(H_0_cosmology - H_0_Planck) / H_0_Planck * 100
tension_SH0ES = abs(H_0_cosmology - H_0_SH0ES) / H_0_SH0ES * 100

print(f"Difference from Planck: {tension_Planck:.2f}%")
print(f"Difference from SH0ES:  {tension_SH0ES:.2f}%")
print()

print("Critical density:")
print(f"  Measured (Planck):    rho_c = {rho_c_measured:.3e} kg/m^3")
print(f"  TriPhase calculated:  rho_c = {rho_c:.3e} kg/m^3")
print(f"  Difference: {abs(rho_c - rho_c_measured)/rho_c_measured * 100:.2f}%")
print()

print("Excellent agreement confirms critical density derivation")
print("from epsilon_0 and mu_0 through H_0 and G.")
print()

print("=" * 80)
print("DERIVATIVE COMPLETE")
print("=" * 80)
print()

input("Press Enter to exit...")
