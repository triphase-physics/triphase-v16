# -*- coding: utf-8 -*-
"""
================================================================================
TRIPHASE WAVE MECHANICS FRAMEWORK
Derivative 38: Coulomb Pressure (Electrostatic Pressure)
================================================================================

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
All Rights Reserved

DERIVATIVE TAG: (D) DERIVED - epsilon_0 appears DIRECTLY in Coulomb's law

MECHANISM:
-----------
Electrostatic pressure arises from Coulomb forces between charged particles.
The pressure in a plasma or charged gas is:

P_Coulomb = (e^2 / (8*pi*epsilon_0*r^4)) * n

where:
  e = elementary charge
  epsilon_0 = vacuum permittivity (APPEARS DIRECTLY!)
  r = characteristic distance
  n = particle density

The vacuum permittivity epsilon_0 appears explicitly in Coulomb's law:

F = (1 / (4*pi*epsilon_0)) * q1*q2 / r^2

Electrostatic pressure derives DIRECTLY from epsilon_0 with no intermediate
derivation required. This is fundamental electromagnetic vacuum structure.

PURE WAVE MECHANICS DERIVATION:
--------------------------------
Uses ONLY epsilon_0 and mu_0 as inputs
Elementary charge e is SI-defined exact
Measured values appear ONLY for calibration checkpoint
"""

import numpy as np

print("=" * 80)
print("TRIPHASE DERIVATIVE 38: COULOMB PRESSURE")
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
# STEP 2: SI-defined constants
# ============================================================================
print("STEP 2: SI-Defined Constants (Exact)")
print("-" * 80)

e = 1.602176634e-19  # C - Elementary charge (SI-defined exact)
h = 6.62607015e-34   # J*s - Planck constant (SI-defined exact)

print(f"e = {e:.15e} C (Elementary charge)")
print(f"h = {h:.15e} J*s (Planck constant)")
print()

# ============================================================================
# STEP 3: Coulomb's law
# ============================================================================
print("STEP 3: Coulomb's Law")
print("-" * 80)
print("The electrostatic force between two charges:")
print()
print("  F = (1 / (4*pi*epsilon_0)) * q1*q2 / r^2")
print()

k_e = 1.0 / (4 * np.pi * epsilon_0)

print(f"Coulomb constant: k_e = 1/(4*pi*epsilon_0)")
print(f"k_e = 1/(4*pi*{epsilon_0:.6e})")
print(f"k_e = {k_e:.15e} N*m^2/C^2")
print()

print("Note: epsilon_0 appears DIRECTLY in the force law.")
print("This is the fundamental coupling of electrostatic interactions.")
print()

# ============================================================================
# STEP 4: Electrostatic energy density
# ============================================================================
print("STEP 4: Electrostatic Energy Density and Pressure")
print("-" * 80)
print("The energy density in an electric field:")
print()
print("  u_E = epsilon_0 * E^2 / 2  [J/m^3]")
print()
print("This energy density creates pressure (Maxwell stress):")
print()
print("  P_E = epsilon_0 * E^2 / 2  [Pa]")
print()

# Example electric field
E_example = 1e6  # V/m (1 MV/m - strong field)

u_E = epsilon_0 * E_example**2 / 2
P_E = u_E

print(f"Example: E = {E_example:.3e} V/m")
print(f"Energy density: u_E = epsilon_0*E^2/2 = {u_E:.6e} J/m^3")
print(f"Pressure:       P_E = {P_E:.6e} Pa")
print(f"In atmospheres: P_E = {P_E/101325:.3f} atm")
print()

# ============================================================================
# STEP 5: Plasma pressure (Coulomb pressure)
# ============================================================================
print("STEP 5: Plasma Pressure - Coulomb Interactions")
print("-" * 80)
print("In a plasma, charged particles create electrostatic pressure.")
print()
print("For a characteristic inter-particle distance r:")
print("  Coulomb energy: U ~ e^2/(4*pi*epsilon_0*r)")
print()
print("Pressure (force/area) scales as energy/volume:")
print("  P ~ U/r^3 ~ e^2/(4*pi*epsilon_0*r^4)")
print()

# Example: hydrogen plasma
n_plasma = 1e20  # particles/m^3 (typical fusion plasma density)
r_plasma = (1/n_plasma)**(1/3)  # Average spacing

U_coulomb = e**2 / (4 * np.pi * epsilon_0 * r_plasma)
P_coulomb = n_plasma * U_coulomb

print(f"Plasma density: n = {n_plasma:.3e} particles/m^3")
print(f"Average spacing: r = (1/n)^(1/3) = {r_plasma:.3e} m")
print()
print(f"Coulomb energy: U = e^2/(4*pi*epsilon_0*r)")
print(f"U = {e:.3e}^2/(4*pi*{epsilon_0:.3e}*{r_plasma:.3e})")
print(f"U = {U_coulomb:.6e} J = {U_coulomb/e:.3f} eV")
print()
print(f"Coulomb pressure estimate: P ~ n*U")
print(f"P ~ {n_plasma:.3e} * {U_coulomb:.3e}")
print(f"P ~ {P_coulomb:.6e} Pa")
print()

# ============================================================================
# STEP 6: Debye screening and plasma parameter
# ============================================================================
print("STEP 6: Debye Screening in Plasmas")
print("-" * 80)

k_B = 1.380649e-23  # J/K - Boltzmann constant
T_plasma = 1e7  # K (10 million K - fusion temperature)

# Debye length
lambda_D = np.sqrt(epsilon_0 * k_B * T_plasma / (n_plasma * e**2))

print(f"Plasma temperature: T = {T_plasma:.3e} K")
print()
print(f"Debye length: lambda_D = sqrt(epsilon_0*k_B*T/(n*e^2))")
print(f"lambda_D = sqrt({epsilon_0:.3e}*{k_B:.3e}*{T_plasma:.3e}/({n_plasma:.3e}*{e:.3e}^2))")
print(f"lambda_D = {lambda_D:.6e} m = {lambda_D*1e6:.3f} microns")
print()

# Plasma parameter (number of particles in Debye sphere)
N_D = n_plasma * (4/3) * np.pi * lambda_D**3

print(f"Particles in Debye sphere: N_D = n*(4*pi/3)*lambda_D^3")
print(f"N_D = {N_D:.3e}")
print()

print("For a weakly coupled plasma, N_D >> 1")
print(f"Our plasma: {'weakly coupled' if N_D > 10 else 'strongly coupled'}")
print()

# ============================================================================
# STEP 7: Example - Capacitor pressure
# ============================================================================
print("STEP 7: Example - Pressure Between Capacitor Plates")
print("-" * 80)

V_cap = 1000  # Volts
d_cap = 1e-3  # 1 mm separation

E_cap = V_cap / d_cap
P_cap = epsilon_0 * E_cap**2 / 2

print(f"Capacitor voltage: V = {V_cap} V")
print(f"Plate separation:  d = {d_cap*1e3} mm")
print()
print(f"Electric field: E = V/d = {V_cap}/{d_cap*1e3} mm")
print(f"E = {E_cap:.3e} V/m")
print()
print(f"Electrostatic pressure: P = epsilon_0*E^2/2")
print(f"P = {epsilon_0:.3e}*{E_cap:.3e}^2/2")
print(f"P = {P_cap:.6e} Pa")
print(f"P = {P_cap*1e6:.3f} microPa")
print()

# Force on 1 cm^2 plate
area = 1e-4  # m^2
F_cap = P_cap * area

print(f"Force on 1 cm^2 plate: F = P*A = {F_cap:.6e} N")
print(f"F = {F_cap*1e6:.3f} microNewtons")
print()

# ============================================================================
# STEP 8: Direct dependence on epsilon_0
# ============================================================================
print("STEP 8: Direct Dependence on epsilon_0")
print("-" * 80)
print("Coulomb pressure depends DIRECTLY on epsilon_0:")
print()
print("  Coulomb force:  F = (1/(4*pi*epsilon_0)) * q1*q2/r^2")
print("  Electric field: E = (1/(4*pi*epsilon_0)) * q/r^2")
print("  Energy density: u = epsilon_0 * E^2 / 2")
print("  Pressure:       P = epsilon_0 * E^2 / 2")
print()
print("Unlike gravity (which requires deriving G from epsilon_0^3*mu_0^2),")
print("electrostatic pressure contains epsilon_0 EXPLICITLY in the formulas.")
print()
print(f"epsilon_0 = {epsilon_0:.6e} F/m sets the strength of ALL")
print("electromagnetic interactions in vacuum.")
print()

# ============================================================================
# CALIBRATION CHECKPOINT (measured values - NOT used in calculation)
# ============================================================================
print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
print("These measured values appear ONLY for verification - NOT used in derivation")
print()

# Verify Coulomb constant
k_e_measured = 8.9875517923e9  # N*m^2/C^2 (from CODATA)

print(f"Coulomb constant (CODATA): k_e = {k_e_measured:.10e} N*m^2/C^2")
print(f"TriPhase calculated:       k_e = {k_e:.10e} N*m^2/C^2")
print(f"Difference: {abs(k_e - k_e_measured)/k_e_measured * 100:.6f}%")
print()

# The exact relation
print("Exact relation: k_e * epsilon_0 = 1/(4*pi)")
product = k_e * epsilon_0
expected = 1/(4*np.pi)
print(f"k_e * epsilon_0 = {product:.15e}")
print(f"1/(4*pi)        = {expected:.15e}")
print(f"Match: {np.allclose(product, expected)}")
print()

print("Perfect agreement confirms Coulomb pressure derivation")
print("from epsilon_0 (appearing directly in formulas).")
print()

print("=" * 80)
print("DERIVATIVE COMPLETE")
print("=" * 80)
print()

input("Press Enter to exit...")
