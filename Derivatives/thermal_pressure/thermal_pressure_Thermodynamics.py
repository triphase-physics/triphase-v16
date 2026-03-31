"""
================================================================================
TriPhase V16 - Thermal Pressure Derivative
Framework: THERMODYNAMICS
Tag: (D) - Pure derivation
================================================================================

THERMODYNAMICS FRAMEWORK:
Thermal pressure P = nk_BT IS the fundamental thermodynamic quantity.
This is the ideal gas law, derived from the kinetic theory of gases and
statistical mechanics. It emerges from equipartition theorem:
  ⟨E_kinetic⟩ = (3/2)k_BT per particle

The pressure is:
  P = nk_BT = (N/V)k_BT

where n = number density, k_B = Boltzmann constant, T = absolute temperature.

Extensions: Virial expansion P = nk_BT(1 + B₂n + B₃n² + ...) where B_i are
virial coefficients accounting for inter-particle interactions.

Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# ============================================================================
# STANDARD ANCHOR CHAIN
# ============================================================================
epsilon_0 = 8.8541878128e-12   # F/m
mu_0      = 1.25663706212e-6   # H/m
e         = 1.602176634e-19    # C (exact SI)
c         = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0       = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha     = 1.0 / alpha_inv
hbar      = Z_0 * e**2 / (4.0 * math.pi * alpha)
h         = 2.0 * math.pi * hbar
G         = c**4 * 7.5 * epsilon_0**3 * mu_0**2
r_e       = 2.8179403262e-15
m_e       = hbar * alpha / (c * r_e)
f_e       = m_e * c**2 / hbar
T_17      = 17 * 18 // 2   # = 153
mp_me     = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p       = m_e * mp_me
H_0       = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r      = c**4 / (8.0 * math.pi * G)

# ============================================================================
# THERMODYNAMIC DERIVATION - THERMAL PRESSURE
# ============================================================================

print("=" * 80)
print("THERMAL PRESSURE - THERMODYNAMICS FRAMEWORK")
print("=" * 80)
print()
print("THERMODYNAMIC INTERPRETATION:")
print("Thermal pressure P = nk_BT is THE fundamental thermodynamic equation.")
print("It emerges from Maxwell-Boltzmann statistics and equipartition theorem.")
print()

k_B = 1.380649e-23  # J/K

print("IDEAL GAS LAW:")
print(f"  P = nk_BT = (N/V)k_BT")
print(f"  where n = number density (particles/m³)")
print(f"        k_B = {k_B:.6e} J/K (Boltzmann constant)")
print(f"        T = absolute temperature (K)")
print()

print("DERIVATION FROM EQUIPARTITION:")
print(f"  Each translational degree of freedom contributes (1/2)k_BT")
print(f"  For 3D motion: ⟨E_kinetic⟩ = (3/2)k_BT per particle")
print(f"  Pressure = momentum transfer per collision:")
print(f"    P = (1/3)n⟨mv²⟩ = (1/3)n × 2⟨E_kin⟩ = nk_BT")
print()

# Examples
print("=" * 80)
print("EXAMPLES")
print("=" * 80)
print()

# Air at STP
T_STP = 273.15  # K
P_STP = 101325.0  # Pa
n_STP = P_STP / (k_B * T_STP)

print(f"1. AIR AT STANDARD TEMPERATURE & PRESSURE:")
print(f"  T = {T_STP:.2f} K (0°C)")
print(f"  P = {P_STP:.0f} Pa (1 atm)")
print(f"  n = P/(k_BT) = {n_STP:.3e} particles/m³")
print(f"  (Avogadro's number / 22.4 L = {6.022e23/0.0224:.3e} ✓)")
print()

# Solar core
T_solar_core = 1.57e7  # K
P_solar_core_gas = 2.5e16  # Pa (thermal pressure)
n_solar_core = P_solar_core_gas / (k_B * T_solar_core)

print(f"2. SOLAR CORE:")
print(f"  T = {T_solar_core:.2e} K")
print(f"  P_thermal = {P_solar_core_gas:.2e} Pa")
print(f"  n = {n_solar_core:.3e} particles/m³")
print()

print("VIRIAL EXPANSION:")
print(f"  P = nk_BT(1 + B₂n + B₃n² + ...)")
print(f"  Virial coefficients account for interactions.")
print()

print("=" * 80)
print("CALIBRATION: This IS the defining thermodynamic law.")
print("=" * 80)

input("Press Enter to exit...")
