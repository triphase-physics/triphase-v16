"""
================================================================================
TriPhase V16 - Coulomb Pressure Derivative
Framework: THERMODYNAMICS
Tag: (D) - Pure derivation
================================================================================

THERMODYNAMICS FRAMEWORK:
Coulomb pressure arises from electrostatic repulsion — an electromagnetic
equation of state. For a charged sphere or plasma:

P_C = e²/(8πε₀r⁴) or more generally from Maxwell stress tensor

In hot plasmas, Debye screening modifies this: λ_D = √(ε₀k_BT/(ne²))
When λ_D < r, thermal screening dominates. When λ_D > r, Coulomb dominates.

Thermodynamic interpretation: Coulomb pressure = electrostatic free energy
gradient. At high T, thermal pressure wins. At low T, Coulomb pressure wins.

Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

epsilon_0 = 8.8541878128e-12
mu_0 = 1.25663706212e-6
e = 1.602176634e-19
c = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0 = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha = 1.0 / alpha_inv
hbar = Z_0 * e**2 / (4.0 * math.pi * alpha)
h = 2.0 * math.pi * hbar
G = c**4 * 7.5 * epsilon_0**3 * mu_0**2
r_e = 2.8179403262e-15
m_e = hbar * alpha / (c * r_e)
f_e = m_e * c**2 / hbar
T_17 = 17 * 18 // 2
mp_me = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p = m_e * mp_me
H_0 = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r = c**4 / (8.0 * math.pi * G)

print("=" * 80)
print("COULOMB PRESSURE - THERMODYNAMICS FRAMEWORK")
print("=" * 80)
print()
print("THERMODYNAMIC INTERPRETATION:")
print("Coulomb pressure = electrostatic equation of state.")
print("P_C competes with thermal pressure P_th = nk_BT in plasmas.")
print()

k_B = 1.380649e-23

print("COULOMB PRESSURE:")
print(f"  From Maxwell stress tensor:")
print(f"  P_C = ε₀E²/2 (electric field pressure)")
print(f"  For point charge: E = e/(4πε₀r²)")
print(f"  → P_C ∝ e²/(ε₀r⁴)")
print()

print("DEBYE SCREENING IN PLASMAS:")
print(f"  Debye length λ_D = √(ε₀k_BT/(ne²))")
print(f"  At r < λ_D: Coulomb force unscreened (strong)")
print(f"  At r > λ_D: Exponential screening (weak)")
print()

# Example: Solar core plasma
T_solar = 1.5e7  # K
n_solar = 1e32  # electrons/m³
lambda_D = math.sqrt(epsilon_0 * k_B * T_solar / (n_solar * e**2))

print(f"EXAMPLE: SOLAR CORE PLASMA")
print(f"  Temperature T = {T_solar:.2e} K")
print(f"  Electron density n = {n_solar:.2e} m⁻³")
print(f"  Debye length λ_D = {lambda_D:.3e} m = {lambda_D*1e10:.2f} Å")
print(f"  Thermal energy k_BT = {k_B*T_solar/e:.2f} keV")
print()

# Coulomb vs thermal
r_test = 1e-10  # 1 Angstrom
P_coulomb = e**2 / (8.0 * math.pi * epsilon_0 * r_test**4)
P_thermal = n_solar * k_B * T_solar

print(f"PRESSURE COMPARISON at r = {r_test*1e10:.1f} Å:")
print(f"  P_Coulomb ~ {P_coulomb:.3e} Pa")
print(f"  P_thermal = {P_thermal:.3e} Pa")
print(f"  Ratio P_th/P_C = {P_thermal/P_coulomb:.3e}")
print(f"  → Thermal pressure dominates (screened plasma)")
print()

print("THERMODYNAMIC REGIMES:")
print(f"  Cold/dense: Coulomb pressure dominates (Wigner crystal)")
print(f"  Hot/dilute: Thermal pressure dominates (ideal plasma)")
print(f"  Transition: Γ = e²/(4πε₀a k_BT) ~ 1 where a = n^(-1/3)")
print()

print("=" * 80)
print("CALIBRATION: Measured in plasmas, white dwarfs, laser fusion.")
print("=" * 80)

input("Press Enter to exit...")
