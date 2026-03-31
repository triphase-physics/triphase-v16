"""
================================================================================
TriPhase V16 - Critical Density Derivative
Framework: THERMODYNAMICS
Tag: (D) - Pure derivation
================================================================================

THERMODYNAMICS FRAMEWORK:
Critical density ρ_c = 3H₀²/(8πG) is the thermodynamic critical point of
the universe. At ρ = ρ_c, the universe is "flat" — balanced between expansion
and collapse, like a phase transition at critical temperature.

ρ < ρ_c: Open universe (expands forever — supercritical phase)
ρ = ρ_c: Flat universe (critical point)
ρ > ρ_c: Closed universe (recollapses — subcritical phase)

Thermodynamic interpretation: ρ_c is the density at which gravitational
binding energy balances kinetic expansion energy.

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
print("CRITICAL DENSITY - THERMODYNAMICS FRAMEWORK")
print("=" * 80)
print()
print("THERMODYNAMIC INTERPRETATION:")
print("ρ_c = 3H₀²/(8πG) is the critical point — flat universe.")
print("Like critical temperature in a phase transition.")
print()

rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)

print("CRITICAL DENSITY:")
print(f"  ρ_c = 3H₀²/(8πG)")
print(f"      = 3 × {H_0:.6e}² / (8π × {G:.6e})")
print(f"      = {rho_c:.6e} kg/m³")
print()

# Convert to useful units
n_H_critical = rho_c / m_p  # Hydrogen atoms per m³

print(f"  ρ_c = {n_H_critical:.3e} hydrogen atoms/m³")
print(f"      ≈ {n_H_critical*1e-6:.1f} atoms/cm³")
print(f"  This is incredibly dilute — better vacuum than lab achieves!")
print()

print("THERMODYNAMIC PHASE DIAGRAM:")
print(f"  ρ > ρ_c: Closed universe (k = +1)")
print(f"           Gravity overcomes expansion → eventual collapse")
print(f"           Thermodynamic analog: subcritical (bound state)")
print()
print(f"  ρ = ρ_c: Flat universe (k = 0)")
print(f"           Perfect balance (critical point)")
print(f"           Thermodynamic analog: phase coexistence")
print()
print(f"  ρ < ρ_c: Open universe (k = -1)")
print(f"           Expansion overcomes gravity → infinite expansion")
print(f"           Thermodynamic analog: supercritical (unbound)")
print()

# Observed universe
Omega_total = 1.00  # Very close to flat (Planck 2018)
print(f"OBSERVED UNIVERSE:")
print(f"  Ω_total = ρ_total/ρ_c = {Omega_total:.2f} (within measurement error)")
print(f"  The universe is almost exactly at the critical point!")
print(f"  This is the \"flatness problem\" — why so finely tuned?")
print()

print("THERMODYNAMIC FREE ENERGY:")
print(f"  E_kinetic = (1/2)Mṙ² (expansion energy)")
print(f"  E_potential = -GMm/r (gravitational binding)")
print(f"  At ρ = ρ_c: E_kinetic + E_potential = 0 (zero total energy)")
print(f"  Like a ball rolling on a flat potential — marginally bound")
print()

print("=" * 80)
print("CALIBRATION: ρ_c derived from measured H₀.")
print("Planck 2018 gives Ω_total = 1.000 ± 0.001 (remarkably flat!)")
print("=" * 80)

input("Press Enter to exit...")
