"""
================================================================================
TriPhase V16 - Matter Density Derivative
Framework: THERMODYNAMICS
Tag: (D*H) - Derived with hypothetical components
================================================================================

THERMODYNAMICS FRAMEWORK:
Matter density ρ_m = Ω_m × ρ_c where Ω_m ≈ 0.315 is determined by the
cosmic thermal history. Two key thermodynamic events:

1. Matter-radiation equality (T_eq ≈ 9400 K): When ρ_m = ρ_rad
2. Big Bang Nucleosynthesis (T ≈ 0.8 MeV): Sets baryon fraction

Baryogenesis and dark matter freeze-out are thermodynamic phase transitions
that determine Ω_m. The observed value encodes the thermal history of the
early universe.

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
print("MATTER DENSITY - THERMODYNAMICS FRAMEWORK")
print("=" * 80)
print()
print("THERMODYNAMIC INTERPRETATION:")
print("ρ_m = Ω_m × ρ_c is set by thermal history of early universe.")
print("Baryogenesis and DM freeze-out are thermodynamic processes.")
print()

Omega_m = 0.315  # Matter density parameter (Planck 2018)
Omega_b = 0.049  # Baryon density parameter
Omega_DM = Omega_m - Omega_b  # Dark matter

rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)
rho_m = Omega_m * rho_c
rho_b = Omega_b * rho_c
rho_DM = Omega_DM * rho_c

print("MATTER DENSITY:")
print(f"  Critical density ρ_c        = {rho_c:.6e} kg/m³")
print(f"  Matter fraction Ω_m         = {Omega_m:.3f}")
print(f"    Baryons Ω_b               = {Omega_b:.3f}")
print(f"    Dark matter Ω_DM          = {Omega_DM:.3f}")
print()
print(f"  ρ_m = Ω_m × ρ_c")
print(f"      = {Omega_m:.3f} × {rho_c:.6e}")
print(f"      = {rho_m:.6e} kg/m³")
print()

k_B = 1.380649e-23

# Matter-radiation equality
Omega_rad_today = 5.4e-5  # Radiation today (negligible)
T_CMB = 2.725  # K
a_eq = Omega_rad_today / Omega_m  # Scale factor at equality
z_eq = 1.0/a_eq - 1.0
T_eq = T_CMB * (1 + z_eq)

print("MATTER-RADIATION EQUALITY:")
print(f"  Redshift z_eq               ≈ {z_eq:.0f}")
print(f"  Temperature T_eq            ≈ {T_eq:.0f} K ≈ {T_eq*k_B/e:.2f} eV")
print(f"  At this point: ρ_m = ρ_rad (thermodynamic transition)")
print(f"  Before: radiation-dominated (P = ρ/3, w = 1/3)")
print(f"  After: matter-dominated (P ≈ 0, w = 0)")
print()

# Big Bang Nucleosynthesis
T_BBN_MeV = 0.8  # MeV
T_BBN_K = T_BBN_MeV * 1e6 * e / k_B

print("BIG BANG NUCLEOSYNTHESIS:")
print(f"  Temperature T_BBN           ≈ {T_BBN_MeV:.1f} MeV ≈ {T_BBN_K:.3e} K")
print(f"  At this T, n/p ratio freezes → determines He-4 abundance")
print(f"  Baryon fraction Ω_b set by:")
print(f"    • Baryon-antibaryon asymmetry (baryogenesis)")
print(f"    • Entropy per baryon (thermal relic)")
print()

# Dark matter freeze-out
print("DARK MATTER FREEZE-OUT:")
print(f"  WIMPs freeze out when Γ_annihilation < H(T)")
print(f"  Typical freeze-out: T_fo ~ m_WIMP/20")
print(f"  For m_WIMP ~ 100 GeV: T_fo ~ 5 GeV ~ 6×10¹³ K")
print(f"  Thermal relic density Ω_DM determined by:")
print(f"    Ω_DM h² ∝ 1/⟨σv⟩ (inverse of annihilation cross section)")
print()

print("THERMODYNAMIC PHASE TRANSITIONS:")
print(f"  1. Electroweak (T ~ 160 GeV): W, Z, Higgs acquire mass")
print(f"  2. QCD (T ~ 170 MeV): Quarks → hadrons (confinement)")
print(f"  3. Recombination (T ~ 0.26 eV): e + p → H (CMB released)")
print()
print(f"  Each transition changes equation of state w = P/ρ")
print(f"  and entropy content — determining Ω_m.")
print()

print("=" * 80)
print("CALIBRATION:")
print(f"Planck 2018: Ω_m = 0.3153 ± 0.0073")
print(f"             Ω_b = 0.0493 ± 0.0006")
print(f"TriPhase uses observed values (thermodynamic relics).")
print("=" * 80)

input("Press Enter to exit...")
