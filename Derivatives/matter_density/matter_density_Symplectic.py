"""
TriPhase V16 — Matter Density (Symplectic Framework)

Copyright © 2025 MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
MIS TAG: (C)

SYMPLECTIC INTERPRETATION:
Matter density ρ_m = ρ_crit × Ω_m is the non-relativistic mass-energy density
in the cosmological phase space (a, ȧ). In the Friedmann-Robertson-Walker (FRW)
metric, matter with equation of state w = 0 (P_m ≈ 0 for cold dark matter and
baryons) contributes to the Hamiltonian constraint as ρ_m(t) = ρ_{m,0} a(t)^{-3},
diluting as the universe expands.

In the symplectic formulation, matter couples to the scale factor a through the
canonical momentum p_a. The density parameter Ω_m = ρ_m/ρ_crit determines the
deceleration of the expansion: matter creates attractive gravity, slowing ȧ.
The phase space trajectory in the (a, ȧ) plane follows a power law a(t) ~ t^{2/3}
during matter domination, transitioning to exponential a(t) ~ e^{Ht} when dark
energy (Ω_Λ) dominates.
"""

import math

# ========== ANCHOR CHAIN (VERBATIM) ==========
epsilon_0 = 8.8541878128e-12   # F/m
mu_0      = 1.25663706212e-6   # H/m
e         = 1.602176634e-19     # C (exact, SI 2019)
c         = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0       = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha     = 1.0 / alpha_inv
hbar      = Z_0 * e**2 / (4.0 * math.pi * alpha)
h         = 2.0 * math.pi * hbar
G         = c**4 * 7.5 * epsilon_0**3 * mu_0**2
r_e       = 2.8179403262e-15   # m (classical electron radius)
m_e       = hbar * alpha / (c * r_e)
f_e       = m_e * c**2 / hbar
T_17      = 17 * 18 // 2
mp_me     = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p       = m_e * mp_me
H_0       = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r      = c**4 / (8.0 * math.pi * G)

print("=" * 70)
print("TRIPHASE V16 — MATTER DENSITY (SYMPLECTIC FRAMEWORK)")
print("=" * 70)
print()

# ========== SYMPLECTIC DERIVATION ==========
print("PHASE SPACE STRUCTURE:")
print("  Cosmological phase space: (a, ȧ) where a(t) = scale factor")
print("  Matter contribution: ρ_m(a) = ρ_{m,0} a^{-3} (conservation of mass)")
print("  Conjugate momentum: p_a = -3ȧa²/(8πG)")
print("  Symplectic 2-form: ω = dp_a ∧ da")
print()

print("HAMILTONIAN FORMULATION:")
print("  Friedmann equation: H² = (ȧ/a)² = (8πG/3)[ρ_m(a) + ρ_r(a) + ρ_Λ]")
print("  Matter (w=0): ρ_m ∝ a^{-3}, P_m ≈ 0 (non-relativistic)")
print("  Radiation (w=1/3): ρ_r ∝ a^{-4}, P_r = ρ_r c²/3")
print("  Dark energy (w=-1): ρ_Λ = const, P_Λ = -ρ_Λ c²")
print("  Density parameters: Ω_i = ρ_i / ρ_crit")
print("  Matter density: ρ_m = ρ_crit × Ω_m")
print()

print("SYMPLECTIC INVARIANT:")
print("  Energy conservation: ρ̇ + 3H(ρ + P/c²) = 0 (1st law of thermodynamics)")
print("  For matter (P_m ≈ 0): ρ_m a³ = const")
print("  Deceleration parameter: q = -ä/(aH²) = Ω_m/2 - Ω_Λ (matter-dominated)")
print("  Phase space trajectory: a(t) ~ t^{2/3} when Ω_m ≈ 1")
print()

print("TRIPHASE DERIVATION:")
print("  Formula: ρ_m = ρ_crit × Ω_m")
print("  Where: ρ_crit = 3 H_0² / (8π G)")
print("         Ω_m = 0.315 (Planck 2018 total matter fraction)")
print("  Breakdown: Ω_b = 0.049 (baryonic matter)")
print("             Ω_CDM = 0.266 (cold dark matter)")
print()

# Compute critical density
rho_crit = 3.0 * H_0**2 / (8.0 * math.pi * G)
print(f"  H_0      = {H_0:.6e} Hz")
print(f"  ρ_crit   = {rho_crit:.6e} kg/m³")
print()

# Compute matter density (Planck 2018 values)
Omega_m = 0.315  # Total matter (baryonic + dark)
Omega_b = 0.049  # Baryonic matter only
Omega_CDM = Omega_m - Omega_b  # Cold dark matter
rho_m = rho_crit * Omega_m
rho_b = rho_crit * Omega_b
rho_CDM = rho_crit * Omega_CDM

print(f"  Matter density parameters (Planck 2018):")
print(f"  Ω_m   = {Omega_m:.3f} (total matter)")
print(f"  Ω_b   = {Omega_b:.3f} (baryonic matter)")
print(f"  Ω_CDM = {Omega_CDM:.3f} (cold dark matter)")
print()
print(f"  ρ_m   = ρ_crit × Ω_m   = {rho_m:.6e} kg/m³")
print(f"  ρ_b   = ρ_crit × Ω_b   = {rho_b:.6e} kg/m³")
print(f"  ρ_CDM = ρ_crit × Ω_CDM = {rho_CDM:.6e} kg/m³")
print()

# Compute number densities
n_proton = rho_b / m_p  # Approximate: all baryonic matter as protons
n_hydrogen = n_proton  # Simplified: all baryons as hydrogen atoms
print(f"  Baryonic number density (assuming all protons):")
print(f"  m_p = {m_p:.6e} kg")
print(f"  n_p = ρ_b / m_p = {n_proton:.6e} protons/m³")
print(f"                  = {n_proton * 1e-6:.6f} protons/cm³")
print()

# Compare to local density
rho_solar = 1.4e3  # kg/m³ (solar average density)
rho_ISM = 1.67e-21  # kg/m³ (interstellar medium, ~1 H atom/cm³)
print(f"  Comparison to local densities:")
print(f"  ρ_m (cosmic average) = {rho_m:.6e} kg/m³")
print(f"  ρ_ISM (interstellar)  = {rho_ISM:.6e} kg/m³")
print(f"  ρ_solar (Sun average) = {rho_solar:.6e} kg/m³")
print(f"  ρ_m / ρ_ISM          = {rho_m / rho_ISM:.3f}")
print()

# ========== CALIBRATION CHECKPOINT ==========
print("CALIBRATION CHECKPOINT:")
print("  Using Planck 2018 + H_0 = 67.4 km/s/Mpc:")
print("    ρ_crit = 8.6e-27 kg/m³")
print("    Ω_m = 0.315")
print("    ρ_m = 2.71e-27 kg/m³")
print()
rho_crit_measured = 8.6e-27
Omega_m_measured = 0.315
rho_m_measured = rho_crit_measured * Omega_m_measured
print(f"  Measured ρ_crit = {rho_crit_measured:.6e} kg/m³")
print(f"  TriPhase ρ_crit = {rho_crit:.6e} kg/m³")
deviation_crit_ppm = abs(rho_crit - rho_crit_measured) / rho_crit_measured * 1e6
print(f"  Deviation: {deviation_crit_ppm:.1f} ppm")
print()
print(f"  Measured ρ_m    = {rho_m_measured:.6e} kg/m³")
print(f"  TriPhase ρ_m    = {rho_m:.6e} kg/m³")
deviation_m_ppm = abs(rho_m - rho_m_measured) / rho_m_measured * 1e6
print(f"  Deviation: {deviation_m_ppm:.1f} ppm")
print()

# Show evolution with scale factor
print(f"  Evolution with scale factor:")
print(f"  ρ_m(a) = ρ_{m,0} a^{-3}")
print(f"  At recombination (z=1100, a=1/1101): ρ_m = {rho_m * 1101**3:.6e} kg/m³")
print(f"  At matter-radiation equality (z≈3400): ρ_m = {rho_m * 3401**3:.6e} kg/m³")
print()

# ========== SYMPLECTIC GEOMETRY INSIGHT ==========
print("SYMPLECTIC GEOMETRY INSIGHT:")
print("  Matter density ρ_m determines the curvature of cosmological phase space.")
print("  In the (a, ȧ) plane, matter creates a decelerating symplectic flow:")
print("  attractive gravity slows the expansion. The trajectory a(t) ~ t^{2/3}")
print("  is a power-law orbit in phase space, characteristic of Ω_m ≈ 1.")
print()
print("  The matter-dominated era (0.3 < z < 3400) is when this symplectic")
print("  structure dominated the universe's evolution. Structure formation")
print("  (galaxies, clusters) occurred via gravitational instability—density")
print("  perturbations growing exponentially in the matter phase space.")
print()
print("  Today, ρ_m ≈ 2.7e-27 kg/m³ is subdominant to dark energy ρ_Λ ≈ 5.9e-27 kg/m³.")
print("  The phase space has transitioned from deceleration (matter-dominated)")
print("  to acceleration (dark energy-dominated), crossing the critical point")
print("  at z ≈ 0.67 where ä changed sign. This is a topological transition in")
print("  the symplectic flow—a bifurcation in the cosmological phase space.")
print("=" * 70)

input("Press Enter to exit...")
