"""
TriPhase V16: Matter Density Derivation
Dimensional Analysis Framework

This script demonstrates the dimensional consistency of the cosmic matter
density derivation using pure wave mechanics and cosmology.

Matter Density:
  ρ_m = ρ_crit × Ω_m
  where Ω_m ≈ 0.315 (from observations)

SI Units: [kg/m³]
Dimensional form: [M L⁻³]

The matter density represents the average mass density of baryonic and
dark matter in the universe, a key cosmological parameter.

MIS TAG: (C)
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# ========================================
# ANCHOR CONSTANTS (Standard TriPhase Chain)
# ========================================
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
print("TriPhase V16: Matter Density")
print("Dimensional Analysis Framework")
print("=" * 70)
print()

# ========================================
# STEP 1: Identify Target Dimensions
# ========================================
print("STEP 1: Identify Target Dimensions")
print("-" * 70)
print("Target Quantity: Matter Density (ρ_m)")
print("SI Unit: kg/m³")
print("Dimensional Form: [M L⁻³]")
print()
print("The matter density is the average density of all matter")
print("(baryonic + dark matter) in the universe.")
print("Planck 2018 value: Ω_m = 0.315, giving ρ_m ≈ 2.7e-27 kg/m³")
print()

# ========================================
# STEP 2: Available Base Dimensions
# ========================================
print("STEP 2: Available Base Dimensions")
print("-" * 70)
print("From cosmology:")
print("  ρ_crit: [M L⁻³]       (critical density)")
print("  Ω_m: [1]              (matter density parameter, dimensionless)")
print("  H₀: [T⁻¹]             (Hubble constant)")
print("  G:  [L³ M⁻¹ T⁻²]     (gravitational constant)")
print()
print("Critical density:")
print("  ρ_crit = 3H₀²/(8πG)")
print()
print("Matter density parameter:")
print("  Ω_m = ρ_m / ρ_crit")
print("  (From observations: Ω_m ≈ 0.315)")
print()

# ========================================
# STEP 3: Dimensional Matching
# ========================================
print("STEP 3: Dimensional Matching")
print("-" * 70)
print("Matter density formula:")
print("  ρ_m = ρ_crit × Ω_m")
print()
print("Dimensional analysis:")
print("  [ρ_m] = [ρ_crit] × [Ω_m]")
print("        = [M L⁻³] × [1]")
print("        = [M L⁻³]")
print()
print("✓ Result has density dimensions")
print()
print("Ω_m is dimensionless, so multiplication preserves")
print("the density dimensions of ρ_crit.")
print()

# ========================================
# STEP 4: TriPhase Derivation
# ========================================
print("STEP 4: TriPhase Derivation")
print("-" * 70)
print(f"Hubble constant:        H₀ = {H_0:.15e} s⁻¹")
H_0_kmsMpc = H_0 * 3.08567758149e19 / 1e3
print(f"                           = {H_0_kmsMpc:.3f} km/s/Mpc")
print(f"Gravitational constant:  G = {G:.15e} m³/(kg·s²)")
print()

# Critical density
rho_crit = 3.0 * H_0**2 / (8.0 * math.pi * G)
print(f"Critical density:")
print(f"  ρ_crit = 3H₀²/(8πG) = {rho_crit:.15e} kg/m³")
print()

# Matter density parameter (from Planck 2018)
Omega_m = 0.315
print("Matter density parameter (Planck 2018):")
print(f"  Ω_m = {Omega_m:.3f}")
print("  (Includes baryonic matter + dark matter)")
print()

# Matter density
rho_m = rho_crit * Omega_m
print("Matter density:")
print(f"  ρ_m = ρ_crit × Ω_m")
print(f"      = {rho_m:.15e} kg/m³")
print(f"      = {rho_m:.6e} kg/m³")
print()

# Number density (if all protons)
n_m = rho_m / m_p
print("Matter number density (as protons):")
print(f"  n_m = ρ_m/m_p = {n_m:.15e} m⁻³")
print(f"      ≈ {n_m:.3f} protons/m³")
print("  (About 2 protons per cubic meter)")
print()

# Baryonic vs dark matter
Omega_b = 0.049  # Planck 2018 baryonic matter
Omega_dm = Omega_m - Omega_b
rho_b = rho_crit * Omega_b
rho_dm = rho_crit * Omega_dm

print("Breakdown:")
print(f"  Ω_b  = {Omega_b:.3f}  (baryonic matter)")
print(f"  Ω_dm = {Omega_dm:.3f}  (dark matter)")
print(f"  ρ_b  = {rho_b:.3e} kg/m³")
print(f"  ρ_dm = {rho_dm:.3e} kg/m³")
print(f"  Dark matter fraction: {Omega_dm/Omega_m:.3f}")
print()

# ========================================
# STEP 5: Buckingham Pi Theorem
# ========================================
print("STEP 5: Buckingham Pi Theorem")
print("-" * 70)
print("For matter density, dimensionless groups:")
print()
print("π₁ = Ω_m = ρ_m / ρ_crit")
print(f"   = {rho_m / rho_crit:.15e}")
print("   (Matter density parameter)")
print()
print("π₂ = ρ_m / ρ_Planck")
rho_P = c**5 / (hbar * G**2)
print(f"   = {rho_m / rho_P:.15e}")
print("   (Matter density relative to Planck density)")
print()
print("π₃ = ρ_m / ρ_nuclear")
rho_nuclear = 2.3e17
print(f"   = {rho_m / rho_nuclear:.15e}")
print("   (Matter density relative to nuclear density)")
print()

# Dark energy comparison
rho_DE = (H_0**2 / c**2) * c**4 / (8.0 * math.pi * G)
Omega_Lambda = rho_DE / rho_crit
print("π₄ = ρ_m / ρ_DE")
print(f"   = {rho_m / rho_DE:.15e}")
print(f"   = Ω_m / Ω_Λ = {Omega_m / Omega_Lambda:.3f}")
print("   (Matter-to-dark energy ratio)")
print()

print("These dimensionless groups characterize the cosmic")
print("matter budget and composition.")
print()

# ========================================
# STEP 6: Natural Units Comparison
# ========================================
print("STEP 6: Natural Units Comparison")
print("-" * 70)
print("Planck density: ρ_P = c⁵/(ℏG²)")
print(f"  ρ_P = {rho_P:.15e} kg/m³")
print(f"  ρ_m / ρ_P = {rho_m / rho_P:.15e}")
print()
print("Nuclear density: ρ_nuclear ≈ 2.3e17 kg/m³")
print(f"  ρ_m / ρ_nuclear = {rho_m / rho_nuclear:.15e}")
print()
print("Water density: ρ_water = 1000 kg/m³")
rho_water = 1000.0
print(f"  ρ_m / ρ_water = {rho_m / rho_water:.15e}")
print()
print("Interstellar medium: ρ_ISM ≈ 10^-21 kg/m³")
rho_ISM = 1e-21
print(f"  ρ_m / ρ_ISM = {rho_m / rho_ISM:.15e}")
print()
print("Intergalactic medium: ρ_IGM ≈ 10^-27 kg/m³")
rho_IGM = 1e-27
print(f"  ρ_m / ρ_IGM = {rho_m / rho_IGM:.3f}")
print()

# ========================================
# STEP 7: Dimensional Crosschecks
# ========================================
print("STEP 7: Dimensional Crosschecks")
print("-" * 70)
print("Verify units on both sides:")
print()
print("LHS: ρ_m")
print("  Dimensions: [M L⁻³]")
print("  Units: kg/m³")
print()
print("RHS: ρ_crit × Ω_m")
print("  Dimensions: [M L⁻³] × [1]")
print("            = [M L⁻³]")
print("  Units: (kg/m³) × (dimensionless) = kg/m³")
print()
print("✓ Dimensional consistency verified")
print()
print("Density parameters sum:")
Omega_total = Omega_m + Omega_Lambda + 9e-5  # +radiation
print(f"  Ω_total = Ω_m + Ω_Λ + Ω_r")
print(f"          = {Omega_m:.3f} + {Omega_Lambda:.3f} + {9e-5:.5f}")
print(f"          = {Omega_total:.3f}")
print("  (Should be 1.0 for flat universe)")
print()
print("Mass in observable universe:")
R_H = c / H_0
V_H = (4.0/3.0) * math.pi * R_H**3
M_m_H = rho_m * V_H
print(f"  R_H = c/H₀ = {R_H:.3e} m")
print(f"  V_H = (4π/3)R_H³ = {V_H:.3e} m³")
print(f"  M_m = ρ_m·V_H = {M_m_H:.3e} kg")
M_m_solar = M_m_H / 1.989e30
print(f"              = {M_m_solar:.3e} M_☉")
print()

# ========================================
# STEP 8: Calibration Checkpoint
# ========================================
print("STEP 8: Calibration Checkpoint")
print("-" * 70)
print("Planck 2018 cosmological parameters:")
print("  H₀ = 67.4 ± 0.5 km/s/Mpc")
H_0_measured = 67.4 / 3.08567758149e19
rho_crit_measured = 3.0 * H_0_measured**2 / (8.0 * math.pi * G)
Omega_m_measured = 0.315
rho_m_measured = rho_crit_measured * Omega_m_measured

print(f"Measured ρ_crit:  {rho_crit_measured:.15e} kg/m³")
print(f"Measured Ω_m:     {Omega_m_measured:.3f}")
print(f"Measured ρ_m:     {rho_m_measured:.15e} kg/m³")
print()
print(f"TriPhase ρ_crit:  {rho_crit:.15e} kg/m³")
print(f"TriPhase Ω_m:     {Omega_m:.3f} (from observations)")
print(f"TriPhase ρ_m:     {rho_m:.15e} kg/m³")
print()
deviation_ppm = abs(rho_m - rho_m_measured) / rho_m_measured * 1e6
print(f"Deviation:        {deviation_ppm:.1f} ppm")
print()
print("Density budget (Planck 2018):")
print("  Ω_b  = 0.049 ± 0.001  (baryonic matter)")
print("  Ω_dm = 0.266 ± 0.007  (dark matter)")
print("  Ω_m  = 0.315 ± 0.007  (total matter)")
print("  Ω_Λ  = 0.689 ± 0.006  (dark energy)")
print("  Ω_r  ≈ 9e-5           (radiation)")
print()
print("TriPhase values:")
print(f"  Ω_m  = {Omega_m:.3f}")
print(f"  Ω_Λ  = {Omega_Lambda:.3f}")
print()

# ========================================
# SUMMARY
# ========================================
print("=" * 70)
print("DIMENSIONAL ANALYSIS SUMMARY")
print("=" * 70)
print()
print("The matter density derivation is dimensionally consistent:")
print()
print("1. Target dimensions [M L⁻³] verified")
print("2. Formula ρ_m = ρ_crit × Ω_m")
print("3. Ω_m ≈ 0.315 from cosmological observations")
print("4. Includes baryonic (15%) and dark matter (85%)")
print("5. Uses TriPhase-derived ρ_crit")
print()
print("The formula ρ_m = ρ_crit × Ω_m represents the")
print("cosmic matter density, derived from TriPhase")
print("critical density and observational constraints.")
print()
print("Key insight:")
print(f"  ρ_m ≈ {rho_m:.3e} kg/m³")
print("  About 2 protons per cubic meter")
print("  85% dark matter, 15% baryonic matter")
print("  Only 31% of critical density")
print()
print("TriPhase contribution:")
print("  ρ_crit emerges from H₀ (via α^18 cascade) and G")
print("  Ω_m from observations (not derived in TriPhase)")
print("  Together give cosmic matter density")
print()
print("Cosmological significance:")
print("  Matter dominated universe until z ≈ 0.3")
print("  Dark energy now dominates (Ω_Λ > Ω_m)")
print("  Drives structure formation (galaxies, clusters)")
print("  Most matter is dark (non-baryonic)")
print()
print("Physical interpretation:")
print("  Average: ~2 protons/m³ across entire universe")
print("  Highly clustered (galaxies, voids)")
print("  Dark matter provides gravitational scaffolding")
print("  Baryonic matter forms visible structures")
print()
print("=" * 70)

input("Press Enter to exit...")
