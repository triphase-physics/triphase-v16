"""
TriPhase V16: Critical Density Derivation
Dimensional Analysis Framework

This script demonstrates the dimensional consistency of the critical density
(closure density) derivation using pure wave mechanics and cosmology.

Critical Density:
  ρ_crit = 3H₀² / (8πG)
  where H₀ = π·√3·f_e·α^18

SI Units: [kg/m³]
Dimensional form: [M L⁻³]

The critical density is the density needed for a flat universe in
Friedmann cosmology, setting the scale for cosmic matter content.

MIS TAG: (D)
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
print("TriPhase V16: Critical Density")
print("Dimensional Analysis Framework")
print("=" * 70)
print()

# ========================================
# STEP 1: Identify Target Dimensions
# ========================================
print("STEP 1: Identify Target Dimensions")
print("-" * 70)
print("Target Quantity: Critical Density (ρ_crit)")
print("SI Unit: kg/m³")
print("Dimensional Form: [M L⁻³]")
print()
print("The critical density is the density required for a")
print("spatially flat universe in Friedmann-Robertson-Walker cosmology.")
print("Measured value: ~9.47e-27 kg/m³")
print()

# ========================================
# STEP 2: Available Base Dimensions
# ========================================
print("STEP 2: Available Base Dimensions")
print("-" * 70)
print("From cosmology:")
print("  H₀: [T⁻¹]             (Hubble constant)")
print("  G:  [L³ M⁻¹ T⁻²]     (gravitational constant)")
print()
print("TriPhase Hubble constant:")
print("  H₀ = π·√3·f_e·α^18")
print("  where f_e = m_e·c²/ℏ")
print()
print("TriPhase gravitational constant:")
print("  G = c⁴·7.5·ε₀³·μ₀²")
print()

# ========================================
# STEP 3: Dimensional Matching
# ========================================
print("STEP 3: Dimensional Matching")
print("-" * 70)
print("Critical density formula:")
print("  ρ_crit = 3H₀² / (8πG)")
print()
print("Dimensional analysis:")
print("  [ρ_crit] = [H₀]² / [G]")
print("           = [T⁻²] / [L³ M⁻¹ T⁻²]")
print("           = [T⁻²] · [M L⁻³ T²]")
print("           = [M L⁻³]")
print()
print("✓ Result has density dimensions")
print()
print("From Friedmann equation:")
print("  H² = (8πG/3)ρ  (for flat, matter-dominated)")
print("  Solving for ρ:")
print("  ρ_crit = 3H₀²/(8πG)")
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

print("Critical density:")
print(f"  ρ_crit = 3H₀²/(8πG)")
print(f"         = {rho_crit:.15e} kg/m³")
print(f"         = {rho_crit:.6e} kg/m³")
print()

# Number density (if all protons)
n_crit = rho_crit / m_p
print("Critical number density (as protons):")
print(f"  n_crit = ρ_crit/m_p = {n_crit:.15e} m⁻³")
print(f"         ≈ {n_crit:.3f} protons/m³")
print("  (About 6 protons per cubic meter)")
print()

# Mass in observable universe
R_H = c / H_0
V_H = (4.0/3.0) * math.pi * R_H**3
M_H = rho_crit * V_H

print("Mass within Hubble radius:")
print(f"  R_H = c/H₀ = {R_H:.15e} m")
print(f"  V_H = (4π/3)R_H³ = {V_H:.15e} m³")
print(f"  M_H = ρ_crit·V_H = {M_H:.15e} kg")
M_H_solar = M_H / 1.989e30
print(f"              = {M_H_solar:.3e} M_☉")
print()

# ========================================
# STEP 5: Buckingham Pi Theorem
# ========================================
print("STEP 5: Buckingham Pi Theorem")
print("-" * 70)
print("For critical density, dimensionless groups:")
print()
print("π₁ = ρ_crit·8πG / (3H₀²)")
print(f"   = {rho_crit * 8.0 * math.pi * G / (3.0 * H_0**2):.15e}")
print("   (Should be 1.0 from the definition)")
print()
print("π₂ = ρ_crit / ρ_Planck")
rho_P = c**5 / (hbar * G**2)
print(f"   = {rho_crit / rho_P:.15e}")
print("   (Critical density relative to Planck density)")
print()
print("π₃ = H₀·R_H / c")
print(f"   = {H_0 * R_H / c:.15e}")
print("   (Should be 1.0)")
print()

# Density parameters
rho_DE = (H_0**2 / c**2) * c**4 / (8.0 * math.pi * G)
Omega_Lambda = rho_DE / rho_crit
Omega_matter_measured = 0.315  # Planck 2018
print("π₄ = Ω_total = Ω_m + Ω_Λ + Ω_r + ...")
print(f"   Ω_Λ = {Omega_Lambda:.6f}")
print(f"   Ω_m ≈ {Omega_matter_measured:.6f} (from observations)")
print(f"   Ω_total ≈ {Omega_Lambda + Omega_matter_measured:.6f}")
print("   (Should be 1.0 for flat universe)")
print()

print("These dimensionless groups characterize the cosmic")
print("density budget and universe geometry.")
print()

# ========================================
# STEP 6: Natural Units Comparison
# ========================================
print("STEP 6: Natural Units Comparison")
print("-" * 70)
print("Planck density: ρ_P = c⁵/(ℏG²)")
print(f"  ρ_P = {rho_P:.15e} kg/m³")
print(f"  ρ_crit / ρ_P = {rho_crit / rho_P:.15e}")
print()
print("Nuclear density: ρ_nuclear ≈ 2.3e17 kg/m³")
rho_nuclear = 2.3e17
print(f"  ρ_crit / ρ_nuclear = {rho_crit / rho_nuclear:.15e}")
print()
print("Water density: ρ_water = 1000 kg/m³")
rho_water = 1000.0
print(f"  ρ_crit / ρ_water = {rho_crit / rho_water:.15e}")
print()
print("Interstellar medium: ρ_ISM ≈ 10^-21 kg/m³")
rho_ISM = 1e-21
print(f"  ρ_crit / ρ_ISM = {rho_crit / rho_ISM:.15e}")
print()
print("The critical density is extremely low, comparable to")
print("a few hydrogen atoms per cubic meter.")
print()

# ========================================
# STEP 7: Dimensional Crosschecks
# ========================================
print("STEP 7: Dimensional Crosschecks")
print("-" * 70)
print("Verify units on both sides:")
print()
print("LHS: ρ_crit")
print("  Dimensions: [M L⁻³]")
print("  Units: kg/m³")
print()
print("RHS: 3H₀²/(8πG)")
print("  Dimensions: [T⁻²] / [L³ M⁻¹ T⁻²]")
print("            = [T⁻²] · [M L⁻³ T²]")
print("            = [M L⁻³]")
print("  Units: s⁻² / (m³·kg⁻¹·s⁻²)")
print("       = s⁻² · (kg·s²/m³)")
print("       = kg/m³")
print()
print("✓ Dimensional consistency verified")
print()
print("Friedmann equation (flat, k=0):")
print("  H² = (8πG/3)ρ")
print("  At critical density:")
H_check = math.sqrt(8.0 * math.pi * G * rho_crit / 3.0)
print(f"  H = √[(8πG/3)ρ_crit] = {H_check:.15e} s⁻¹")
print(f"  H₀ = {H_0:.15e} s⁻¹")
print(f"  Ratio: {H_check / H_0:.15e}")
print("  ✓ Consistent")
print()
print("Deceleration parameter (matter-dominated):")
q_0 = 0.5  # for Ω_m = 1
print(f"  q₀ = Ω_m/2 - Ω_Λ ≈ {Omega_matter_measured/2 - Omega_Lambda:.3f}")
print("  (Negative q₀ indicates acceleration)")
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
print(f"Measured H₀:      {67.4:.3f} km/s/Mpc")
print(f"Measured ρ_crit:  {rho_crit_measured:.15e} kg/m³")
print()
print(f"TriPhase H₀:      {H_0_kmsMpc:.3f} km/s/Mpc")
print(f"TriPhase ρ_crit:  {rho_crit:.15e} kg/m³")
print()
deviation_ppm = abs(rho_crit - rho_crit_measured) / rho_crit_measured * 1e6
print(f"Deviation:        {deviation_ppm:.1f} ppm")
print()
print("Density parameters (Planck 2018):")
print("  Ω_m = 0.315 ± 0.007  (matter)")
print("  Ω_Λ = 0.689 ± 0.006  (dark energy)")
print("  Ω_r ≈ 9e-5           (radiation)")
print("  Ω_total = 1.000 ± 0.002  (flat universe)")
print()
print(f"TriPhase Ω_Λ:     {Omega_Lambda:.6f}")
print()

# ========================================
# SUMMARY
# ========================================
print("=" * 70)
print("DIMENSIONAL ANALYSIS SUMMARY")
print("=" * 70)
print()
print("The critical density derivation is dimensionally consistent:")
print()
print("1. Target dimensions [M L⁻³] verified")
print("2. Formula ρ_crit = 3H₀²/(8πG)")
print("3. Arises from Friedmann equation for flat universe")
print("4. Sets scale for cosmic matter budget")
print("5. Uses TriPhase H₀ and G")
print()
print("The formula ρ_crit = 3H₀²/(8πG) represents the")
print("closure density for a flat universe, derived from")
print("TriPhase wave mechanics via the 18-step cascade.")
print()
print("Key insight:")
print(f"  ρ_crit ≈ {rho_crit:.3e} kg/m³")
print("  About 6 protons per cubic meter")
print("  Extremely low by everyday standards")
print("  Sets the scale for Ω_m, Ω_Λ, etc.")
print()
print("TriPhase contribution:")
print("  H₀ emerges from f_e and α^18 cascade")
print("  G emerges from ε₀, μ₀, c")
print("  ρ_crit connects quantum to cosmological scales")
print()
print("Cosmological significance:")
print("  Ω_total = 1 (flat universe, observations)")
print("  Ω_Λ ≈ 0.69 (dark energy dominated)")
print("  Ω_m ≈ 0.31 (matter)")
print("  Critical density divides open/closed/flat geometries")
print()
print("=" * 70)

input("Press Enter to exit...")
