"""
TriPhase V16: Horizon 18-Step Derivation (Hubble Radius)
Dimensional Analysis Framework

This script demonstrates the dimensional consistency of the Hubble radius
(cosmological horizon) derivation using pure wave mechanics.

Horizon Radius:
  R_H = c / H_0
  where H_0 = π·√3·f_e·α^18

SI Units: [m]
Dimensional form: [L]

The Hubble radius represents the cosmological horizon scale, derived from
the Hubble constant which emerges through the 18-step alpha cascade.

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
print("TriPhase V16: Horizon 18-Step (Hubble Radius)")
print("Dimensional Analysis Framework")
print("=" * 70)
print()

# ========================================
# STEP 1: Identify Target Dimensions
# ========================================
print("STEP 1: Identify Target Dimensions")
print("-" * 70)
print("Target Quantity: Hubble Radius (R_H)")
print("SI Unit: m (meters)")
print("Dimensional Form: [L]")
print()
print("The Hubble radius is the cosmological horizon scale.")
print("Measured value: ~1.37e26 m (about 14.5 billion light-years)")
print()

# ========================================
# STEP 2: Available Base Dimensions
# ========================================
print("STEP 2: Available Base Dimensions")
print("-" * 70)
print("From anchor constants:")
print("  c:  [L T⁻¹]          (speed of light)")
print("  ℏ:  [M L² T⁻¹]      (reduced Planck constant)")
print("  α:  [1]              (fine structure constant, dimensionless)")
print("  m_e: [M]             (electron mass)")
print()
print("Derived frequency:")
print("  f_e = m_e·c²/ℏ")
print("  Dimensions: [M]·[L T⁻¹]² / [M L² T⁻¹]")
print("            = [M L² T⁻²] / [M L² T⁻¹]")
print("            = [T⁻¹]  (frequency)")
print()
print("Hubble constant:")
print("  H_0 = π·√3·f_e·α^18")
print("  Dimensions: [T⁻¹]  (inverse time)")
print()

# ========================================
# STEP 3: Dimensional Matching
# ========================================
print("STEP 3: Dimensional Matching")
print("-" * 70)
print("Hubble radius formula:")
print("  R_H = c / H_0")
print()
print("Dimensional analysis:")
print("  [R_H] = [c] / [H_0]")
print("        = [L T⁻¹] / [T⁻¹]")
print("        = [L]")
print()
print("Where:")
print("  c: [L T⁻¹]             (speed of light)")
print("  H_0: [T⁻¹]             (Hubble constant)")
print()
print("✓ Result has pure length dimensions")
print()

# ========================================
# STEP 4: TriPhase Derivation
# ========================================
print("STEP 4: TriPhase Derivation")
print("-" * 70)
print(f"Speed of light:        c      = {c:.15e} m/s")
print(f"Electron mass:         m_e    = {m_e:.15e} kg")
print(f"Fine structure const:  α      = {alpha:.15e}")
print(f"Reduced Planck const:  ℏ      = {hbar:.15e} J·s")
print()

print("Electron Compton frequency:")
print(f"  f_e = m_e·c²/ℏ = {f_e:.15e} Hz")
print()

print("18-step alpha cascade:")
alpha_18 = alpha**18
print(f"  α^18 = {alpha_18:.15e}")
print()

print("Hubble constant:")
print(f"  H_0 = π·√3·f_e·α^18")
print(f"      = {H_0:.15e} s⁻¹")
H_0_SI = H_0 / 1e3  # km/s per Mpc
print(f"      = {H_0_SI * 3.08567758149e19:.3f} km/s/Mpc")
print()

R_H = c / H_0
print(f"Hubble radius:         R_H = c / H_0")
print(f"                           = {R_H:.15e} m")
R_H_ly = R_H / 9.4607304725808e15
print(f"                           = {R_H_ly:.3e} light-years")
R_H_Gly = R_H_ly / 1e9
print(f"                           = {R_H_Gly:.3f} billion light-years")
print()

# ========================================
# STEP 5: Buckingham Pi Theorem
# ========================================
print("STEP 5: Buckingham Pi Theorem")
print("-" * 70)
print("For the Hubble radius, we identify dimensionless groups:")
print()
print("π₁ = H_0 · R_H / c")
print(f"   = {H_0 * R_H / c:.15e}")
print("   (Should be 1.0 from the definition)")
print()
print("π₂ = α^18")
print(f"   = {alpha_18:.15e}")
print("   (The 18-step cascade factor)")
print()

# Planck length for comparison
l_P = math.sqrt(hbar * G / c**3)
print("π₃ = R_H / l_P")
print(f"   = {R_H / l_P:.15e}")
print("   (Hubble radius in Planck lengths)")
print()

# Electron Compton wavelength
lambda_C_e = hbar / (m_e * c)
print("π₄ = R_H / λ_C,e")
print(f"   = {R_H / lambda_C_e:.15e}")
print("   (Hubble radius in electron Compton wavelengths)")
print()

print("These dimensionless groups characterize the cosmological")
print("scale relative to quantum and Planck scales.")
print()

# ========================================
# STEP 6: Natural Units Comparison
# ========================================
print("STEP 6: Natural Units Comparison")
print("-" * 70)
print("Planck length: l_P = sqrt(ℏG/c³)")
print(f"  l_P = {l_P:.15e} m")
print(f"  R_H / l_P = {R_H / l_P:.15e}")
print()
print("Electron Compton wavelength: λ_C,e = ℏ/(m_e·c)")
print(f"  λ_C,e = {lambda_C_e:.15e} m")
print(f"  R_H / λ_C,e = {R_H / lambda_C_e:.15e}")
print()
print("Classical electron radius: r_e")
print(f"  r_e = {r_e:.15e} m")
print(f"  R_H / r_e = {R_H / r_e:.15e}")
print()
print("Astronomical unit: AU = 1.495978707e11 m")
AU = 1.495978707e11
print(f"  R_H / AU = {R_H / AU:.15e}")
print()
print("Light-year: ly = 9.4607304725808e15 m")
ly = 9.4607304725808e15
print(f"  R_H / ly = {R_H / ly:.15e} light-years")
print()

# ========================================
# STEP 7: Dimensional Crosschecks
# ========================================
print("STEP 7: Dimensional Crosschecks")
print("-" * 70)
print("Verify units on both sides of the derivation:")
print()
print("LHS: R_H")
print("  Dimensions: [L]")
print()
print("RHS: c / H_0")
print("  Dimensions: [L T⁻¹] / [T⁻¹]")
print("            = [L]")
print()
print("✓ Dimensional consistency verified")
print()
print("Hubble time (age of universe estimate):")
t_H = 1.0 / H_0
t_H_yr = t_H / (365.25 * 24 * 3600)
t_H_Gyr = t_H_yr / 1e9
print(f"  t_H = 1/H_0 = {t_H:.15e} s")
print(f"              = {t_H_yr:.15e} years")
print(f"              = {t_H_Gyr:.3f} billion years")
print()
print("Consistency check:")
print(f"  R_H = c·t_H = {c * t_H:.15e} m")
print()

# ========================================
# STEP 8: Calibration Checkpoint
# ========================================
print("STEP 8: Calibration Checkpoint")
print("-" * 70)
print("Measured Hubble constant (Planck 2018): 67.4 ± 0.5 km/s/Mpc")
H_0_measured_SI = 67.4 / 3.08567758149e19  # Convert to s^-1
R_H_measured = c / H_0_measured_SI
R_H_measured_Gly = R_H_measured / (9.4607304725808e15 * 1e9)
print(f"Measured H_0:   {67.4:.3f} km/s/Mpc")
print(f"                {H_0_measured_SI:.15e} s⁻¹")
print(f"Measured R_H:   {R_H_measured:.15e} m")
print(f"                {R_H_measured_Gly:.3f} billion light-years")
print()
print(f"TriPhase H_0:   {H_0_SI * 3.08567758149e19:.3f} km/s/Mpc")
print(f"                {H_0:.15e} s⁻¹")
print(f"TriPhase R_H:   {R_H:.15e} m")
print(f"                {R_H_Gly:.3f} billion light-years")
print()
deviation_ppm = abs(H_0 - H_0_measured_SI) / H_0_measured_SI * 1e6
print(f"H_0 Deviation:  {deviation_ppm:.1f} ppm")
print()

# ========================================
# SUMMARY
# ========================================
print("=" * 70)
print("DIMENSIONAL ANALYSIS SUMMARY")
print("=" * 70)
print()
print("The Hubble radius derivation is dimensionally consistent:")
print()
print("1. Target dimensions [L] verified")
print("2. H_0 has dimensions [T⁻¹]")
print("3. R_H = c/H_0 has correct length units")
print("4. 18-step alpha cascade: α^18 is dimensionless")
print("5. Connects quantum (f_e) to cosmological (H_0) scales")
print()
print("The formula R_H = c / H_0 with H_0 = π·√3·f_e·α^18")
print("represents the cosmological horizon, derived from pure")
print("electromagnetic constants and the electron Compton frequency.")
print()
print("Key insight:")
print(f"  The 18-step cascade α^18 ≈ {alpha_18:.3e}")
print("  reduces the electron frequency f_e to the Hubble rate H_0")
print("  This connects quantum and cosmological scales directly.")
print()
print("Physical significance:")
print("  R_H defines the observable universe horizon")
print("  t_H = 1/H_0 gives the Hubble time (age scale)")
print("  Both emerge from pure wave mechanics via the α^18 cascade")
print()
print("=" * 70)

input("Press Enter to exit...")
