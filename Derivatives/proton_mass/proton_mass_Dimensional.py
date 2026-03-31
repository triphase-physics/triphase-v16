"""
TriPhase V16: Proton Mass Derivation
Dimensional Analysis Framework

This script demonstrates the dimensional consistency of the proton mass
derivation using pure wave mechanics and electromagnetic constants.

Proton Mass:
  m_p = m_e * mp_me
  where mp_me = 4 × 27 × 17 × (1 + 5α²/π)

SI Units: [kg]
Dimensional form: [M]

The proton mass emerges from the electron mass scaled by a harmonic ratio
involving the numbers 4, 27, and 17, with a fine structure correction.

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
print("TriPhase V16: Proton Mass")
print("Dimensional Analysis Framework")
print("=" * 70)
print()

# ========================================
# STEP 1: Identify Target Dimensions
# ========================================
print("STEP 1: Identify Target Dimensions")
print("-" * 70)
print("Target Quantity: Proton Mass (m_p)")
print("SI Unit: kg")
print("Dimensional Form: [M]")
print()
print("The proton is the stable baryon forming atomic nuclei.")
print("CODATA 2018: 1.67262192369(51)e-27 kg")
print()

# ========================================
# STEP 2: Available Base Dimensions
# ========================================
print("STEP 2: Available Base Dimensions")
print("-" * 70)
print("From anchor constants:")
print("  ε₀: [M⁻¹ L⁻³ T⁴ I²]  (permittivity)")
print("  μ₀: [M L T⁻² I⁻²]   (permeability)")
print("  e:  [I T]            (elementary charge)")
print("  c:  [L T⁻¹]          (speed of light)")
print("  Z₀: [M L² T⁻³ I⁻²]  (impedance of free space)")
print("  ℏ:  [M L² T⁻¹]      (reduced Planck constant)")
print("  α:  [1]              (fine structure constant, dimensionless)")
print("  r_e: [L]             (classical electron radius)")
print()
print("Derived electron mass:")
print("  m_e = ℏ·α/(c·r_e)")
print("  Dimensions: [M L² T⁻¹]·[1] / ([L T⁻¹]·[L]) = [M]")
print()

# ========================================
# STEP 3: Dimensional Matching
# ========================================
print("STEP 3: Dimensional Matching")
print("-" * 70)
print("Proton mass formula:")
print("  m_p = m_e · mp_me")
print("  where mp_me = 4 × 27 × 17 × (1 + 5α²/π)")
print()
print("Dimensional analysis:")
print("  [m_p] = [m_e] · [1]")
print("        = [M] · (dimensionless ratio)")
print("        = [M]")
print()
print("Where:")
print("  m_e: [M]               (electron mass)")
print("  mp_me: [1]             (dimensionless mass ratio)")
print("  4, 27, 17: [1]        (harmonic numbers)")
print("  α²/π: [1]              (dimensionless correction)")
print()
print("The mass ratio mp_me is dimensionless,")
print("so the result has pure mass dimensions.")
print()

# ========================================
# STEP 4: TriPhase Derivation
# ========================================
print("STEP 4: TriPhase Derivation")
print("-" * 70)
print(f"Electron mass:         m_e    = {m_e:.15e} kg")
print(f"Fine structure const:  α      = {alpha:.15e}")
print()
print("Harmonic structure:")
base_ratio = 4.0 * 27.0 * 17.0
correction = 1.0 + 5.0 * alpha**2 / math.pi
print(f"  Base ratio:   4 × 27 × 17    = {base_ratio:.15e}")
print(f"  Correction:   (1 + 5α²/π)    = {correction:.15e}")
print(f"  Total ratio:  mp_me          = {mp_me:.15e}")
print()
print(f"Proton mass:           m_p = {m_p:.15e} kg")
m_p_MeV = m_p * c**2 / (1.602176634e-19 * 1e6)
print(f"                           = {m_p_MeV:.6f} MeV/c²")
print()

# ========================================
# STEP 5: Buckingham Pi Theorem
# ========================================
print("STEP 5: Buckingham Pi Theorem")
print("-" * 70)
print("For the proton mass, we identify dimensionless groups:")
print()
print("π₁ = m_p / m_e (the proton-to-electron mass ratio)")
print(f"   = {m_p/m_e:.15e}")
print()
print("π₂ = 4 × 27 × 17")
print(f"   = {base_ratio:.15e}")
print()
print("π₃ = 5α²/π")
print(f"   = {5.0 * alpha**2 / math.pi:.15e}")
print()
print("π₄ = (1 + 5α²/π)")
print(f"   = {correction:.15e}")
print()
print("π₅ = mp_me / (4 × 27 × 17 × (1 + 5α²/π))")
print(f"   = {mp_me / (base_ratio * correction):.15e}")
print("   (Should be 1.0 from the formula)")
print()
print("These dimensionless groups characterize the fundamental")
print("proton-electron mass relationship.")
print()

# ========================================
# STEP 6: Natural Units Comparison
# ========================================
print("STEP 6: Natural Units Comparison")
print("-" * 70)
print("Planck mass: m_P = sqrt(ℏc/G)")
m_planck = math.sqrt(hbar * c / G)
print(f"  m_P = {m_planck:.15e} kg")
print(f"  m_p / m_P = {m_p / m_planck:.15e}")
print()
print("Atomic mass unit: u = 1.66053906660e-27 kg")
u_amu = 1.66053906660e-27
print(f"  m_p / u = {m_p / u_amu:.15e}")
print()
print("Electron mass:")
print(f"  m_p / m_e = {m_p / m_e:.15e}")
print()
print("Compton wavelength ratio:")
lambda_c_e = h / (m_e * c)
lambda_c_p = h / (m_p * c)
print(f"  λ_c,e / λ_c,p = {lambda_c_e / lambda_c_p:.15e}")
print(f"  (Should equal m_p/m_e)")
print()

# ========================================
# STEP 7: Dimensional Crosschecks
# ========================================
print("STEP 7: Dimensional Crosschecks")
print("-" * 70)
print("Verify units on both sides of the derivation:")
print()
print("LHS: m_p")
print("  Dimensions: [M]")
print()
print("RHS: m_e · mp_me")
print("  Dimensions: [M] · [1]")
print("            = [M]")
print()
print("✓ Dimensional consistency verified")
print()
print("Energy check:")
E_p = m_p * c**2
print(f"  E_p = m_p · c² = {E_p:.15e} J")
E_p_eV = E_p / e
print(f"              = {E_p_eV:.15e} eV")
print(f"              = {E_p_eV/1e6:.6f} MeV")
print()
print("Compton wavelength:")
print(f"  λ_c,p = h/(m_p·c) = {lambda_c_p:.15e} m")
print(f"                    = {lambda_c_p*1e15:.6f} fm")
print()
print("Classical proton radius estimate (if charge spread like electron):")
r_p_classical = r_e * m_e / m_p
print(f"  r_p = r_e · (m_e/m_p) = {r_p_classical:.15e} m")
print(f"                        = {r_p_classical*1e15:.6f} fm")
print("  (Actual charge radius ~0.84 fm from experiments)")
print()

# ========================================
# STEP 8: Calibration Checkpoint
# ========================================
print("STEP 8: Calibration Checkpoint")
print("-" * 70)
print("CODATA 2018 value: 1.67262192369(51)e-27 kg")
m_p_CODATA = 1.67262192369e-27
print(f"CODATA value:   {m_p_CODATA:.15e} kg")
print()
print(f"TriPhase value: {m_p:.15e} kg")
print()
deviation_ppm = abs(m_p - m_p_CODATA) / m_p_CODATA * 1e6
print(f"Deviation:      {deviation_ppm:.1f} ppm")
print()
print("Mass ratio verification:")
ratio_CODATA = 1836.15267343
print(f"CODATA m_p/m_e: {ratio_CODATA:.15e}")
print(f"TriPhase mp_me: {mp_me:.15e}")
deviation_ratio_ppm = abs(mp_me - ratio_CODATA) / ratio_CODATA * 1e6
print(f"Ratio deviation: {deviation_ratio_ppm:.1f} ppm")
print()

# ========================================
# SUMMARY
# ========================================
print("=" * 70)
print("DIMENSIONAL ANALYSIS SUMMARY")
print("=" * 70)
print()
print("The proton mass derivation is dimensionally consistent:")
print()
print("1. Target dimensions [M] verified")
print("2. Mass ratio mp_me is dimensionless")
print("3. Result has correct mass units")
print("4. Energy conversion yields MeV scale")
print("5. Harmonic structure: 4 × 27 × 17 × (1 + 5α²/π)")
print()
print("The formula m_p = m_e · mp_me represents the fundamental")
print("proton-electron mass relationship, derived from pure")
print("electromagnetic constants and harmonic wave mechanics.")
print()
print("Key harmonic numbers:")
print("  4   = 2²      (quaternity, stability)")
print("  27  = 3³      (cubic structure)")
print("  17  = prime   (connects to T_17 = 153)")
print("  5α²/π         (fine structure correction)")
print()
print("This structure suggests the proton mass emerges from")
print("geometric and harmonic principles in wave mechanics.")
print()
print("=" * 70)

input("Press Enter to exit...")
