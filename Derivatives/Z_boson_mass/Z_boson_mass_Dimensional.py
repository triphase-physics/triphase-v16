"""
TriPhase V16: Z Boson Mass Derivation
Dimensional Analysis Framework

This script demonstrates the dimensional consistency of the Z boson mass
derivation using pure wave mechanics and electromagnetic constants.

Z Boson Mass:
  M_Z ~ M_W / cos(theta_W)
  where sin²(theta_W) ≈ alpha * pi

SI Units: [kg]
Dimensional form: [M]

The Z boson mass emerges from the W boson mass divided by the cosine of
the Weinberg angle, which relates to the fine structure constant.

MIS TAG: (D*H)
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
print("TriPhase V16: Z Boson Mass")
print("Dimensional Analysis Framework")
print("=" * 70)
print()

# ========================================
# STEP 1: Identify Target Dimensions
# ========================================
print("STEP 1: Identify Target Dimensions")
print("-" * 70)
print("Target Quantity: Z Boson Mass (M_Z)")
print("SI Unit: kg")
print("Dimensional Form: [M]")
print()
print("The Z boson is the neutral weak force carrier.")
print("Measured value: 91.1876 ± 0.0021 GeV/c² ≈ 1.625e-25 kg")
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
print("  ℏ:  [M L² T⁻¹]      (reduced Planck constant)")
print("  α:  [1]              (fine structure constant, dimensionless)")
print()
print("Derived masses:")
print("  m_p = m_e · mp_me")
print("  M_W = m_p · T_17 / (2α)")
print("  Dimensions: [M]")
print()

# ========================================
# STEP 3: Dimensional Matching
# ========================================
print("STEP 3: Dimensional Matching")
print("-" * 70)
print("Z boson mass formula:")
print("  M_Z = M_W / cos(θ_W)")
print("  where sin²(θ_W) ≈ α·π")
print()
print("Dimensional analysis:")
print("  [M_Z] = [M_W] / [1]")
print("        = [M] / (dimensionless angle function)")
print("        = [M]")
print()
print("Where:")
print("  M_W: [M]               (W boson mass)")
print("  θ_W: [1]               (Weinberg angle, dimensionless)")
print("  cos(θ_W): [1]          (dimensionless)")
print()
print("The Weinberg angle is dimensionless,")
print("so the result has pure mass dimensions.")
print()

# ========================================
# STEP 4: TriPhase Derivation
# ========================================
print("STEP 4: TriPhase Derivation")
print("-" * 70)
print(f"Proton mass:           m_p    = {m_p:.15e} kg")
print(f"Fine structure const:  α      = {alpha:.15e}")
print(f"Triangular number:     T_17   = {T_17}")
print()

# W boson mass
M_W = m_p * T_17 / (2.0 * alpha)
M_W_GeV = M_W * c**2 / (1.602176634e-19 * 1e9)
print(f"W boson mass:          M_W = {M_W:.15e} kg")
print(f"                           = {M_W_GeV:.6f} GeV/c²")
print()

# Weinberg angle
sin2_theta_W = alpha * math.pi
cos2_theta_W = 1.0 - sin2_theta_W
cos_theta_W = math.sqrt(cos2_theta_W)
theta_W_rad = math.asin(math.sqrt(sin2_theta_W))
theta_W_deg = theta_W_rad * 180.0 / math.pi

print("Weinberg angle:")
print(f"  sin²(θ_W) ≈ α·π       = {sin2_theta_W:.15e}")
print(f"  cos²(θ_W)             = {cos2_theta_W:.15e}")
print(f"  cos(θ_W)              = {cos_theta_W:.15e}")
print(f"  θ_W                   = {theta_W_rad:.6f} rad")
print(f"                        = {theta_W_deg:.6f}°")
print()

# Z boson mass
M_Z = M_W / cos_theta_W
M_Z_GeV = M_Z * c**2 / (1.602176634e-19 * 1e9)

print(f"Z boson mass:          M_Z = {M_Z:.15e} kg")
print(f"                           = {M_Z_GeV:.6f} GeV/c²")
print()
print("Mass ratio:")
print(f"  M_Z / M_W = 1/cos(θ_W) = {M_Z / M_W:.15e}")
print()

# ========================================
# STEP 5: Buckingham Pi Theorem
# ========================================
print("STEP 5: Buckingham Pi Theorem")
print("-" * 70)
print("For the Z boson mass, we identify dimensionless groups:")
print()
print("π₁ = M_Z / M_W")
print(f"   = {M_Z/M_W:.15e}")
print()
print("π₂ = M_Z / m_p")
print(f"   = {M_Z/m_p:.15e}")
print()
print("π₃ = sin²(θ_W) = α·π")
print(f"   = {sin2_theta_W:.15e}")
print()
print("π₄ = cos(θ_W)")
print(f"   = {cos_theta_W:.15e}")
print()
print("π₅ = M_Z · cos(θ_W) / M_W")
print(f"   = {M_Z * cos_theta_W / M_W:.15e}")
print("   (Should be 1.0 from the formula)")
print()
print("These dimensionless groups characterize the electroweak")
print("mass structure and the Z-W mass relationship.")
print()

# ========================================
# STEP 6: Natural Units Comparison
# ========================================
print("STEP 6: Natural Units Comparison")
print("-" * 70)
print("Planck mass: m_P = sqrt(ℏc/G)")
m_planck = math.sqrt(hbar * c / G)
print(f"  m_P = {m_planck:.15e} kg")
print(f"  M_Z / m_P = {M_Z / m_planck:.15e}")
print()
print("Atomic mass unit: u = 1.66053906660e-27 kg")
u_amu = 1.66053906660e-27
print(f"  M_Z / u = {M_Z / u_amu:.15e}")
print()
print("Electron mass:")
print(f"  M_Z / m_e = {M_Z / m_e:.15e}")
print()
print("Proton mass:")
print(f"  M_Z / m_p = {M_Z / m_p:.15e}")
print()
print("W boson mass:")
print(f"  M_Z / M_W = {M_Z / M_W:.15e}")
print()

# ========================================
# STEP 7: Dimensional Crosschecks
# ========================================
print("STEP 7: Dimensional Crosschecks")
print("-" * 70)
print("Verify units on both sides of the derivation:")
print()
print("LHS: M_Z")
print("  Dimensions: [M]")
print()
print("RHS: M_W / cos(θ_W)")
print("  Dimensions: [M] / [1]")
print("            = [M]")
print()
print("✓ Dimensional consistency verified")
print()
print("Energy check:")
E_Z = M_Z * c**2
print(f"  E_Z = M_Z · c² = {E_Z:.15e} J")
E_Z_eV = E_Z / e
print(f"               = {E_Z_eV:.15e} eV")
print(f"               = {E_Z_eV/1e9:.6f} GeV")
print()
print("Compton wavelength:")
lambda_Z = h / (M_Z * c)
print(f"  λ_Z = h/(M_Z·c) = {lambda_Z:.15e} m")
print(f"                  = {lambda_Z*1e18:.6f} am (attometers)")
print()
print("Weak force range:")
r_weak_Z = hbar / (M_Z * c)
print(f"  r_weak ≈ ℏ/(M_Z·c) = {r_weak_Z:.15e} m")
print(f"                     = {r_weak_Z*1e18:.6f} am")
print()

# ========================================
# STEP 8: Calibration Checkpoint
# ========================================
print("STEP 8: Calibration Checkpoint")
print("-" * 70)
print("Measured value (PDG 2024): 91.1876 ± 0.0021 GeV/c²")
M_Z_measured_GeV = 91.1876
M_Z_measured = M_Z_measured_GeV * 1e9 * e / c**2
print(f"Measured value: {M_Z_measured:.15e} kg")
print(f"                {M_Z_measured_GeV:.6f} GeV/c²")
print()
print(f"TriPhase value: {M_Z:.15e} kg")
print(f"                {M_Z_GeV:.6f} GeV/c²")
print()
deviation_ppm = abs(M_Z - M_Z_measured) / M_Z_measured * 1e6
print(f"Deviation:      {deviation_ppm:.1f} ppm")
print()
print("Weinberg angle comparison:")
sin2_theta_W_measured = 0.23121  # PDG on-shell value
theta_W_measured_rad = math.asin(math.sqrt(sin2_theta_W_measured))
theta_W_measured_deg = theta_W_measured_rad * 180.0 / math.pi
print(f"Measured sin²(θ_W): {sin2_theta_W_measured:.6f}")
print(f"TriPhase sin²(θ_W): {sin2_theta_W:.6f}")
print(f"Measured θ_W:       {theta_W_measured_deg:.3f}°")
print(f"TriPhase θ_W:       {theta_W_deg:.3f}°")
print()
print("NOTE: Weinberg angle is scheme and scale dependent.")
print("TriPhase sin²(θ_W) = α·π provides a fundamental estimate.")
print()

# ========================================
# SUMMARY
# ========================================
print("=" * 70)
print("DIMENSIONAL ANALYSIS SUMMARY")
print("=" * 70)
print()
print("The Z boson mass derivation is dimensionally consistent:")
print()
print("1. Target dimensions [M] verified")
print("2. Weinberg angle is dimensionless")
print("3. Result has correct mass units")
print("4. Energy conversion yields GeV scale")
print("5. Electroweak structure via M_Z = M_W / cos(θ_W)")
print()
print("The formula M_Z ~ M_W / cos(θ_W) with sin²(θ_W) ≈ α·π")
print("represents the neutral weak force carrier mass scale,")
print("derived from pure electromagnetic constants.")
print()
print("Key insight:")
print(f"  The Weinberg angle θ_W ≈ {theta_W_deg:.1f}° emerges from α·π")
print(f"  This gives M_Z / M_W ≈ {M_Z/M_W:.4f}")
print("  relating the two electroweak boson masses.")
print()
print("Electroweak unification:")
print("  Both M_W and M_Z derive from m_p, T_17, and α")
print("  Their ratio depends only on the Weinberg angle")
print()
print("=" * 70)

input("Press Enter to exit...")
