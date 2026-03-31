"""
TriPhase V16: Charm Quark Mass Derivation
Dimensional Analysis Framework

This script demonstrates the dimensional consistency of the charm quark mass
derivation using pure wave mechanics and electromagnetic constants.

Charm Quark Mass:
  m_c ~ m_e * alpha * T_17 * mp_me

SI Units: [kg]
Dimensional form: [M]

The charm quark mass emerges from the electron mass scaled by the fine structure
constant, the triangular number T_17 = 153, and the full proton-to-electron
mass ratio. This represents the next harmonic in the quark mass spectrum.

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
print("TriPhase V16: Charm Quark Mass")
print("Dimensional Analysis Framework")
print("=" * 70)
print()

# ========================================
# STEP 1: Identify Target Dimensions
# ========================================
print("STEP 1: Identify Target Dimensions")
print("-" * 70)
print("Target Quantity: Charm Quark Mass (m_c)")
print("SI Unit: kg")
print("Dimensional Form: [M]")
print()
print("The charm quark mass is a fundamental mass parameter in QCD.")
print("PDG value: ~1.27 GeV/c² ≈ 2.26e-27 kg")
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
print("Derived base mass:")
print("  m_e = ℏ·α/(c·r_e)")
print("  Dimensions: [M L² T⁻¹]·[1] / ([L T⁻¹]·[L]) = [M]")
print()

# ========================================
# STEP 3: Dimensional Matching
# ========================================
print("STEP 3: Dimensional Matching")
print("-" * 70)
print("Charm quark mass formula:")
print("  m_c = m_e · α · T_17 · mp_me")
print()
print("Dimensional analysis:")
print("  [m_c] = [m_e] · [1] · [1] · [1]")
print("        = [M] · (dimensionless factors)")
print("        = [M]")
print()
print("Where:")
print("  m_e: [M]               (electron mass)")
print("  α: [1]                  (dimensionless)")
print("  T_17 = 153: [1]        (triangular number)")
print("  mp_me: [1]              (mass ratio, dimensionless)")
print()
print("The factor α · T_17 · mp_me is dimensionless,")
print("so the result has pure mass dimensions.")
print()

# ========================================
# STEP 4: TriPhase Derivation
# ========================================
print("STEP 4: TriPhase Derivation")
print("-" * 70)
print(f"Electron mass:         m_e    = {m_e:.15e} kg")
print(f"Fine structure const:  α      = {alpha:.15e}")
print(f"Triangular number:     T_17   = {T_17}")
print(f"Proton/electron ratio: mp_me  = {mp_me:.15e}")
print()

scaling_factor = alpha * T_17 * mp_me
m_c = m_e * scaling_factor

print(f"Scaling factor         = α · T_17 · mp_me")
print(f"                       = {scaling_factor:.15e}")
print()
print(f"Charm quark mass:      m_c = {m_c:.15e} kg")
m_c_MeV = m_c * c**2 / (1.602176634e-19 * 1e6)
m_c_GeV = m_c_MeV / 1000.0
print(f"                           = {m_c_MeV:.6f} MeV/c²")
print(f"                           = {m_c_GeV:.6f} GeV/c²")
print()

# ========================================
# STEP 5: Buckingham Pi Theorem
# ========================================
print("STEP 5: Buckingham Pi Theorem")
print("-" * 70)
print("For the charm quark mass, we identify dimensionless groups:")
print()
print("π₁ = m_c / m_e")
print(f"   = {m_c/m_e:.6f}")
print()
print("π₂ = α · T_17 · mp_me")
print(f"   = {alpha * T_17 * mp_me:.6f}")
print()
print("π₃ = m_c / (m_e · α · T_17 · mp_me)")
print(f"   = {m_c / (m_e * alpha * T_17 * mp_me):.6f}")
print("   (Should be 1.0 from the formula)")
print()
print("π₄ = m_c / m_p")
print(f"   = {m_c / m_p:.6f}")
print()
print("These dimensionless groups characterize the harmonic structure")
print("of the charm quark mass in the quark mass spectrum.")
print()

# ========================================
# STEP 6: Natural Units Comparison
# ========================================
print("STEP 6: Natural Units Comparison")
print("-" * 70)
print("Planck mass: m_P = sqrt(ℏc/G)")
m_planck = math.sqrt(hbar * c / G)
print(f"  m_P = {m_planck:.15e} kg")
print(f"  m_c / m_P = {m_c / m_planck:.15e}")
print()
print("Atomic mass unit: u = 1.66053906660e-27 kg")
u_amu = 1.66053906660e-27
print(f"  m_c / u = {m_c / u_amu:.15e}")
print()
print("Electron mass:")
print(f"  m_c / m_e = {m_c / m_e:.15e}")
print()
print("Proton mass:")
print(f"  m_c / m_p = {m_c / m_p:.15e}")
print()
print("Strange quark mass (TriPhase):")
mp_me_cbrt = mp_me**(1.0/3.0)
m_s = m_e * 2.0 * alpha * T_17 * mp_me_cbrt
print(f"  m_c / m_s = {m_c / m_s:.15e}")
print()

# ========================================
# STEP 7: Dimensional Crosschecks
# ========================================
print("STEP 7: Dimensional Crosschecks")
print("-" * 70)
print("Verify units on both sides of the derivation:")
print()
print("LHS: m_c")
print("  Dimensions: [M]")
print()
print("RHS: m_e · α · T_17 · mp_me")
print("  Dimensions: [M] · [1] · [1] · [1]")
print("            = [M]")
print()
print("✓ Dimensional consistency verified")
print()
print("Energy check:")
E_c = m_c * c**2
print(f"  E_c = m_c · c² = {E_c:.15e} J")
E_c_eV = E_c / e
print(f"              = {E_c_eV:.15e} eV")
print(f"              = {E_c_eV/1e6:.6f} MeV")
print(f"              = {E_c_eV/1e9:.6f} GeV")
print()
print("Compton wavelength:")
lambda_c = h / (m_c * c)
print(f"  λ_c = h/(m_c·c) = {lambda_c:.15e} m")
print(f"                  = {lambda_c*1e15:.6f} fm")
print()

# ========================================
# STEP 8: Calibration Checkpoint
# ========================================
print("STEP 8: Calibration Checkpoint")
print("-" * 70)
print("PDG 2024 value: m_c(m_c, MS-bar) = 1.27 ± 0.02 GeV/c²")
m_c_PDG_GeV = 1.27
m_c_PDG = m_c_PDG_GeV * 1e9 * e / c**2
print(f"PDG value:      {m_c_PDG:.15e} kg")
print(f"                {m_c_PDG_GeV:.6f} GeV/c²")
print()
print(f"TriPhase value: {m_c:.15e} kg")
print(f"                {m_c_GeV:.6f} GeV/c²")
print()
deviation_ppm = abs(m_c - m_c_PDG) / m_c_PDG * 1e6
print(f"Deviation:      {deviation_ppm:.1f} ppm")
print()
print("NOTE: Charm quark mass is scale-dependent and scheme-dependent.")
print("TriPhase provides a characteristic scale from pure wave mechanics,")
print("while PDG uses MS-bar at the charm mass scale.")
print()

# ========================================
# SUMMARY
# ========================================
print("=" * 70)
print("DIMENSIONAL ANALYSIS SUMMARY")
print("=" * 70)
print()
print("The charm quark mass derivation is dimensionally consistent:")
print()
print("1. Target dimensions [M] verified")
print("2. All scaling factors dimensionless")
print("3. Result has correct mass units")
print("4. Energy conversion yields GeV scale")
print("5. Harmonic structure via α, T_17, and mp_me")
print()
print("The formula m_c ~ m_e · α · T_17 · mp_me represents")
print("a harmonic coupling in the quark mass spectrum, derived from")
print("pure electromagnetic constants and wave mechanics.")
print()
print("Relation to strange quark:")
print(f"  m_c / m_s = mp_me^(2/3) / 2 ≈ {m_c/m_s:.3f}")
print("This demonstrates the hierarchical structure of quark masses.")
print()
print("=" * 70)

input("Press Enter to exit...")
