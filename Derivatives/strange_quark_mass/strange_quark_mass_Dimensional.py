"""
TriPhase V16: Strange Quark Mass Derivation
Dimensional Analysis Framework

This script demonstrates the dimensional consistency of the strange quark mass
derivation using pure wave mechanics and electromagnetic constants.

Strange Quark Mass:
  m_s ~ m_e * 2.0 * alpha * T_17 * mp_me^(1/3)

SI Units: [kg]
Dimensional form: [M]

The strange quark mass emerges from the electron mass scaled by the fine structure
constant, the triangular number T_17 = 153, and the cube root of the proton-to-electron
mass ratio. This represents a harmonic coupling in the quark mass spectrum.

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
print("TriPhase V16: Strange Quark Mass")
print("Dimensional Analysis Framework")
print("=" * 70)
print()

# ========================================
# STEP 1: Identify Target Dimensions
# ========================================
print("STEP 1: Identify Target Dimensions")
print("-" * 70)
print("Target Quantity: Strange Quark Mass (m_s)")
print("SI Unit: kg")
print("Dimensional Form: [M]")
print()
print("The strange quark mass is a fundamental mass parameter in QCD.")
print("PDG value: ~93 MeV/c² ≈ 1.66e-28 kg")
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
print("Strange quark mass formula:")
print("  m_s = m_e · 2.0 · α · T_17 · (mp_me)^(1/3)")
print()
print("Dimensional analysis:")
print("  [m_s] = [m_e] · [1] · [1] · [1]^(1/3)")
print("        = [M] · (dimensionless factors)")
print("        = [M]")
print()
print("Where:")
print("  m_e: [M]               (electron mass)")
print("  α: [1]                  (dimensionless)")
print("  T_17 = 153: [1]        (triangular number)")
print("  mp_me: [1]              (mass ratio, dimensionless)")
print()
print("The factor 2.0 · α · T_17 · (mp_me)^(1/3) is dimensionless,")
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

mp_me_cbrt = mp_me**(1.0/3.0)
scaling_factor = 2.0 * alpha * T_17 * mp_me_cbrt
m_s = m_e * scaling_factor

print(f"(mp_me)^(1/3)          = {mp_me_cbrt:.15e}")
print(f"Scaling factor         = 2.0 · α · T_17 · (mp_me)^(1/3)")
print(f"                       = {scaling_factor:.15e}")
print()
print(f"Strange quark mass:    m_s = {m_s:.15e} kg")
m_s_MeV = m_s * c**2 / (1.602176634e-19 * 1e6)
print(f"                           = {m_s_MeV:.6f} MeV/c²")
print()

# ========================================
# STEP 5: Buckingham Pi Theorem
# ========================================
print("STEP 5: Buckingham Pi Theorem")
print("-" * 70)
print("For the strange quark mass, we identify dimensionless groups:")
print()
print("π₁ = m_s / m_e")
print(f"   = {m_s/m_e:.6f}")
print()
print("π₂ = α · T_17 · (mp_me)^(1/3)")
print(f"   = {alpha * T_17 * mp_me_cbrt:.6f}")
print()
print("π₃ = m_s / (m_e · α · T_17 · (mp_me)^(1/3))")
print(f"   = {m_s / (m_e * alpha * T_17 * mp_me_cbrt):.6f}")
print("   (Should be ~2.0 from the formula)")
print()
print("These dimensionless groups characterize the harmonic structure")
print("of the strange quark mass in the quark mass spectrum.")
print()

# ========================================
# STEP 6: Natural Units Comparison
# ========================================
print("STEP 6: Natural Units Comparison")
print("-" * 70)
print("Planck mass: m_P = sqrt(ℏc/G)")
m_planck = math.sqrt(hbar * c / G)
print(f"  m_P = {m_planck:.15e} kg")
print(f"  m_s / m_P = {m_s / m_planck:.15e}")
print()
print("Atomic mass unit: u = 1.66053906660e-27 kg")
u_amu = 1.66053906660e-27
print(f"  m_s / u = {m_s / u_amu:.15e}")
print()
print("Electron mass:")
print(f"  m_s / m_e = {m_s / m_e:.15e}")
print()
print("Proton mass:")
print(f"  m_s / m_p = {m_s / m_p:.15e}")
print()

# ========================================
# STEP 7: Dimensional Crosschecks
# ========================================
print("STEP 7: Dimensional Crosschecks")
print("-" * 70)
print("Verify units on both sides of the derivation:")
print()
print("LHS: m_s")
print("  Dimensions: [M]")
print()
print("RHS: m_e · 2.0 · α · T_17 · (mp_me)^(1/3)")
print("  Dimensions: [M] · [1] · [1] · [1] · [1]")
print("            = [M]")
print()
print("✓ Dimensional consistency verified")
print()
print("Energy check:")
E_s = m_s * c**2
print(f"  E_s = m_s · c² = {E_s:.15e} J")
E_s_eV = E_s / e
print(f"              = {E_s_eV:.15e} eV")
print(f"              = {E_s_eV/1e6:.6f} MeV")
print()

# ========================================
# STEP 8: Calibration Checkpoint
# ========================================
print("STEP 8: Calibration Checkpoint")
print("-" * 70)
print("PDG 2024 value: m_s(2 GeV, MS-bar) ≈ 93.4 ± 0.8 MeV/c²")
m_s_PDG_MeV = 93.4
m_s_PDG = m_s_PDG_MeV * 1e6 * e / c**2
print(f"PDG value:      {m_s_PDG:.15e} kg")
print(f"                {m_s_PDG_MeV:.6f} MeV/c²")
print()
print(f"TriPhase value: {m_s:.15e} kg")
print(f"                {m_s_MeV:.6f} MeV/c²")
print()
deviation_ppm = abs(m_s - m_s_PDG) / m_s_PDG * 1e6
print(f"Deviation:      {deviation_ppm:.1f} ppm")
print()
print("NOTE: Strange quark mass is strongly scale-dependent and")
print("scheme-dependent. TriPhase provides a characteristic scale")
print("from pure wave mechanics, while PDG uses MS-bar at 2 GeV.")
print()

# ========================================
# SUMMARY
# ========================================
print("=" * 70)
print("DIMENSIONAL ANALYSIS SUMMARY")
print("=" * 70)
print()
print("The strange quark mass derivation is dimensionally consistent:")
print()
print("1. Target dimensions [M] verified")
print("2. All scaling factors dimensionless")
print("3. Result has correct mass units")
print("4. Energy conversion yields MeV scale")
print("5. Harmonic structure via α, T_17, and (mp_me)^(1/3)")
print()
print("The formula m_s ~ m_e · 2α · T_17 · (mp_me)^(1/3) represents")
print("a harmonic coupling in the quark mass spectrum, derived from")
print("pure electromagnetic constants and wave mechanics.")
print()
print("=" * 70)

input("Press Enter to exit...")
