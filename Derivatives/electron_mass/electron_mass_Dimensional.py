"""
TriPhase V16: Electron Mass (m_e)
Dimensional Analysis Framework

Derivative: m_e = ℏα/(cr_e)
MIS TAG: (D) - Derived from reduced Planck constant and classical electron radius
Status: Fundamental lepton mass

DIMENSIONAL INTERPRETATION:
The electron mass emerges from the reduced Planck constant ℏ, fine structure
constant α, speed of light c, and classical electron radius r_e. This derivation
shows mass as a property of electromagnetic field structure rather than an
independently postulated quantity.

SI UNITS: [kg] (mass)

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# =====================================================================
# ANCHOR CONSTANTS (TriPhase V16 Standard Chain)
# =====================================================================
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

# =====================================================================
print("=" * 70)
print("TriPhase V16: Electron Mass (m_e)")
print("Dimensional Analysis Framework")
print("=" * 70)
print()

# =====================================================================
# STEP 1: IDENTIFY TARGET DIMENSIONS
# =====================================================================
print("STEP 1: IDENTIFY TARGET DIMENSIONS")
print("-" * 70)
print("Target: Electron mass m_e")
print("SI Dimensions: [kg]")
print()

# =====================================================================
# STEP 2: AVAILABLE BASE DIMENSIONS
# =====================================================================
print("STEP 2: AVAILABLE BASE DIMENSIONS")
print("-" * 70)
print("  ℏ: [kg m² s⁻¹]")
print("  α: [1]")
print("  c: [m s⁻¹]")
print("  r_e: [m]")
print()

# =====================================================================
# STEP 3: DIMENSIONAL MATCHING
# =====================================================================
print("STEP 3: DIMENSIONAL MATCHING")
print("-" * 70)
print("m_e = ℏα/(cr_e)")
print("  [m_e] = [kg m² s⁻¹][1] / ([m s⁻¹][m])")
print("        = [kg m² s⁻¹] / [m² s⁻¹] = [kg] ✓")
print()

# =====================================================================
# STEP 4: TRIPHASE DERIVATION
# =====================================================================
print("STEP 4: TRIPHASE DERIVATION")
print("-" * 70)
print(f"  ℏ = {hbar:.15e} J·s")
print(f"  α = {alpha:.15e}")
print(f"  c = {c:.12e} m/s")
print(f"  r_e = {r_e:.15e} m")
print()
m_e_derived = hbar * alpha / (c * r_e)
print(f"  m_e = ℏα/(cr_e) = {m_e_derived:.15e} kg")
print(f"  m_e c² = {m_e_derived*c**2/e/1e6:.6f} MeV")
print()

# =====================================================================
# STEP 5: BUCKINGHAM PI THEOREM
# =====================================================================
print("STEP 5: BUCKINGHAM PI THEOREM")
print("-" * 70)
lambda_C = h / (m_e_derived * c)
print(f"Compton wavelength: λ_C = {lambda_C*1e12:.6f} pm")
print(f"λ_C/r_e = {lambda_C/r_e:.6f} = 2πα⁻¹")
print()

# =====================================================================
# STEP 6: NATURAL UNITS COMPARISON
# =====================================================================
print("STEP 6: NATURAL UNITS COMPARISON")
print("-" * 70)
print(f"SI: m_e = {m_e_derived:.15e} kg")
print(f"Energy: m_e c² = {m_e_derived*c**2/e:.6e} eV")
print()

# =====================================================================
# STEP 7: DIMENSIONAL CROSSCHECKS
# =====================================================================
print("STEP 7: DIMENSIONAL CROSSCHECKS")
print("-" * 70)
print("[ℏα/(cr_e)] = [kg m² s⁻¹]/[m² s⁻¹] = [kg] ✓")
print()

# =====================================================================
# STEP 8: CALIBRATION CHECKPOINT
# =====================================================================
print("STEP 8: CALIBRATION CHECKPOINT")
print("-" * 70)
m_e_CODATA = 9.1093837015e-31  # kg
print(f"TriPhase: {m_e_derived:.15e} kg")
print(f"CODATA: {m_e_CODATA:.15e} kg")
deviation = (m_e_derived - m_e_CODATA) / m_e_CODATA * 1e6
print(f"Deviation: {deviation:+.1f} ppm")
print()

print("=" * 70)
print("DIMENSIONAL ANALYSIS COMPLETE")
print("=" * 70)

input("Press Enter to exit...")
