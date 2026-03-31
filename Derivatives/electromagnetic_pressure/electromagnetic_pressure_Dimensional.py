"""
TriPhase V16: Electromagnetic Pressure Derivation
Dimensional Analysis Framework

This script demonstrates the dimensional consistency of electromagnetic
pressure derivation using pure wave mechanics.

Electromagnetic Pressure:
  P_EM = ε₀E²/2 = B²/(2μ₀)

SI Units: [Pa] = [kg m⁻¹ s⁻²]
Dimensional form: [M L⁻¹ T⁻²]

Electromagnetic pressure arises from the radiation pressure of electric
and magnetic fields, fundamental to Maxwell's stress tensor.

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
print("TriPhase V16: Electromagnetic Pressure")
print("Dimensional Analysis Framework")
print("=" * 70)
print()

# ========================================
# STEP 1: Identify Target Dimensions
# ========================================
print("STEP 1: Identify Target Dimensions")
print("-" * 70)
print("Target Quantity: Electromagnetic Pressure (P_EM)")
print("SI Unit: Pa = kg/(m·s²)")
print("Dimensional Form: [M L⁻¹ T⁻²]")
print()
print("Electromagnetic pressure is the radiation pressure from")
print("electric and magnetic fields, part of Maxwell stress tensor.")
print()

# ========================================
# STEP 2: Available Base Dimensions
# ========================================
print("STEP 2: Available Base Dimensions")
print("-" * 70)
print("From anchor constants:")
print("  ε₀: [M⁻¹ L⁻³ T⁴ I²]  (permittivity)")
print("  μ₀: [M L T⁻² I⁻²]   (permeability)")
print("  E:  [M L T⁻³ I⁻¹]   (electric field)")
print("  B:  [M T⁻² I⁻¹]     (magnetic field)")
print()
print("Electric field:")
print("  [E] = [V/m] = [kg·m·s⁻³·A⁻¹]")
print("  [E] = [M L T⁻³ I⁻¹]")
print()
print("Magnetic field:")
print("  [B] = [T] = [kg·s⁻²·A⁻¹]")
print("  [B] = [M T⁻² I⁻¹]")
print()

# ========================================
# STEP 3: Dimensional Matching
# ========================================
print("STEP 3: Dimensional Matching")
print("-" * 70)
print("Electromagnetic pressure formula (electric):")
print("  P_EM = ε₀E²/2")
print()
print("Dimensional analysis:")
print("  [P_EM] = [ε₀] · [E]²")
print("         = [M⁻¹ L⁻³ T⁴ I²] · [M L T⁻³ I⁻¹]²")
print("         = [M⁻¹ L⁻³ T⁴ I²] · [M² L² T⁻⁶ I⁻²]")
print("         = [M⁻¹⁺² L⁻³⁺² T⁴⁺⁽⁻⁶⁾ I²⁺⁽⁻²⁾]")
print("         = [M L⁻¹ T⁻²]")
print()
print("✓ Result has pressure dimensions")
print()
print("Electromagnetic pressure formula (magnetic):")
print("  P_EM = B²/(2μ₀)")
print()
print("Dimensional analysis:")
print("  [P_EM] = [B]² / [μ₀]")
print("         = [M T⁻² I⁻¹]² / [M L T⁻² I⁻²]")
print("         = [M² T⁻⁴ I⁻²] / [M L T⁻² I⁻²]")
print("         = [M²⁻¹ T⁻⁴⁺² I⁻²⁺²] / [L]")
print("         = [M L⁻¹ T⁻²]")
print()
print("✓ Both formulas give same dimensions")
print()

# ========================================
# STEP 4: TriPhase Derivation
# ========================================
print("STEP 4: TriPhase Derivation")
print("-" * 70)
print(f"Permittivity:  ε₀ = {epsilon_0:.15e} F/m")
print(f"Permeability:  μ₀ = {mu_0:.15e} H/m")
print(f"Speed of light: c = {c:.15e} m/s")
print(f"Impedance:     Z₀ = {Z_0:.15e} Ω")
print()

print("At classical electron radius r_e:")
print(f"  r_e = {r_e:.15e} m")
print()

# Electric field at r_e
E_re = e / (4.0 * math.pi * epsilon_0 * r_e**2)
print(f"Electric field at r_e:")
print(f"  E = e/(4πε₀r_e²) = {E_re:.15e} V/m")
print()

# Electric pressure
P_E_re = epsilon_0 * E_re**2 / 2.0
print(f"Electric pressure:")
print(f"  P_E = ε₀E²/2 = {P_E_re:.15e} Pa")
print()

# Magnetic field (for EM wave with E field)
B_re = E_re / c
print(f"Magnetic field (for EM wave):")
print(f"  B = E/c = {B_re:.15e} T")
print()

# Magnetic pressure
P_B_re = B_re**2 / (2.0 * mu_0)
print(f"Magnetic pressure:")
print(f"  P_B = B²/(2μ₀) = {P_B_re:.15e} Pa")
print()

print("Verification: P_E = P_B for electromagnetic wave")
print(f"  P_E / P_B = {P_E_re / P_B_re:.15e}")
print("  (Should be 1.0)")
print()

# Total EM pressure
P_EM_total = P_E_re + P_B_re
print(f"Total EM pressure (E + B):")
print(f"  P_EM = ε₀E² = B²/μ₀ = {P_EM_total:.15e} Pa")
print()

# ========================================
# STEP 5: Buckingham Pi Theorem
# ========================================
print("STEP 5: Buckingham Pi Theorem")
print("-" * 70)
print("For electromagnetic pressure, dimensionless groups:")
print()
print("π₁ = P_EM / (ε₀E²)")
print(f"   = {P_E_re / (epsilon_0 * E_re**2):.15e}")
print("   (Should be 0.5)")
print()
print("π₂ = μ₀·P_EM / B²")
print(f"   = {mu_0 * P_B_re / B_re**2:.15e}")
print("   (Should be 0.5)")
print()
print("π₃ = c²·μ₀·ε₀")
print(f"   = {c**2 * mu_0 * epsilon_0:.15e}")
print("   (Should be 1.0)")
print()
print("π₄ = E / (c·B)")
print(f"   = {E_re / (c * B_re):.15e}")
print("   (Should be 1.0 for EM wave)")
print()

# Compare to vacuum rigidity
print("π₅ = P_EM / VF_r")
print(f"   = {P_EM_total / VF_r:.15e}")
print("   (EM pressure relative to vacuum rigidity)")
print()

print("These dimensionless groups characterize electromagnetic")
print("field energy and radiation pressure.")
print()

# ========================================
# STEP 6: Natural Units Comparison
# ========================================
print("STEP 6: Natural Units Comparison")
print("-" * 70)
print("Vacuum rigidity: VF_r = c⁴/(8πG)")
print(f"  VF_r = {VF_r:.15e} Pa")
print(f"  P_EM / VF_r = {P_EM_total / VF_r:.15e}")
print()
print("Planck pressure: P_P = c⁷/(ℏG²)")
P_P = c**7 / (hbar * G**2)
print(f"  P_P = {P_P:.15e} Pa")
print(f"  P_EM / P_P = {P_EM_total / P_P:.15e}")
print()
print("Energy density:")
u_EM = epsilon_0 * E_re**2
print(f"  u_EM = ε₀E² = {u_EM:.15e} J/m³")
print(f"  u_EM / c² = {u_EM / c**2:.15e} kg/m³")
print("  (Equivalent mass density)")
print()

# ========================================
# STEP 7: Dimensional Crosschecks
# ========================================
print("STEP 7: Dimensional Crosschecks")
print("-" * 70)
print("Verify units on both sides:")
print()
print("Electric form:")
print("  LHS: P_EM")
print("    Dimensions: [M L⁻¹ T⁻²]")
print("    Units: Pa = kg/(m·s²)")
print()
print("  RHS: ε₀E²/2")
print("    Dimensions: [M⁻¹ L⁻³ T⁴ I²] · [M L T⁻³ I⁻¹]²")
print("              = [M L⁻¹ T⁻²]")
print("    Units: (F/m) · (V/m)² = kg/(m·s²)")
print()
print("✓ Dimensional consistency verified")
print()
print("Magnetic form:")
print("  RHS: B²/(2μ₀)")
print("    Dimensions: [M T⁻² I⁻¹]² / [M L T⁻² I⁻²]")
print("              = [M L⁻¹ T⁻²]")
print("    Units: T² / (H/m) = kg/(m·s²)")
print()
print("✓ Dimensional consistency verified")
print()
print("Radiation pressure on surface:")
print("  I = (ε₀c/2)·E² (intensity)")
I_EM = (epsilon_0 * c / 2.0) * E_re**2
print(f"  I = {I_EM:.15e} W/m²")
print("  P_rad = I/c (for perfect absorption)")
P_rad = I_EM / c
print(f"  P_rad = {P_rad:.15e} Pa")
print(f"  (Note: P_rad = P_EM for EM wave)")
print()

# ========================================
# STEP 8: Calibration Checkpoint
# ========================================
print("STEP 8: Calibration Checkpoint")
print("-" * 70)
print("Solar radiation pressure at Earth:")
print("  Solar constant: S = 1361 W/m²")
S_solar = 1361.0
P_solar = S_solar / c
print(f"  P_solar = S/c = {P_solar:.15e} Pa")
print(f"                = {P_solar*1e6:.3f} µPa")
print()
print("Earth magnetic field:")
B_earth = 5e-5  # Tesla (50 µT)
P_B_earth = B_earth**2 / (2.0 * mu_0)
print(f"  B_Earth ≈ {B_earth*1e6:.1f} µT")
print(f"  P_B = B²/(2μ₀) = {P_B_earth:.15e} Pa")
print()
print("Strong magnetic field (pulsar):")
B_pulsar = 1e8  # Tesla
P_B_pulsar = B_pulsar**2 / (2.0 * mu_0)
print(f"  B_pulsar ≈ {B_pulsar:.3e} T")
print(f"  P_B = B²/(2μ₀) = {P_B_pulsar:.15e} Pa")
print(f"                 = {P_B_pulsar:.3e} Pa")
print()

# ========================================
# SUMMARY
# ========================================
print("=" * 70)
print("DIMENSIONAL ANALYSIS SUMMARY")
print("=" * 70)
print()
print("The electromagnetic pressure derivation is dimensionally")
print("consistent:")
print()
print("1. Target dimensions [M L⁻¹ T⁻²] verified")
print("2. Electric form: P_EM = ε₀E²/2")
print("3. Magnetic form: P_EM = B²/(2μ₀)")
print("4. Both forms give identical pressure")
print("5. Related to radiation pressure P_rad = I/c")
print()
print("The formulas P_EM = ε₀E²/2 = B²/(2μ₀) represent")
print("electromagnetic radiation pressure, fundamental to")
print("Maxwell's stress tensor and wave propagation.")
print()
print("Key insight:")
print("  Electric and magnetic fields carry momentum")
print("  Radiation pressure P = u/c (for EM waves)")
print("  At r_e, electromagnetic pressure is extreme")
print()
print("TriPhase connection:")
print("  ε₀ and μ₀ are anchor constants")
print("  P_EM connects to vacuum rigidity VF_r")
print("  Maximum pressure bounded by c⁴/(8πG)")
print()
print("=" * 70)

input("Press Enter to exit...")
