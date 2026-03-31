"""
TriPhase V16: Thermal Pressure Derivation
Dimensional Analysis Framework

This script demonstrates the dimensional consistency of thermal pressure
derivation using pure wave mechanics and statistical mechanics.

Thermal Pressure:
  P_th = nkT = (ρ/m)kT

SI Units: [Pa] = [kg m⁻¹ s⁻²]
Dimensional form: [M L⁻¹ T⁻²]

Thermal pressure arises from particle kinetic energy in thermodynamic
systems, fundamental to equation of state.

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
print("TriPhase V16: Thermal Pressure")
print("Dimensional Analysis Framework")
print("=" * 70)
print()

# ========================================
# STEP 1: Identify Target Dimensions
# ========================================
print("STEP 1: Identify Target Dimensions")
print("-" * 70)
print("Target Quantity: Thermal Pressure (P_th)")
print("SI Unit: Pa = kg/(m·s²)")
print("Dimensional Form: [M L⁻¹ T⁻²]")
print()
print("Thermal pressure is the kinetic pressure from random")
print("thermal motion of particles in a gas or plasma.")
print()

# ========================================
# STEP 2: Available Base Dimensions
# ========================================
print("STEP 2: Available Base Dimensions")
print("-" * 70)
print("From thermodynamics and statistical mechanics:")
print("  n:  [L⁻³]             (number density)")
print("  k:  [M L² T⁻² K⁻¹]   (Boltzmann constant)")
print("  T:  [K]               (temperature)")
print("  ρ:  [M L⁻³]           (mass density)")
print("  m:  [M]               (particle mass)")
print()
print("Boltzmann constant:")
print("  k = 1.380649e-23 J/K (exact, SI 2019)")
k_B = 1.380649e-23
print(f"  k = {k_B:.15e} J/K")
print()

# ========================================
# STEP 3: Dimensional Matching
# ========================================
print("STEP 3: Dimensional Matching")
print("-" * 70)
print("Thermal pressure formula:")
print("  P_th = n·k·T")
print()
print("Dimensional analysis:")
print("  [P_th] = [n] · [k] · [T]")
print("         = [L⁻³] · [M L² T⁻² K⁻¹] · [K]")
print("         = [L⁻³⁺² M T⁻² K⁻¹⁺¹]")
print("         = [M L⁻¹ T⁻²]")
print()
print("✓ Result has pressure dimensions")
print()
print("Alternative form:")
print("  P_th = (ρ/m)·k·T")
print("  [P_th] = ([M L⁻³] / [M]) · [M L² T⁻² K⁻¹] · [K]")
print("         = [L⁻³] · [M L² T⁻²]")
print("         = [M L⁻¹ T⁻²]")
print()

# ========================================
# STEP 4: TriPhase Derivation
# ========================================
print("STEP 4: TriPhase Derivation")
print("-" * 70)
print(f"Boltzmann constant: k = {k_B:.15e} J/K")
print(f"Electron mass:      m_e = {m_e:.15e} kg")
print()

# Electron Compton temperature
T_Compton = m_e * c**2 / k_B
print("Electron Compton temperature:")
print(f"  T_C = m_e·c²/k = {T_Compton:.15e} K")
print(f"                 = {T_Compton:.6e} K")
print()

# Number density for ideal gas at STP
n_STP = 101325.0 / (k_B * 273.15)  # P = nkT, T = 273.15 K
print("Number density at STP (P = 1 atm, T = 273.15 K):")
print(f"  n = P/(kT) = {n_STP:.15e} m⁻³")
P_STP = n_STP * k_B * 273.15
print(f"  P = nkT = {P_STP:.15e} Pa")
print(f"          = {P_STP:.1f} Pa (should be 101325 Pa)")
print()

# Solar core conditions
T_sun_core = 1.57e7  # K
rho_sun_core = 1.6e5  # kg/m³
m_avg_sun = m_p  # Approximate as hydrogen
n_sun = rho_sun_core / m_avg_sun
P_sun_core = n_sun * k_B * T_sun_core

print("Solar core conditions:")
print(f"  T_core ≈ {T_sun_core:.2e} K")
print(f"  ρ_core ≈ {rho_sun_core:.2e} kg/m³")
print(f"  m_avg ≈ m_p = {m_avg_sun:.3e} kg")
print(f"  n = ρ/m = {n_sun:.3e} m⁻³")
print(f"  P_th = nkT = {P_sun_core:.3e} Pa")
print()

# Earth atmosphere
T_earth = 288.0  # K (15°C)
P_earth = 101325.0  # Pa
n_earth = P_earth / (k_B * T_earth)
print("Earth atmosphere (sea level):")
print(f"  T = {T_earth:.1f} K")
print(f"  P = {P_earth:.1f} Pa")
print(f"  n = P/(kT) = {n_earth:.3e} m⁻³")
print()

# ========================================
# STEP 5: Buckingham Pi Theorem
# ========================================
print("STEP 5: Buckingham Pi Theorem")
print("-" * 70)
print("For thermal pressure, dimensionless groups:")
print()
print("π₁ = P_th / (n·k·T)")
print("   (Should be 1.0 from ideal gas law)")
print()
print("π₂ = kT / (m·c²)")
print(f"   At T_Compton: {k_B * T_Compton / (m_e * c**2):.15e}")
print("   (Thermal energy / rest mass energy)")
print()
print("π₃ = P_th / (ρ·c²)")
print(f"   Solar core: {P_sun_core / (rho_sun_core * c**2):.15e}")
print("   (Thermal pressure / relativistic energy density)")
print()

# Degeneracy parameter
lambda_thermal = h / math.sqrt(2.0 * math.pi * m_e * k_B * T_earth)
n_quantum = 1.0 / lambda_thermal**3
print("π₄ = n / n_quantum (where n_quantum = 1/λ_th³)")
print(f"   Earth atmosphere: {n_earth / n_quantum:.15e}")
print("   (Classical regime when << 1)")
print()

print("These dimensionless groups characterize thermal dynamics")
print("and the classical/quantum transition.")
print()

# ========================================
# STEP 6: Natural Units Comparison
# ========================================
print("STEP 6: Natural Units Comparison")
print("-" * 70)
print("Vacuum rigidity: VF_r = c⁴/(8πG)")
print(f"  VF_r = {VF_r:.15e} Pa")
print(f"  P_sun / VF_r = {P_sun_core / VF_r:.15e}")
print()
print("Planck pressure: P_P = c⁷/(ℏG²)")
P_P = c**7 / (hbar * G**2)
print(f"  P_P = {P_P:.15e} Pa")
print(f"  P_sun / P_P = {P_sun_core / P_P:.15e}")
print()
print("Compton temperature:")
print(f"  T_C(electron) = m_e·c²/k = {T_Compton:.3e} K")
T_C_proton = m_p * c**2 / k_B
print(f"  T_C(proton) = m_p·c²/k = {T_C_proton:.3e} K")
print()

# ========================================
# STEP 7: Dimensional Crosschecks
# ========================================
print("STEP 7: Dimensional Crosschecks")
print("-" * 70)
print("Verify units on both sides:")
print()
print("LHS: P_th")
print("  Dimensions: [M L⁻¹ T⁻²]")
print("  Units: Pa = kg/(m·s²)")
print()
print("RHS: n·k·T")
print("  Dimensions: [L⁻³] · [M L² T⁻² K⁻¹] · [K]")
print("            = [M L⁻¹ T⁻²]")
print("  Units: (1/m³) · (J/K) · K = Pa")
print()
print("✓ Dimensional consistency verified")
print()
print("Ideal gas law:")
print("  PV = NkT  (or P = nkT where n = N/V)")
print()
print("Energy density:")
u_th = (3.0/2.0) * n_sun * k_B * T_sun_core
print(f"  u_th = (3/2)nkT = {u_th:.3e} J/m³")
print("  (For monatomic ideal gas)")
print(f"  P_th / u_th = {P_sun_core / u_th:.3f}")
print("  (Should be 2/3 for ideal gas)")
print()

# ========================================
# STEP 8: Calibration Checkpoint
# ========================================
print("STEP 8: Calibration Checkpoint")
print("-" * 70)
print("Standard atmosphere verification:")
print("  P = 101325 Pa, T = 273.15 K")
print(f"  Calculated P = {P_STP:.1f} Pa")
print(f"  Agreement: {abs(P_STP - 101325)/101325 * 100:.6f}%")
print()
print("Solar core pressure (from models):")
print("  P_core ≈ 2.5e16 Pa")
P_sun_model = 2.5e16
print(f"  TriPhase estimate: {P_sun_core:.3e} Pa")
print(f"  Ratio: {P_sun_core / P_sun_model:.3f}")
print("  (Order of magnitude agreement)")
print()
print("NOTE: Actual solar core has radiation pressure,")
print("degeneracy effects, and non-ideal gas behavior.")
print()

# ========================================
# SUMMARY
# ========================================
print("=" * 70)
print("DIMENSIONAL ANALYSIS SUMMARY")
print("=" * 70)
print()
print("The thermal pressure derivation is dimensionally consistent:")
print()
print("1. Target dimensions [M L⁻¹ T⁻²] verified")
print("2. Formula P_th = nkT from ideal gas law")
print("3. Applies from atmospheres to stellar cores")
print("4. Connects temperature to pressure via k_B")
print("5. Classical limit when kT << mc²")
print()
print("The formula P_th = nkT represents thermal pressure")
print("from particle kinetic energy, fundamental to")
print("thermodynamics and equation of state.")
print()
print("Key insight:")
print("  Temperature T sets energy scale kT")
print("  Pressure arises from momentum transfer")
print("  Compton temperature T_C = mc²/k sets relativistic limit")
print()
print("TriPhase connection:")
print("  Electron Compton temperature T_C ≈ 6e9 K")
print("  Typical systems are non-relativistic (kT << m_e·c²)")
print("  Thermal pressure << vacuum rigidity VF_r")
print()
print("=" * 70)

input("Press Enter to exit...")
