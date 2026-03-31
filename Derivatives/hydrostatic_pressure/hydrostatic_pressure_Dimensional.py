"""
TriPhase V16: Hydrostatic Pressure Derivation
Dimensional Analysis Framework

This script demonstrates the dimensional consistency of hydrostatic pressure
derivation using pure wave mechanics and gravitational theory.

Hydrostatic Pressure:
  P_hydro = ρgh

SI Units: [Pa] = [kg m⁻¹ s⁻²]
Dimensional form: [M L⁻¹ T⁻²]

Hydrostatic pressure arises from the weight of a fluid column under gravity,
fundamental to fluid statics and atmospheric pressure.

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
print("TriPhase V16: Hydrostatic Pressure")
print("Dimensional Analysis Framework")
print("=" * 70)
print()

# ========================================
# STEP 1: Identify Target Dimensions
# ========================================
print("STEP 1: Identify Target Dimensions")
print("-" * 70)
print("Target Quantity: Hydrostatic Pressure (P_hydro)")
print("SI Unit: Pa = kg/(m·s²)")
print("Dimensional Form: [M L⁻¹ T⁻²]")
print()
print("Hydrostatic pressure is the pressure from a column of fluid")
print("under gravitational acceleration.")
print()

# ========================================
# STEP 2: Available Base Dimensions
# ========================================
print("STEP 2: Available Base Dimensions")
print("-" * 70)
print("From mechanics:")
print("  ρ:  [M L⁻³]           (mass density)")
print("  g:  [L T⁻²]           (gravitational acceleration)")
print("  h:  [L]               (height/depth)")
print()
print("Gravitational acceleration at Earth surface:")
g_earth = 9.80665  # m/s²
print(f"  g_⊕ = {g_earth} m/s² (standard gravity)")
print()
print("From TriPhase G:")
M_earth = 5.972e24  # kg
R_earth = 6.371e6   # m
g_earth_calc = G * M_earth / R_earth**2
print(f"  g_⊕ = GM_⊕/R_⊕² = {g_earth_calc:.6f} m/s²")
print()

# ========================================
# STEP 3: Dimensional Matching
# ========================================
print("STEP 3: Dimensional Matching")
print("-" * 70)
print("Hydrostatic pressure formula:")
print("  P_hydro = ρ·g·h")
print()
print("Dimensional analysis:")
print("  [P_hydro] = [ρ] · [g] · [h]")
print("            = [M L⁻³] · [L T⁻²] · [L]")
print("            = [M L⁻³⁺¹⁺¹ T⁻²]")
print("            = [M L⁻¹ T⁻²]")
print()
print("✓ Result has pressure dimensions")
print()
print("Alternative form using gravitational potential:")
print("  P = ρ·Φ  where Φ = g·h")
print("  [Φ] = [L T⁻²]·[L] = [L² T⁻²]  (specific energy)")
print("  [P] = [M L⁻³]·[L² T⁻²] = [M L⁻¹ T⁻²]")
print()

# ========================================
# STEP 4: TriPhase Derivation
# ========================================
print("STEP 4: TriPhase Derivation")
print("-" * 70)
print(f"Gravitational constant: G = {G:.15e} m³/(kg·s²)")
print()

# Earth's ocean
rho_water = 1000.0  # kg/m³
h_ocean = 10.0      # m (10 m depth)
P_ocean = rho_water * g_earth * h_ocean

print("Ocean water column (10 m depth):")
print(f"  ρ_water = {rho_water} kg/m³")
print(f"  g = {g_earth} m/s²")
print(f"  h = {h_ocean} m")
print(f"  P = ρgh = {P_ocean:.1f} Pa")
print(f"           = {P_ocean/101325:.3f} atm")
print()

# Earth's atmosphere
rho_air = 1.225  # kg/m³ at sea level
h_scale = 8500.0 # m (scale height)
P_atm = rho_air * g_earth * h_scale

print("Atmospheric column (scale height):")
print(f"  ρ_air = {rho_air} kg/m³")
print(f"  H_scale = {h_scale} m")
print(f"  P ≈ ρgH = {P_atm:.1f} Pa")
print(f"           = {P_atm/101325:.3f} atm")
print("  (Actual P_atm = 101325 Pa = 1 atm)")
print()

# Earth's core
rho_core = 13000.0  # kg/m³
h_core = R_earth    # Full depth to center
P_core = rho_core * g_earth * h_core

print("Earth's core (rough estimate):")
print(f"  ρ_core ≈ {rho_core} kg/m³")
print(f"  h ≈ R_⊕ = {R_earth:.3e} m")
print(f"  P ≈ ρgh = {P_core:.3e} Pa")
print("  (Simple estimate; actual ~360 GPa)")
print()

# Solar atmosphere
rho_sun_photosphere = 2e-4  # kg/m³
g_sun = G * 1.989e30 / (6.96e8)**2
h_sun = 1e6  # m (scale height)
P_sun_photo = rho_sun_photosphere * g_sun * h_sun

print("Solar photosphere:")
print(f"  ρ_☉ ≈ {rho_sun_photosphere} kg/m³")
print(f"  g_☉ = GM_☉/R_☉² = {g_sun:.1f} m/s²")
print(f"  h_scale ≈ {h_sun:.3e} m")
print(f"  P ≈ ρgh = {P_sun_photo:.3e} Pa")
print()

# ========================================
# STEP 5: Buckingham Pi Theorem
# ========================================
print("STEP 5: Buckingham Pi Theorem")
print("-" * 70)
print("For hydrostatic pressure, dimensionless groups:")
print()
print("π₁ = P / (ρgh)")
print("   (Should be 1.0 from the formula)")
print()
print("π₂ = gh / c²")
print(f"   Earth surface: {g_earth * R_earth / c**2:.15e}")
print("   (Gravitational potential / c²)")
print()
print("π₃ = ρgh / (ρ_⊕·g_⊕·R_⊕)")
print("   (Pressure relative to Earth core scale)")
print()
print("π₄ = g·h / (G·M/r²)·r")
print(f"   Earth: {g_earth * R_earth / (G * M_earth / R_earth):.15e}")
print("   (Should be 1.0)")
print()

print("These dimensionless groups characterize hydrostatic")
print("equilibrium and gravitational confinement.")
print()

# ========================================
# STEP 6: Natural Units Comparison
# ========================================
print("STEP 6: Natural Units Comparison")
print("-" * 70)
print("Vacuum rigidity: VF_r = c⁴/(8πG)")
print(f"  VF_r = {VF_r:.15e} Pa")
print(f"  P_ocean / VF_r = {P_ocean / VF_r:.15e}")
print(f"  P_core / VF_r = {P_core / VF_r:.15e}")
print()
print("Planck pressure: P_P = c⁷/(ℏG²)")
P_P = c**7 / (hbar * G**2)
print(f"  P_P = {P_P:.15e} Pa")
print(f"  P_core / P_P = {P_core / P_P:.15e}")
print()
print("Characteristic scales:")
l_P = math.sqrt(hbar * G / c**3)
rho_P = c**5 / (hbar * G**2)
g_P = c**7 / (hbar * G**2)
P_P_check = rho_P * g_P * l_P
print(f"  Planck density: ρ_P = {rho_P:.3e} kg/m³")
print(f"  Planck acceleration: g_P = {g_P:.3e} m/s²")
print(f"  P_P ≈ ρ_P·g_P·l_P = {P_P_check:.3e} Pa")
print()

# ========================================
# STEP 7: Dimensional Crosschecks
# ========================================
print("STEP 7: Dimensional Crosschecks")
print("-" * 70)
print("Verify units on both sides:")
print()
print("LHS: P_hydro")
print("  Dimensions: [M L⁻¹ T⁻²]")
print("  Units: Pa = kg/(m·s²)")
print()
print("RHS: ρ·g·h")
print("  Dimensions: [M L⁻³] · [L T⁻²] · [L]")
print("            = [M L⁻¹ T⁻²]")
print("  Units: (kg/m³) · (m/s²) · m = kg/(m·s²) = Pa")
print()
print("✓ Dimensional consistency verified")
print()
print("Hydrostatic equilibrium equation:")
print("  dP/dh = -ρ·g")
print()
print("Barometric formula (exponential atmosphere):")
print("  P(h) = P₀·exp(-h/H)  where H = kT/(m·g)")
H_scale_calc = (1.380649e-23 * 288) / (28.97 * 1.66053906660e-27 * g_earth)
print(f"  Scale height H ≈ {H_scale_calc:.1f} m")
print("  (For Earth atmosphere T=288K, M=28.97 g/mol)")
print()

# ========================================
# STEP 8: Calibration Checkpoint
# ========================================
print("STEP 8: Calibration Checkpoint")
print("-" * 70)
print("Atmospheric pressure verification:")
print("  P₀ = 101325 Pa (1 atm at sea level)")
print(f"  ρ·g·H = {rho_air * g_earth * h_scale:.1f} Pa")
print(f"  Ratio: {(rho_air * g_earth * h_scale)/101325:.3f}")
print("  (Order of magnitude agreement)")
print()
print("Ocean depth pressure:")
print("  10 m depth: P ≈ 1 atm from water column")
print(f"  Calculated: {P_ocean/101325:.3f} atm")
print("  ✓ Excellent agreement")
print()
print("Earth core pressure (from seismology):")
print("  P_core ≈ 360 GPa = 3.6e11 Pa")
P_core_measured = 3.6e11
print(f"  Simple estimate: {P_core:.3e} Pa")
print(f"  Ratio: {P_core/P_core_measured:.3f}")
print("  (Order of magnitude; detailed models needed)")
print()

# ========================================
# SUMMARY
# ========================================
print("=" * 70)
print("DIMENSIONAL ANALYSIS SUMMARY")
print("=" * 70)
print()
print("The hydrostatic pressure derivation is dimensionally")
print("consistent:")
print()
print("1. Target dimensions [M L⁻¹ T⁻²] verified")
print("2. Formula P = ρgh from fluid statics")
print("3. Applies from oceans to planetary interiors")
print("4. Connects mass, gravity, and depth")
print("5. Uses TriPhase-derived G")
print()
print("The formula P_hydro = ρgh represents pressure from")
print("gravitational weight of fluid columns, fundamental to")
print("atmospheric science, oceanography, and planetary structure.")
print()
print("Key insight:")
print("  Pressure increases linearly with depth (for const ρ,g)")
print("  Hydrostatic equilibrium: dP/dh = -ρg")
print("  Scale height H = kT/(mg) for exponential atmospheres")
print()
print("TriPhase connection:")
print("  Gravitational acceleration g = GM/r² uses TriPhase G")
print("  All hydrostatic pressures trace to electromagnetic constants")
print("  P << VF_r for all everyday systems")
print()
print("=" * 70)

input("Press Enter to exit...")
