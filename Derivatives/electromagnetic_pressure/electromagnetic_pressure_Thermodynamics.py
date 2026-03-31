"""
================================================================================
TriPhase V16 - Electromagnetic Pressure Derivative
Framework: THERMODYNAMICS
Tag: (D) - Pure derivation
================================================================================

THERMODYNAMICS FRAMEWORK:
Interprets each physical quantity through statistical mechanics / thermodynamic
concepts: partition functions Z and free energy F = -k_BT ln(Z), entropy
S = -∂F/∂T, equipartition theorem (½k_BT per degree of freedom), Boltzmann
distributions, Maxwell-Boltzmann statistics, phase transitions, order parameters,
critical phenomena, equations of state, thermodynamic potentials (U, H, F, G),
degrees of freedom counting, mode counting, heat capacity, specific heat,
adiabatic processes, black-body radiation, Planck distribution, thermodynamic
stability conditions, Stefan-Boltzmann law, Wien displacement, chemical potential,
Gibbs free energy.

PHYSICAL DERIVATION:
Electromagnetic radiation exerts pressure — this is radiation pressure, a
fundamental thermodynamic quantity. The EM stress tensor gives:

    P_EM = u_EM = ε₀E²/2 + B²/(2μ₀)

For isotropic radiation in thermal equilibrium, the pressure is:
    P_rad = u_rad / 3

where u_rad is the radiation energy density. This 1/3 factor comes from
averaging the momentum flux over all directions.

The Stefan-Boltzmann law gives:
    u_rad = aT⁴ where a = 4σ/c = (π²k_B⁴)/(15ℏ³c³)

Therefore:
    P_rad = (aT⁴)/3 = (π²/45)(k_BT/ℏc)³ × k_BT

This is the thermodynamic equation of state for a photon gas.

Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# ============================================================================
# STANDARD ANCHOR CHAIN
# ============================================================================
epsilon_0 = 8.8541878128e-12   # F/m
mu_0      = 1.25663706212e-6   # H/m
e         = 1.602176634e-19    # C (exact SI)
c         = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0       = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha     = 1.0 / alpha_inv
hbar      = Z_0 * e**2 / (4.0 * math.pi * alpha)
h         = 2.0 * math.pi * hbar
G         = c**4 * 7.5 * epsilon_0**3 * mu_0**2
r_e       = 2.8179403262e-15
m_e       = hbar * alpha / (c * r_e)
f_e       = m_e * c**2 / hbar
T_17      = 17 * 18 // 2   # = 153
mp_me     = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p       = m_e * mp_me
H_0       = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r      = c**4 / (8.0 * math.pi * G)

# ============================================================================
# THERMODYNAMIC DERIVATION - ELECTROMAGNETIC PRESSURE
# ============================================================================

print("=" * 80)
print("ELECTROMAGNETIC PRESSURE - THERMODYNAMICS FRAMEWORK")
print("=" * 80)
print()
print("THERMODYNAMIC INTERPRETATION:")
print("EM radiation pressure P_rad = u/3 is the equation of state for a")
print("photon gas. This is a fundamental thermodynamic relation for massless")
print("bosons in thermal equilibrium.")
print()

# EM energy density
print("ELECTROMAGNETIC ENERGY DENSITY:")
print(f"  u_EM = ε₀E²/2 + B²/(2μ₀)")
print(f"  For plane waves: ε₀E² = B²/μ₀ → u_EM = ε₀E² = B²/μ₀")
print()

# Radiation pressure (general)
print("RADIATION PRESSURE (INSTANTANEOUS):")
print(f"  P_EM = u_EM (pressure = energy density for EM fields)")
print(f"  But for isotropic radiation averaged over directions:")
print(f"  P_rad = u_rad / 3")
print()

# Stefan-Boltzmann constant
k_B = 1.380649e-23  # J/K
sigma_SB = (2.0 * math.pi**5 * k_B**4) / (15.0 * h**3 * c**2)
a_rad = 4.0 * sigma_SB / c  # Radiation constant

print("STEFAN-BOLTZMANN LAW:")
print(f"  u_rad = aT⁴")
print(f"  where a = 4σ/c = {a_rad:.6e} J/(m³·K⁴)")
print(f"  σ_SB = {sigma_SB:.6e} W/(m²·K⁴)")
print()

# ============================================================================
# THERMODYNAMIC DERIVATION FROM PARTITION FUNCTION
# ============================================================================

print("=" * 80)
print("THERMODYNAMIC DERIVATION FROM PARTITION FUNCTION")
print("=" * 80)
print()

print("PHOTON GAS PARTITION FUNCTION:")
print(f"  For photons (massless bosons), the partition function is:")
print(f"    ln(Z) = -(V/π²ℏ³c³) ∫₀^∞ k² ln(1 - e^(-ℏck/k_BT)) dk")
print()
print(f"  This integral gives:")
print(f"    F = -k_BT ln(Z) = -(π²V/45)(k_BT)⁴/(ℏc)³")
print()

print("THERMODYNAMIC QUANTITIES:")
print(f"  Energy: U = -∂(ln Z)/∂β|_V = (π²V/15)(k_BT)⁴/(ℏc)³")
print(f"  Pressure: P = -∂F/∂V|_T = U/(3V) = aT⁴/3")
print(f"  Entropy: S = -∂F/∂T|_V = (4/3)aVT³")
print()

# Temperature dependence
T_CMB = 2.725  # K (CMB temperature today)
u_CMB = a_rad * T_CMB**4
P_CMB = u_CMB / 3.0

print("COSMIC MICROWAVE BACKGROUND:")
print(f"  Temperature T_CMB           = {T_CMB:.3f} K")
print(f"  Energy density u_CMB        = {u_CMB:.6e} J/m³")
print(f"  Radiation pressure P_CMB    = {P_CMB:.6e} Pa")
print()

# Solar radiation
T_sun_surface = 5778.0  # K
u_sun = a_rad * T_sun_surface**4
P_sun_rad = u_sun / 3.0

print("SOLAR SURFACE RADIATION:")
print(f"  Temperature T_sun           = {T_sun_surface:.0f} K")
print(f"  Energy density u_sun        = {u_sun:.6e} J/m³")
print(f"  Radiation pressure P_sun    = {P_sun_rad:.6e} Pa")
print()

# Solar radiation pressure at Earth
solar_constant = 1361.0  # W/m² (solar irradiance at Earth)
P_solar_Earth = solar_constant / c

print("SOLAR RADIATION PRESSURE AT EARTH:")
print(f"  Solar constant S            = {solar_constant:.1f} W/m²")
print(f"  Radiation pressure P_rad    = S/c")
print(f"                              = {P_solar_Earth:.6e} Pa")
print(f"                              ≈ {P_solar_Earth*1e6:.2f} μPa")
print()
print(f"  This is the force per area on a solar sail!")
print()

# ============================================================================
# PLANCK DISTRIBUTION
# ============================================================================

print("=" * 80)
print("PLANCK DISTRIBUTION (THERMODYNAMIC SPECTRUM)")
print("=" * 80)
print()

print("ENERGY DENSITY PER FREQUENCY:")
print(f"  u_ν(ν,T) = (8πh/c³) × ν³ / (e^(hν/k_BT) - 1)")
print()
print(f"  This is the Planck blackbody spectrum — the thermodynamic")
print(f"  equilibrium distribution for photons at temperature T.")
print()

# Wien displacement law
b_Wien = h * c / (4.965 * k_B)  # Wien constant
lambda_max_CMB = b_Wien / T_CMB
lambda_max_sun = b_Wien / T_sun_surface

print("WIEN DISPLACEMENT LAW:")
print(f"  λ_max × T = constant ≈ {b_Wien*1000:.3f} mm·K")
print()
print(f"  CMB peak wavelength λ_max   = {lambda_max_CMB*1000:.2f} mm")
print(f"  Solar peak wavelength λ_max = {lambda_max_sun*1e9:.0f} nm (green)")
print()

# ============================================================================
# APPLICATIONS
# ============================================================================

print("=" * 80)
print("THERMODYNAMIC APPLICATIONS")
print("=" * 80)
print()

print("1. STELLAR INTERIORS:")
print(f"  In hot stars, radiation pressure P_rad can exceed gas pressure P_gas.")
print(f"  Eddington limit: Maximum luminosity when P_rad balances gravity.")
print(f"  For T = 10⁷ K (solar core):")
T_solar_core = 1.57e7  # K
P_rad_core = a_rad * T_solar_core**4 / 3.0
print(f"    P_rad ≈ {P_rad_core:.3e} Pa")
print()

print("2. EARLY UNIVERSE:")
print(f"  In the radiation-dominated era (T > 3000 K):")
print(f"    ρ_rad ∝ T⁴, P_rad = ρ_rad/3")
print(f"  Friedmann equation: H² ∝ ρ_rad ∝ T⁴")
print(f"  → T ∝ 1/a (temperature scales inversely with scale factor)")
print()

print("3. COSMIC ACCELERATION:")
print(f"  Radiation has w = P/ρ = 1/3 → decelerates expansion")
print(f"  Matter has w = 0 → also decelerates")
print(f"  Dark energy has w = -1 → accelerates expansion!")
print()

# ============================================================================
# MOMENTUM TRANSFER
# ============================================================================

print("=" * 80)
print("RADIATION PRESSURE AS MOMENTUM TRANSFER")
print("=" * 80)
print()

print("PHOTON MOMENTUM:")
print(f"  Each photon carries momentum p = E/c = hν/c")
print(f"  For N photons per second per area:")
print(f"    Momentum flux = N × (hν/c)")
print()

print("FORCE ON ABSORBING SURFACE:")
print(f"  F/A = (energy flux) / c = I/c")
print(f"  where I is intensity (W/m²)")
print()

print("FORCE ON REFLECTING SURFACE:")
print(f"  F/A = 2I/c (factor of 2 from reflection)")
print()

# Solar sail example
area_sail = 100.0  # m²
F_solar_sail = P_solar_Earth * area_sail * 2.0  # Reflecting surface
mass_sail = 10.0  # kg
accel_sail = F_solar_sail / mass_sail

print("SOLAR SAIL EXAMPLE:")
print(f"  Sail area A                 = {area_sail:.0f} m²")
print(f"  Sail mass m                 = {mass_sail:.0f} kg")
print(f"  Force F = 2 × P_rad × A     = {F_solar_sail:.6e} N")
print(f"  Acceleration a = F/m        = {accel_sail:.6e} m/s²")
print(f"                              = {accel_sail*1e3:.3f} mm/s²")
print()
print(f"  After 1 year:")
delta_v = accel_sail * 365.25 * 24 * 3600
print(f"    Δv ≈ {delta_v:.0f} m/s ≈ {delta_v/1000:.1f} km/s")
print()

# ============================================================================
# CALIBRATION CHECKPOINT
# ============================================================================

print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
print()

sigma_SB_PDG = 5.670374419e-8  # W/(m²·K⁴)
a_rad_PDG = 4.0 * sigma_SB_PDG / c

deviation_sigma = ((sigma_SB - sigma_SB_PDG) / sigma_SB_PDG) * 100
deviation_a = ((a_rad - a_rad_PDG) / a_rad_PDG) * 100

print(f"TriPhase σ_SB               = {sigma_SB:.6e} W/(m²·K⁴)")
print(f"CODATA 2018 σ_SB            = {sigma_SB_PDG:.9e} W/(m²·K⁴)")
print(f"Deviation                   = {deviation_sigma:+.4f}%")
print()
print(f"TriPhase a_rad              = {a_rad:.6e} J/(m³·K⁴)")
print(f"CODATA a_rad                = {a_rad_PDG:.6e} J/(m³·K⁴)")
print(f"Deviation                   = {deviation_a:+.4f}%")
print()

print("MEASURED RADIATION PRESSURE:")
print(f"  Solar constant at Earth     = {solar_constant:.1f} W/m² (satellite measured)")
print(f"  CMB temperature             = {T_CMB:.3f} K (COBE/WMAP/Planck)")
print(f"  Both confirm P_rad = u/3 equation of state.")
print()

print("THERMODYNAMIC SIGNIFICANCE:")
print("  Radiation pressure is a pure thermodynamic quantity:")
print("    - Emerges from photon partition function")
print("    - Equation of state: P = u/3 (massless bosons)")
print("    - Stefan-Boltzmann: u ∝ T⁴ (blackbody thermodynamics)")
print("  This drives cosmic expansion in early universe.")
print()

print("=" * 80)
print("END ELECTROMAGNETIC PRESSURE THERMODYNAMICS DERIVATIVE")
print("=" * 80)

input("Press Enter to exit...")
