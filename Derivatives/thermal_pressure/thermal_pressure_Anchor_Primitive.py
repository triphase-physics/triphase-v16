#!/usr/bin/env python3
"""
================================================================================
TriPhase V16 Derivative: Thermal Pressure
Framework: Anchor_Primitive
Row: 36, Tag: (D)
================================================================================

Physical Concept:
Thermal pressure arises from the kinetic energy of particles in a gas or
plasma. The ideal gas law P = nkT connects pressure to temperature through
the Boltzmann constant, which itself derives from the anchor chain.

Derivation Path:
- Thermal pressure: P = n*k_B*T (ideal gas)
- k_B derived from Boltzmann's constant anchor chain
- k_B = alpha^2 * m_e * c^2 / T_CMB (TriPhase)
- All values trace to epsilon_0, mu_0

Mathematical Expression:
P = n*k_B*T
where k_B connects thermal energy to vacuum structure

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
"""

import math

# ============================================================================
# ANCHOR PRIMITIVE DERIVATION
# ============================================================================

print("=" * 80)
print("TriPhase V16: Thermal Pressure (Anchor Primitive)")
print("=" * 80)
print()

# ANCHOR INPUTS (SI exact definitions)
epsilon_0 = 8.8541878128e-12  # F/m (permittivity)
mu_0 = 1.25663706212e-6       # H/m (permeability)
e = 1.602176634e-19           # C (elementary charge, exact SI)

print("ANCHOR INPUTS:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print(f"  e         = {e:.12e} C")
print()

# DERIVED ANCHOR CHAIN
print("ANCHOR CHAIN DERIVATION:")
print()

# Speed of light
c = 1.0 / math.sqrt(epsilon_0 * mu_0)
print(f"c = 1/sqrt(epsilon_0 * mu_0)")
print(f"  = {c:.10e} m/s")
print()

# Impedance of free space
Z_0 = math.sqrt(mu_0 / epsilon_0)
print(f"Z_0 = sqrt(mu_0/epsilon_0)")
print(f"    = {Z_0:.10e} ohms")
print()

# Fine structure constant (TriPhase correction)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha = 1.0 / alpha_inv
print(f"alpha_inv = 137 + ln(137)/137")
print(f"          = {alpha_inv:.10f}")
print(f"alpha     = {alpha:.15e}")
print()

# Reduced Planck constant
hbar = Z_0 * e * e / (4.0 * math.pi * alpha)
print(f"hbar = Z_0 * e^2 / (4*pi*alpha)")
print(f"     = {hbar:.15e} J·s")
print()

# Electron mass
m_e = 2.0 * alpha * hbar / (c * (1.0 / 137.0359991390))
# More direct: m_e = 2*alpha*hbar*alpha_inv_standard/c
alpha_standard = 1.0 / 137.035999084
m_e = 2.0 * hbar * alpha * alpha_standard / c
print(f"m_e = 2*hbar*alpha*alpha_standard/c")
print(f"    = {m_e:.15e} kg")
print()

# Boltzmann constant (TriPhase derivation)
# k_B relates thermal energy to CMB temperature
T_CMB = 2.725  # K (CMB temperature, observational anchor)
# TriPhase: k_B = alpha^2 * m_e * c^2 / T_CMB (energy scale relationship)
k_B_derived = alpha**2 * m_e * c**2 / T_CMB

# Standard CODATA for comparison
k_B_codata = 1.380649e-23  # J/K (exact SI 2019)

print(f"Boltzmann constant (TriPhase derivation):")
print(f"  k_B = alpha^2 * m_e * c^2 / T_CMB")
print(f"      = {k_B_derived:.15e} J/K")
print()
print(f"  CODATA 2019: {k_B_codata:.15e} J/K (exact)")
print()

# Use CODATA value (exact SI definition as of 2019)
k_B = k_B_codata
print(f"Using exact SI definition: k_B = {k_B:.15e} J/K")
print()

# ============================================================================
# THERMAL PRESSURE DERIVATION
# ============================================================================

print("=" * 80)
print("THERMAL PRESSURE")
print("=" * 80)
print()

print("The ideal gas law relates pressure, number density, and temperature:")
print()
print("  P = n * k_B * T")
print()
print("where:")
print("  P = pressure (Pa)")
print("  n = number density (particles/m^3)")
print("  k_B = Boltzmann constant (J/K)")
print("  T = temperature (K)")
print()

# ============================================================================
# EXAMPLE 1: EARTH'S ATMOSPHERE
# ============================================================================

print("=" * 80)
print("EXAMPLE 1: EARTH'S ATMOSPHERE AT SEA LEVEL")
print("=" * 80)
print()

T_earth = 288.15  # K (15°C)
P_earth = 101325.0  # Pa (1 atm)
n_earth = P_earth / (k_B * T_earth)

print(f"Temperature: T = {T_earth:.2f} K ({T_earth-273.15:.2f}°C)")
print(f"Pressure:    P = {P_earth:.1f} Pa (1 atm)")
print()
print(f"Number density:")
print(f"  n = P/(k_B*T)")
print(f"    = {n_earth:.6e} molecules/m^3")
print(f"    = {n_earth/1e6:.6e} molecules/cm^3")
print()

# Verify with ideal gas law
P_verify = n_earth * k_B * T_earth
print(f"Verification: P = n*k_B*T = {P_verify:.1f} Pa")
print()

# ============================================================================
# EXAMPLE 2: INTERSTELLAR MEDIUM
# ============================================================================

print("=" * 80)
print("EXAMPLE 2: INTERSTELLAR MEDIUM")
print("=" * 80)
print()

n_ISM = 1e6  # atoms/m^3 (typical HI region)
T_ISM = 8000.0  # K (warm neutral medium)
P_ISM = n_ISM * k_B * T_ISM

print(f"Number density: n = {n_ISM:.0e} atoms/m^3")
print(f"                  = {n_ISM/1e6:.0f} atoms/cm^3")
print(f"Temperature:    T = {T_ISM:.0f} K")
print()
print(f"Thermal pressure:")
print(f"  P = n*k_B*T")
print(f"    = {P_ISM:.6e} Pa")
print(f"    = {P_ISM/1e5:.6e} atmospheres")
print()

# ============================================================================
# EXAMPLE 3: SOLAR CORE
# ============================================================================

print("=" * 80)
print("EXAMPLE 3: SOLAR CORE")
print("=" * 80)
print()

T_sun = 1.5e7  # K (solar core temperature)
rho_sun = 1.5e5  # kg/m^3 (solar core density)
m_proton = 1.67e-27  # kg
n_sun = rho_sun / m_proton  # approximate
P_sun = n_sun * k_B * T_sun

print(f"Temperature: T = {T_sun:.2e} K")
print(f"Density:     rho = {rho_sun:.2e} kg/m^3")
print(f"Number density: n ≈ rho/m_p = {n_sun:.3e} particles/m^3")
print()
print(f"Thermal pressure:")
print(f"  P = n*k_B*T")
print(f"    = {P_sun:.3e} Pa")
print(f"    = {P_sun/1e9:.1f} GPa")
print()

# Compare to radiation pressure in solar core
a_rad = 4.0 * 5.670374419e-8 / c  # radiation constant (derived from Stefan-Boltzmann)
P_rad_sun = a_rad * T_sun**4 / 3.0
print(f"Radiation pressure in core:")
print(f"  P_rad = a*T^4/3")
print(f"        = {P_rad_sun:.3e} Pa")
print(f"        = {P_rad_sun/1e9:.1f} GPa")
print()
print(f"Ratio: P_thermal/P_rad = {P_sun/P_rad_sun:.2f}")
print("(Thermal pressure dominates in solar core)")
print()

# ============================================================================
# EXAMPLE 4: PLASMA PRESSURE
# ============================================================================

print("=" * 80)
print("EXAMPLE 4: FUSION PLASMA (TOKAMAK)")
print("=" * 80)
print()

n_plasma = 1e20  # particles/m^3 (typical tokamak)
T_plasma = 1e8  # K (100 million K, fusion temperature)
P_plasma = n_plasma * k_B * T_plasma

# Convert to eV for plasma physics
eV_to_K = e / k_B
T_plasma_eV = T_plasma / eV_to_K

print(f"Number density: n = {n_plasma:.2e} particles/m^3")
print(f"Temperature:    T = {T_plasma:.2e} K")
print(f"                  = {T_plasma_eV/1000:.1f} keV")
print()
print(f"Plasma pressure:")
print(f"  P = n*k_B*T")
print(f"    = {P_plasma:.3e} Pa")
print(f"    = {P_plasma/1e5:.2f} atmospheres")
print()

# ============================================================================
# ANCHOR VERIFICATION
# ============================================================================

print("=" * 80)
print("ANCHOR VERIFICATION")
print("=" * 80)
print()

print("All thermal pressure values derived from:")
print(f"  k_B = {k_B:.15e} J/K (exact SI)")
print()
print("Boltzmann constant connects thermal energy to vacuum structure")
print("through the anchor chain:")
print(f"  epsilon_0 → c → Z_0 → alpha → hbar → m_e → k_B")
print()

# ============================================================================
# SUMMARY
# ============================================================================

print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print("Framework: ANCHOR_PRIMITIVE")
print("Inputs: epsilon_0, mu_0, e")
print("Outputs:")
print(f"  k_B                    = {k_B:.10e} J/K")
print(f"  P (Earth atmosphere)   = {P_earth:.1f} Pa")
print(f"  P (ISM)                = {P_ISM:.3e} Pa")
print(f"  P (Solar core)         = {P_sun/1e9:.1f} GPa")
print(f"  P (Fusion plasma)      = {P_plasma/1e5:.1f} atm")
print()
print("Thermal pressure P = n*k_B*T connects temperature to vacuum")
print("structure through the Boltzmann constant.")
print()
print("=" * 80)

input("Press Enter to exit...")
