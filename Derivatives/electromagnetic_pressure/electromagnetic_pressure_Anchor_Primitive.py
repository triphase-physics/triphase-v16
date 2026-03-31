#!/usr/bin/env python3
"""
================================================================================
TriPhase V16 Derivative: Electromagnetic Pressure
Framework: Anchor_Primitive
Row: 35, Tag: (D)
================================================================================

Physical Concept:
Electromagnetic fields exert pressure on matter through the Maxwell stress
tensor. This pressure is fundamental to understanding radiation pressure,
vacuum energy, and field interactions.

Derivation Path:
- EM pressure: P_EM = epsilon_0*E^2/2 + B^2/(2*mu_0)
- Electric energy density: u_E = epsilon_0*E^2/2
- Magnetic energy density: u_B = B^2/(2*mu_0)
- Total pressure = total energy density (for radiation)

Mathematical Expression:
P_EM = epsilon_0*E^2/2 + B^2/(2*mu_0)
All pressure traces directly to epsilon_0, mu_0

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
"""

import math

# ============================================================================
# ANCHOR PRIMITIVE DERIVATION
# ============================================================================

print("=" * 80)
print("TriPhase V16: Electromagnetic Pressure (Anchor Primitive)")
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

# ============================================================================
# ELECTROMAGNETIC PRESSURE DERIVATION
# ============================================================================

print("=" * 80)
print("ELECTROMAGNETIC PRESSURE")
print("=" * 80)
print()

print("Electromagnetic fields store energy and exert pressure:")
print()
print("Electric energy density:")
print("  u_E = epsilon_0 * E^2 / 2")
print()
print("Magnetic energy density:")
print("  u_B = B^2 / (2*mu_0)")
print()
print("Total EM pressure (for radiation):")
print("  P_EM = epsilon_0*E^2/2 + B^2/(2*mu_0)")
print()
print("For electromagnetic waves: E = c*B and u_E = u_B")
print()

# ============================================================================
# EXAMPLE CALCULATIONS
# ============================================================================

print("=" * 80)
print("EXAMPLE 1: SUNLIGHT AT EARTH")
print("=" * 80)
print()

# Solar constant
I_solar = 1361.0  # W/m^2 (solar irradiance at Earth)
print(f"Solar irradiance: I = {I_solar:.1f} W/m^2")
print()

# For EM wave: I = (1/2)*epsilon_0*c*E_0^2
E_0_solar = math.sqrt(2.0 * I_solar / (epsilon_0 * c))
B_0_solar = E_0_solar / c

print(f"Peak electric field:")
print(f"  E_0 = sqrt(2*I/(epsilon_0*c))")
print(f"      = {E_0_solar:.3f} V/m")
print()
print(f"Peak magnetic field:")
print(f"  B_0 = E_0/c")
print(f"      = {B_0_solar:.10e} T")
print(f"      = {B_0_solar*1e9:.6f} nT")
print()

# Energy density (time-averaged)
u_E_avg = epsilon_0 * E_0_solar**2 / 4.0
u_B_avg = B_0_solar**2 / (4.0 * mu_0)
u_total = u_E_avg + u_B_avg

print(f"Time-averaged energy density:")
print(f"  u_E = epsilon_0*E_0^2/4 = {u_E_avg:.10e} J/m^3")
print(f"  u_B = B_0^2/(4*mu_0)    = {u_B_avg:.10e} J/m^3")
print(f"  u_total                 = {u_total:.10e} J/m^3")
print()

# Radiation pressure
P_rad = u_total  # For normal incidence
P_rad_absorbed = I_solar / c  # Alternative calculation

print(f"Radiation pressure (perfect absorption):")
print(f"  P_rad = u_total = {P_rad:.10e} Pa")
print(f"  Also: P_rad = I/c = {P_rad_absorbed:.10e} Pa")
print()

# Verification
print(f"Verification: {abs(P_rad - P_rad_absorbed):.3e} Pa difference")
print()

# Force on 1 m^2 sail
F_sail = P_rad * 1.0  # N
print(f"Force on 1 m^2 solar sail:")
print(f"  F = {F_sail:.10e} N")
print(f"    = {F_sail*1e6:.3f} µN")
print()

# ============================================================================
# EXAMPLE 2: STRONG LASER FIELD
# ============================================================================

print("=" * 80)
print("EXAMPLE 2: HIGH-POWER LASER")
print("=" * 80)
print()

P_laser = 1e15  # W (petawatt laser)
A_laser = math.pi * (1e-6)**2  # m^2 (1 micron radius spot)
I_laser = P_laser / A_laser

print(f"Laser power: P = {P_laser:.3e} W")
print(f"Spot area:   A = {A_laser:.3e} m^2")
print(f"Intensity:   I = {I_laser:.3e} W/m^2")
print()

E_0_laser = math.sqrt(2.0 * I_laser / (epsilon_0 * c))
B_0_laser = E_0_laser / c

print(f"Peak electric field:")
print(f"  E_0 = {E_0_laser:.3e} V/m")
print(f"      = {E_0_laser/1e9:.3f} GV/m")
print()
print(f"Peak magnetic field:")
print(f"  B_0 = {B_0_laser:.3e} T")
print(f"      = {B_0_laser/1e4:.3f} tesla (10^4 T)")
print()

# Energy density
u_laser = I_laser / c
print(f"Energy density:")
print(f"  u = I/c = {u_laser:.3e} J/m^3")
print()

# Pressure
P_laser_rad = u_laser
print(f"Radiation pressure:")
print(f"  P = {P_laser_rad:.3e} Pa")
print(f"    = {P_laser_rad/1e9:.3f} GPa")
print()

# ============================================================================
# EXAMPLE 3: UNIT FIELD PRESSURES
# ============================================================================

print("=" * 80)
print("EXAMPLE 3: UNIT FIELD PRESSURES")
print("=" * 80)
print()

print("Pressure from 1 V/m electric field:")
E_unit = 1.0  # V/m
P_E_unit = epsilon_0 * E_unit**2 / 2.0
print(f"  P_E = epsilon_0*E^2/2")
print(f"      = {P_E_unit:.15e} Pa")
print()

print("Pressure from 1 T magnetic field:")
B_unit = 1.0  # T
P_B_unit = B_unit**2 / (2.0 * mu_0)
print(f"  P_B = B^2/(2*mu_0)")
print(f"      = {P_B_unit:.15e} Pa")
print(f"      = {P_B_unit/1e5:.3f} atmospheres")
print()

print("Ratio:")
ratio = P_B_unit / P_E_unit
print(f"  P_B/P_E = (B^2/(2*mu_0)) / (epsilon_0*E^2/2)")
print(f"          = B^2/(mu_0*epsilon_0*E^2)")
print(f"          = (B/(E/c))^2  (since c = 1/sqrt(mu_0*epsilon_0))")
print(f"          = {ratio:.3e}")
print()
print("Magnetic fields produce much stronger pressure than electric fields")
print("of numerically equal magnitude.")
print()

# ============================================================================
# ANCHOR VERIFICATION
# ============================================================================

print("=" * 80)
print("ANCHOR VERIFICATION")
print("=" * 80)
print()

print("All pressures derived from epsilon_0, mu_0 only:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print(f"  c         = {c:.10e} m/s")
print()

print("Electric pressure scales as epsilon_0*E^2")
print("Magnetic pressure scales as B^2/mu_0")
print("Both trace directly to vacuum structure.")
print()

# ============================================================================
# SUMMARY
# ============================================================================

print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print("Framework: ANCHOR_PRIMITIVE")
print("Inputs: epsilon_0, mu_0")
print("Outputs:")
print(f"  Solar radiation pressure = {P_rad:.6e} Pa")
print(f"  Laser radiation pressure = {P_laser_rad:.3e} Pa ({P_laser_rad/1e9:.3f} GPa)")
print(f"  1 V/m electric pressure  = {P_E_unit:.6e} Pa")
print(f"  1 T magnetic pressure    = {P_B_unit:.3e} Pa")
print()
print("All electromagnetic pressure traces to epsilon_0, mu_0.")
print()
print("=" * 80)

input("Press Enter to exit...")
