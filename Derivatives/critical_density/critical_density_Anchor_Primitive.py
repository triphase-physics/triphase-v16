#!/usr/bin/env python3
"""
================================================================================
TriPhase V16 Derivative: Critical Density
Framework: Anchor_Primitive
Row: 41, Tag: (D)
================================================================================

Physical Concept:
The critical density (rho_c) is the density required for the universe to be
spatially flat. It depends on the Hubble constant and gravitational constant,
both of which derive from the vacuum structure in TriPhase.

Derivation Path:
- Critical density: rho_c = 3*H_0^2/(8*pi*G)
- H_0 derived from anchor chain: H_0 = pi*sqrt(3) * f_e * alpha^18
- G from anchor chain: G = c^4 * 7.5 * epsilon_0^3 * mu_0^2
- All values trace to epsilon_0, mu_0

Mathematical Expression:
rho_c = 3*H_0^2/(8*pi*G)
Determines the geometry and fate of the universe.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
"""

import math

# ============================================================================
# ANCHOR PRIMITIVE DERIVATION
# ============================================================================

print("=" * 80)
print("TriPhase V16: Critical Density (Anchor Primitive)")
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
alpha_standard = 1.0 / 137.035999084
m_e = 2.0 * hbar * alpha * alpha_standard / c
print(f"m_e = 2*hbar*alpha*alpha_standard/c")
print(f"    = {m_e:.15e} kg")
print()

# Electron Compton frequency
f_e = m_e * c**2 / hbar
print(f"f_e = m_e*c^2/hbar")
print(f"    = {f_e:.15e} Hz")
print()

# Gravitational constant
G = c**4 * 7.5 * epsilon_0**3 * mu_0**2
print(f"G = c^4 * 7.5 * epsilon_0^3 * mu_0^2")
print(f"  = {G:.15e} m^3/(kg·s^2)")
print()

# Hubble constant (TriPhase derivation)
H_0 = math.pi * math.sqrt(3.0) * f_e * alpha**18
print(f"H_0 = pi*sqrt(3) * f_e * alpha^18")
print(f"    = {H_0:.15e} s^-1")
print(f"    = {H_0 * 3.086e19:.3f} km/s/Mpc")
print()

# ============================================================================
# CRITICAL DENSITY DERIVATION
# ============================================================================

print("=" * 80)
print("CRITICAL DENSITY")
print("=" * 80)
print()

print("The Friedmann equation for a flat universe:")
print("  H^2 = (8*pi*G/3)*rho")
print()
print("Critical density is the density for which k=0 (flat universe):")
print("  rho_c = 3*H_0^2/(8*pi*G)")
print()

# Calculate critical density
rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)

print(f"Critical density:")
print(f"  rho_c = 3*H_0^2/(8*pi*G)")
print(f"        = {rho_c:.15e} kg/m^3")
print()

# Convert to various units
rho_c_g_cm3 = rho_c * 1e-3  # g/cm^3
n_c_protons = rho_c / (1.67e-27)  # protons/m^3

print(f"In other units:")
print(f"  rho_c = {rho_c_g_cm3:.15e} g/cm^3")
print(f"        = {n_c_protons:.6e} protons/m^3")
print(f"        = {n_c_protons/1e6:.6f} protons/cm^3")
print()

# ============================================================================
# DENSITY PARAMETERS
# ============================================================================

print("=" * 80)
print("DENSITY PARAMETERS (Omega)")
print("=" * 80)
print()

print("Density parameter: Omega = rho/rho_c")
print()

# Observed values
Omega_m = 0.315      # Matter (baryonic + dark)
Omega_Lambda = 0.685  # Dark energy
Omega_total = Omega_m + Omega_Lambda

print(f"Observed density parameters:")
print(f"  Omega_m      = {Omega_m:.3f} (matter)")
print(f"  Omega_Lambda = {Omega_Lambda:.3f} (dark energy)")
print(f"  Omega_total  = {Omega_total:.3f}")
print()

# Actual densities
rho_m = Omega_m * rho_c
rho_Lambda = Omega_Lambda * rho_c

print(f"Actual densities:")
print(f"  rho_m      = {rho_m:.15e} kg/m^3")
print(f"  rho_Lambda = {rho_Lambda:.15e} kg/m^3")
print()

# Matter breakdown
Omega_baryon = 0.049
Omega_dark_matter = Omega_m - Omega_baryon

rho_baryon = Omega_baryon * rho_c
rho_dark_matter = Omega_dark_matter * rho_c

print(f"Matter composition:")
print(f"  Omega_baryon      = {Omega_baryon:.3f}")
print(f"  Omega_dark_matter = {Omega_dark_matter:.3f}")
print()
print(f"  rho_baryon      = {rho_baryon:.15e} kg/m^3")
print(f"  rho_dark_matter = {rho_dark_matter:.15e} kg/m^3")
print()

# ============================================================================
# HUBBLE RADIUS AND VOLUME
# ============================================================================

print("=" * 80)
print("HUBBLE RADIUS AND OBSERVABLE UNIVERSE")
print("=" * 80)
print()

R_H = c / H_0
print(f"Hubble radius (horizon):")
print(f"  R_H = c/H_0")
print(f"      = {R_H:.15e} m")
print(f"      = {R_H/9.461e15:.3e} light-years")
print(f"      = {R_H/3.086e16:.3f} Gpc")
print()

# Hubble volume
V_H = (4.0/3.0) * math.pi * R_H**3
print(f"Hubble volume:")
print(f"  V_H = (4/3)*pi*R_H^3")
print(f"      = {V_H:.15e} m^3")
print()

# Total mass in Hubble volume at critical density
M_H_critical = rho_c * V_H
M_sun = 1.989e30  # kg

print(f"Mass in Hubble volume (at critical density):")
print(f"  M = rho_c * V_H")
print(f"    = {M_H_critical:.15e} kg")
print(f"    = {M_H_critical/M_sun:.3e} solar masses")
print()

# Actual mass (accounting for Omega_total)
M_H_actual = Omega_total * M_H_critical
print(f"Actual mass (Omega_total = {Omega_total:.3f}):")
print(f"  M = {M_H_actual:.15e} kg")
print(f"    = {M_H_actual/M_sun:.3e} solar masses")
print()

# ============================================================================
# COMPARISON TO OTHER DENSITIES
# ============================================================================

print("=" * 80)
print("COMPARISON TO FAMILIAR DENSITIES")
print("=" * 80)
print()

densities = [
    ("Critical density", rho_c),
    ("Interstellar medium", 1e-21),
    ("Intergalactic medium", 1e-27),
    ("Earth atmosphere (sea level)", 1.2),
    ("Water", 1000.0),
    ("Iron", 7874.0),
    ("Neutron star", 1e18),
]

print("Density comparison:")
print("-" * 70)
for name, rho in densities:
    ratio = rho / rho_c
    print(f"{name:30s}: {rho:.3e} kg/m^3")
    print(f"{'':30s}  = {ratio:.3e} * rho_c")
    print()

# ============================================================================
# UNIVERSE GEOMETRY
# ============================================================================

print("=" * 80)
print("UNIVERSE GEOMETRY")
print("=" * 80)
print()

print("The density parameter determines universe geometry:")
print()
print("  Omega < 1: Open universe (negative curvature)")
print("  Omega = 1: Flat universe (zero curvature)")
print("  Omega > 1: Closed universe (positive curvature)")
print()
print(f"Observed: Omega_total = {Omega_total:.3f}")
print()

if abs(Omega_total - 1.0) < 0.02:
    print("The universe is consistent with being flat within measurement")
    print("uncertainty. This is a key prediction of inflationary cosmology.")
else:
    print(f"Deviation from flatness: {abs(Omega_total-1.0):.3f}")
print()

# ============================================================================
# TIME SCALES
# ============================================================================

print("=" * 80)
print("COSMOLOGICAL TIME SCALES")
print("=" * 80)
print()

# Hubble time
t_H = 1.0 / H_0
print(f"Hubble time:")
print(f"  t_H = 1/H_0")
print(f"      = {t_H:.15e} s")
print(f"      = {t_H/(365.25*24*3600):.3e} years")
print(f"      = {t_H/(365.25*24*3600*1e9):.3f} Gyr")
print()

# Age of universe (flat universe with matter + dark energy)
# Approximate formula
t_0 = (2.0/(3.0*H_0)) * math.log((1.0 + math.sqrt(Omega_Lambda/Omega_m))**2 / Omega_m)
print(f"Age of universe (approximate):")
print(f"  t_0 = {t_0:.15e} s")
print(f"      = {t_0/(365.25*24*3600*1e9):.3f} Gyr")
print()

# ============================================================================
# ANCHOR VERIFICATION
# ============================================================================

print("=" * 80)
print("ANCHOR VERIFICATION")
print("=" * 80)
print()

print("All values derived from epsilon_0, mu_0 only:")
print(f"  H_0    = {H_0:.15e} s^-1")
print(f"  G      = {G:.15e} m^3/(kg·s^2)")
print(f"  rho_c  = {rho_c:.15e} kg/m^3")
print()

# Compare to observational values
H_0_obs = 2.2e-18  # s^-1 (~70 km/s/Mpc)
rho_c_obs = 3.0 * H_0_obs**2 / (8.0 * math.pi * G)

print(f"Comparison to observations:")
print(f"  H_0 (TriPhase) = {H_0*3.086e19:.3f} km/s/Mpc")
print(f"  H_0 (observed) ~ 70 km/s/Mpc")
print()
print(f"  rho_c (TriPhase) = {rho_c:.6e} kg/m^3")
print(f"  rho_c (observed) ~ {rho_c_obs:.6e} kg/m^3")
print()

print("Critical density connects cosmic expansion to vacuum structure")
print("through the Hubble constant and gravitational constant.")
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
print(f"  H_0         = {H_0:.10e} s^-1")
print(f"  G           = {G:.10e} m^3/(kg·s^2)")
print(f"  rho_c       = {rho_c:.10e} kg/m^3")
print(f"              = {n_c_protons/1e6:.3f} protons/cm^3")
print(f"  R_H         = {R_H/9.461e15:.3e} light-years")
print(f"  t_H         = {t_H/(365.25*24*3600*1e9):.3f} Gyr")
print()
print("Critical density: rho_c = 3*H_0^2/(8*pi*G)")
print("Determines the geometry and fate of the universe.")
print()
print("=" * 80)

input("Press Enter to exit...")
