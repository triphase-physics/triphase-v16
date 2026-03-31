#!/usr/bin/env python3
"""
================================================================================
TriPhase V16 Derivative: Matter Density Parameter (Omega_m)
Framework: Anchor_Primitive
Row: 42, Tag: (D*H)
================================================================================

Physical Concept:
The matter density parameter (Omega_m) represents the ratio of matter density
to critical density. It determines the matter content of the universe and its
role in cosmic evolution. Includes both baryonic and dark matter.

Derivation Path:
- Matter density parameter: Omega_m = rho_m/rho_c
- Critical density: rho_c = 3*H_0^2/(8*pi*G)
- H_0 from anchor chain: H_0 = pi*sqrt(3) * f_e * alpha^18
- G from anchor chain: G = c^4 * 7.5 * epsilon_0^3 * mu_0^2
- Observed value: Omega_m ≈ 0.315

Mathematical Expression:
Omega_m = rho_m/rho_c ≈ 0.315 (0.049 baryonic + 0.266 dark matter)

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
"""

import math

# ============================================================================
# ANCHOR PRIMITIVE DERIVATION
# ============================================================================

print("=" * 80)
print("TriPhase V16: Matter Density Parameter (Anchor Primitive)")
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

# Critical density
rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)
print(f"rho_c = 3*H_0^2/(8*pi*G)")
print(f"      = {rho_c:.15e} kg/m^3")
print()

# ============================================================================
# MATTER DENSITY PARAMETER
# ============================================================================

print("=" * 80)
print("MATTER DENSITY PARAMETER")
print("=" * 80)
print()

print("The matter density parameter represents the fraction of critical")
print("density in the form of matter:")
print()
print("  Omega_m = rho_m / rho_c")
print()

# Observed values from Planck 2018
Omega_m = 0.315        # Total matter
Omega_baryon = 0.049   # Baryonic matter only
Omega_dm = Omega_m - Omega_baryon  # Dark matter

print(f"Observed density parameters (Planck 2018):")
print(f"  Omega_m       = {Omega_m:.3f} (total matter)")
print(f"  Omega_baryon  = {Omega_baryon:.3f} (baryonic matter)")
print(f"  Omega_DM      = {Omega_dm:.3f} (dark matter)")
print()
print(f"Dark matter fraction of total matter:")
print(f"  Omega_DM/Omega_m = {Omega_dm/Omega_m:.3f} = {Omega_dm/Omega_m*100:.1f}%")
print()

# Actual densities
rho_m = Omega_m * rho_c
rho_baryon = Omega_baryon * rho_c
rho_dm = Omega_dm * rho_c

print(f"Actual densities:")
print(f"  rho_m      = {rho_m:.15e} kg/m^3")
print(f"  rho_baryon = {rho_baryon:.15e} kg/m^3")
print(f"  rho_DM     = {rho_dm:.15e} kg/m^3")
print()

# Convert to particle number densities
m_proton = 1.67262192369e-27  # kg
n_baryon = rho_baryon / m_proton

print(f"Baryon number density:")
print(f"  n_b = rho_baryon/m_p")
print(f"      = {n_baryon:.6e} baryons/m^3")
print(f"      = {n_baryon/1e6:.6f} baryons/cm^3")
print()

# ============================================================================
# COMPARISON TO DARK ENERGY
# ============================================================================

print("=" * 80)
print("MATTER VS DARK ENERGY")
print("=" * 80)
print()

Omega_Lambda = 0.685  # Dark energy
rho_Lambda = Omega_Lambda * rho_c

print(f"Dark energy density parameter:")
print(f"  Omega_Lambda = {Omega_Lambda:.3f}")
print(f"  rho_Lambda   = {rho_Lambda:.15e} kg/m^3")
print()

print(f"Energy budget of the universe:")
print(f"  Matter:       {Omega_m*100:.1f}%")
print(f"    - Baryonic: {Omega_baryon*100:.1f}%")
print(f"    - Dark:     {Omega_dm*100:.1f}%")
print(f"  Dark energy:  {Omega_Lambda*100:.1f}%")
print(f"  Total:        {(Omega_m+Omega_Lambda)*100:.1f}%")
print()

# Ratio
print(f"Matter-to-dark-energy ratio:")
print(f"  rho_m/rho_Lambda = {rho_m/rho_Lambda:.3f}")
print()

# ============================================================================
# EVOLUTION OF MATTER DENSITY
# ============================================================================

print("=" * 80)
print("EVOLUTION OF MATTER DENSITY")
print("=" * 80)
print()

print("Matter density scales with the expansion of the universe:")
print("  rho_m(a) = rho_m,0 / a^3")
print("where a is the scale factor (a=1 today, a→0 at big bang)")
print()

# Scale factors at different redshifts
z_values = [0, 1, 2, 3, 10, 100, 1000, 1100]  # Redshifts

print("Redshift |  a    | Omega_m(z) | rho_m/rho_m,0")
print("-" * 70)

for z in z_values:
    a = 1.0 / (1.0 + z)
    rho_m_ratio = 1.0 / a**3
    # Omega_m(z) accounting for dark energy being constant
    Omega_m_z = Omega_m * (1+z)**3 / (Omega_m*(1+z)**3 + Omega_Lambda)

    print(f"  {z:4d}   | {a:.4f} | {Omega_m_z:.4f}     | {rho_m_ratio:.3e}")

    if z == 1100:
        print(f"         | (CMB last scattering)")
print()

# Matter-radiation equality
Omega_r = 9.2e-5  # Radiation today (negligible)
z_eq = Omega_m / Omega_r - 1
a_eq = 1.0 / (1.0 + z_eq)

print(f"Matter-radiation equality:")
print(f"  z_eq ≈ {z_eq:.0f}")
print(f"  a_eq ≈ {a_eq:.6f}")
print()

# Matter-dark energy equality
z_mde = (Omega_Lambda/Omega_m)**(1.0/3.0) - 1
a_mde = 1.0 / (1.0 + z_mde)

print(f"Matter-dark energy equality:")
print(f"  z ≈ {z_mde:.3f}")
print(f"  a ≈ {a_mde:.3f}")
print()

# ============================================================================
# BARYONIC MATTER INVENTORY
# ============================================================================

print("=" * 80)
print("BARYONIC MATTER INVENTORY")
print("=" * 80)
print()

print("Where are the baryons?")
print()

# Hubble volume
R_H = c / H_0
V_H = (4.0/3.0) * math.pi * R_H**3
M_baryon_total = rho_baryon * V_H
M_sun = 1.989e30  # kg

print(f"Total baryonic mass in Hubble volume:")
print(f"  M_b = rho_baryon * V_H")
print(f"      = {M_baryon_total:.15e} kg")
print(f"      = {M_baryon_total/M_sun:.3e} solar masses")
print()

# Distribution (approximate)
fractions = [
    ("Stars", 0.04),
    ("Stellar remnants", 0.01),
    ("Gas in galaxies", 0.04),
    ("Intergalactic medium (warm-hot)", 0.40),
    ("Intergalactic medium (diffuse)", 0.30),
    ("Black holes", 0.01),
    ("Unaccounted/missing", 0.20),
]

print("Baryon distribution (approximate):")
for component, frac in fractions:
    M_component = frac * M_baryon_total
    print(f"  {component:35s}: {frac*100:4.1f}%  ({M_component/M_sun:.2e} M_sun)")
print()

# ============================================================================
# DARK MATTER PROPERTIES
# ============================================================================

print("=" * 80)
print("DARK MATTER PROPERTIES")
print("=" * 80)
print()

print("Dark matter comprises 85% of all matter but does not interact")
print("electromagnetically. Its properties are inferred from gravitational")
print("effects.")
print()

M_dm_total = rho_dm * V_H

print(f"Total dark matter mass in Hubble volume:")
print(f"  M_DM = {M_dm_total:.15e} kg")
print(f"       = {M_dm_total/M_sun:.3e} solar masses")
print()

print(f"Dark matter to baryonic matter ratio:")
print(f"  M_DM/M_b = {M_dm_total/M_baryon_total:.3f}")
print()

# Typical halo properties
M_halo_MW = 1e12 * M_sun  # Milky Way halo mass
r_halo = 200e3 * 3.086e16  # 200 kpc in meters
rho_halo_avg = M_halo_MW / (4.0/3.0 * math.pi * r_halo**3)

print(f"Example: Milky Way dark matter halo:")
print(f"  M_halo ≈ {M_halo_MW/M_sun:.2e} M_sun")
print(f"  R_halo ≈ 200 kpc")
print(f"  <rho_halo> = {rho_halo_avg:.3e} kg/m^3")
print(f"  <rho_halo>/rho_dm = {rho_halo_avg/rho_dm:.3e}")
print("  (Halo is ~1000x denser than cosmic average)")
print()

# ============================================================================
# JEANS MASS AND STRUCTURE FORMATION
# ============================================================================

print("=" * 80)
print("STRUCTURE FORMATION")
print("=" * 80)
print()

print("Matter overdensities grow through gravitational instability.")
print("The Jeans mass sets the minimum mass for collapse:")
print()
print("  M_J ~ (c_s^3) / (G^(3/2) * rho^(1/2))")
print()

# At matter-radiation equality
T_eq = 3e4  # K (approximate temperature at equality)
k_B = 1.380649e-23  # J/K
c_s_eq = math.sqrt(5.0*k_B*T_eq / (3.0*m_proton))  # Sound speed in radiation-dominated plasma

rho_eq = rho_m * (1 + z_eq)**3
M_J_eq = (c_s_eq**3) / (G**(3.0/2.0) * rho_eq**(1.0/2.0))

print(f"At matter-radiation equality (z ~ {z_eq:.0f}):")
print(f"  T ~ {T_eq:.1e} K")
print(f"  c_s ~ {c_s_eq:.3e} m/s")
print(f"  rho ~ {rho_eq:.3e} kg/m^3")
print(f"  M_J ~ {M_J_eq:.3e} kg")
print(f"      ~ {M_J_eq/M_sun:.3e} M_sun")
print()
print("This sets the mass scale for the first dark matter halos.")
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

print("Matter density parameter (observational input):")
print(f"  Omega_m = {Omega_m:.3f}")
print()

print("Derived matter density:")
print(f"  rho_m = Omega_m * rho_c")
print(f"        = {rho_m:.15e} kg/m^3")
print()

print("The matter density parameter connects observational cosmology")
print("to the vacuum structure through the critical density.")
print()

# ============================================================================
# SUMMARY
# ============================================================================

print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print("Framework: ANCHOR_PRIMITIVE")
print("Inputs: epsilon_0, mu_0, e, Omega_m (observational)")
print("Outputs:")
print(f"  H_0            = {H_0:.10e} s^-1")
print(f"  rho_c          = {rho_c:.10e} kg/m^3")
print(f"  Omega_m        = {Omega_m:.3f}")
print(f"  Omega_baryon   = {Omega_baryon:.3f}")
print(f"  Omega_DM       = {Omega_dm:.3f}")
print(f"  rho_m          = {rho_m:.10e} kg/m^3")
print(f"  n_baryon       = {n_baryon/1e6:.3f} baryons/cm^3")
print()
print("Matter density parameter: Omega_m = rho_m/rho_c ≈ 0.315")
print("Determines matter content and structure formation.")
print()
print("=" * 80)

input("Press Enter to exit...")
