"""
critical_density_Anchor_Condensed.py

TriPhase V16 - Critical Density
Row 41 - Tag: (D) DERIVED

Derives critical density rho_c = 3*H_0^2/(8*pi*G)
Both H_0 and G fully derived from epsilon_0, mu_0 via anchor chain.

Critical density is the density needed for flat universe (k=0).
Shows various unit conversions.

All derived from epsilon_0 and mu_0 via the anchor chain.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# Anchor chain
epsilon_0 = 8.8541878128e-12  # F/m
mu_0      = 1.25663706212e-6   # H/m
e         = 1.602176634e-19    # C (exact SI)

c     = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0   = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha = 1.0 / alpha_inv
hbar  = Z_0 * e**2 / (4.0 * math.pi * alpha)
h     = 2.0 * math.pi * hbar
G     = c**4 * 7.5 * epsilon_0**3 * mu_0**2
m_e   = hbar * alpha / (c * 2.8179403262e-15)
f_e   = m_e * c**2 / hbar
mp_me = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p   = m_e * mp_me
H_0   = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r  = c**4 / (8.0 * math.pi * G)

print("=" * 70)
print("TriPhase V16 - Critical Density")
print("Row 41 - DERIVED from epsilon_0, mu_0")
print("=" * 70)
print()

# Critical density
rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)

print("CRITICAL DENSITY:")
print(f"  rho_c = 3*H_0^2/(8*pi*G)")
print(f"        = {rho_c:.15e} kg/m^3")
print()
print("  This is the density required for a flat universe (k=0)")
print("  Omega_total = rho_total/rho_c")
print("    Omega < 1 → open (hyperbolic)")
print("    Omega = 1 → flat (Euclidean)")
print("    Omega > 1 → closed (spherical)")
print()

# Unit conversions
print("UNIT CONVERSIONS:")
print()

# Protons per cubic meter
n_protons = rho_c / m_p
print(f"  Proton number density:")
print(f"    n_p = rho_c/m_p = {n_protons:.6e} protons/m^3")
print(f"        = {n_protons*1e-6:.6f} protons/cm^3")
print()

# GeV per cubic meter
u_GeV = rho_c * c**2 / (e * 1e9)
print(f"  Energy density:")
print(f"    u = rho_c*c^2 = {rho_c*c**2:.6e} J/m^3")
print(f"      = {u_GeV:.6e} GeV/m^3")
print()

# eV per cubic centimeter
u_eV_cm3 = rho_c * c**2 / e * 1e-6
print(f"    u = {u_eV_cm3:.6f} eV/cm^3")
print()

# grams per cubic centimeter
rho_c_cgs = rho_c * 1e-3
print(f"  CGS units:")
print(f"    rho_c = {rho_c_cgs:.6e} g/cm^3")
print()

# Compare to familiar densities
print("COMPARISON TO FAMILIAR DENSITIES:")
densities = [
    ("Best lab vacuum", 1e-17),
    ("Interstellar medium", 1e-21),
    ("Critical density", rho_c),
    ("Intergalactic medium", 1e-27),
    ("Air (sea level)", 1.225),
    ("Water", 1000.0),
]

for name, rho in sorted(densities, key=lambda x: x[1]):
    ratio = rho / rho_c
    print(f"  {name:25s}: {rho:.6e} kg/m^3 (rho/rho_c = {ratio:.2e})")
print()

# Hubble parameter and critical density
print("HUBBLE PARAMETER:")
print(f"  H_0 = {H_0:.15e} Hz")
print(f"      = {H_0*1e3:.6f} milliHz")
print(f"  H_0/(2*pi) = {H_0/(2.0*math.pi):.6e} Hz")
print()

# Hubble time and length
t_H = 1.0 / H_0
d_H = c / H_0

print("HUBBLE SCALES:")
print(f"  Hubble time: t_H = 1/H_0")
print(f"             = {t_H:.6e} s")
print(f"             = {t_H/(365.25*24*3600):.3e} years")
print(f"             = {t_H/(365.25*24*3600*1e9):.3f} Gyr")
print()
print(f"  Hubble length: d_H = c/H_0")
print(f"               = {d_H:.6e} m")
print(f"               = {d_H/9.461e15:.3e} ly")
print(f"               = {d_H/(9.461e15*1e9):.3f} Gly")
print()

# Friedmann equation
print("FRIEDMANN EQUATION:")
print("  H^2 = (8*pi*G/3)*rho - k*c^2/a^2 + Lambda*c^2/3")
print()
print("  For flat universe (k=0, Lambda=0):")
print("    H^2 = (8*pi*G/3)*rho")
print("    rho_c = 3*H^2/(8*pi*G)")
print()

# Current energy budget
Omega_m = 0.315   # Planck 2018
Omega_b = 0.049   # baryons
Omega_DM = Omega_m - Omega_b  # dark matter
Omega_DE = 0.685  # dark energy
Omega_total = Omega_m + Omega_DE

print("CURRENT ENERGY BUDGET (Planck 2018):")
print(f"  Omega_baryons     = {Omega_b:.3f} ({100.0*Omega_b:.1f}%)")
print(f"  Omega_dark_matter = {Omega_DM:.3f} ({100.0*Omega_DM:.1f}%)")
print(f"  Omega_matter      = {Omega_m:.3f} ({100.0*Omega_m:.1f}%)")
print(f"  Omega_dark_energy = {Omega_DE:.3f} ({100.0*Omega_DE:.1f}%)")
print(f"  Omega_total       = {Omega_total:.3f}")
print()

# Actual densities
rho_b = Omega_b * rho_c
rho_DM = Omega_DM * rho_c
rho_m = Omega_m * rho_c
rho_DE = Omega_DE * rho_c

print("ACTUAL DENSITIES:")
print(f"  rho_baryons     = {rho_b:.6e} kg/m^3")
print(f"                  = {rho_b/m_p:.6e} protons/m^3")
print(f"                  = {rho_b/m_p*1e-6:.6f} protons/cm^3")
print()
print(f"  rho_dark_matter = {rho_DM:.6e} kg/m^3")
print()
print(f"  rho_matter      = {rho_m:.6e} kg/m^3")
print()
print(f"  rho_dark_energy = {rho_DE:.6e} kg/m^3")
print()

# Compare to vacuum frame rigidity
P_c = rho_c * c**2

print("COMPARISON TO VACUUM FRAME RIGIDITY:")
print(f"  P_c = rho_c*c^2 = {P_c:.6e} Pa")
print(f"  VF_r = {VF_r:.6e} Pa")
print(f"  P_c/VF_r = {P_c/VF_r:.6e}")
print()
print("  Critical density pressure is ~52 orders below vacuum rigidity")
print()

print("ANCHOR CHAIN DERIVATION:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print(f"  c         = {c:.10e} m/s")
print(f"  alpha     = {alpha:.15e}")
print(f"  hbar      = {hbar:.15e} J·s")
print(f"  m_e       = {m_e:.15e} kg")
print(f"  f_e       = m_e*c^2/hbar = {f_e:.15e} Hz")
print(f"  m_p       = {m_p:.15e} kg")
print()
print(f"  G         = c^4 * 7.5 * epsilon_0^3 * mu_0^2")
print(f"            = {G:.14e} m^3/(kg·s^2)")
print()
print(f"  H_0       = pi*sqrt(3)*f_e*alpha^18")
print(f"            = {H_0:.15e} Hz")
print()
print(f"  rho_c     = 3*H_0^2/(8*pi*G)")
print(f"            = {rho_c:.15e} kg/m^3")
print()

# CODATA comparison
H_0_Planck = 67.4e3 / (3.086e22)  # 67.4 km/s/Mpc in Hz
G_CODATA = 6.67430e-11
rho_c_Planck = 3.0 * H_0_Planck**2 / (8.0 * math.pi * G_CODATA)

print("PLANCK 2018 COMPARISON (calibration only):")
print(f"  H_0 (Planck) = 67.4 km/s/Mpc = {H_0_Planck:.6e} Hz")
print(f"  H_0 (derived) = {H_0:.6e} Hz")
print(f"  Error: {100.0*abs(H_0 - H_0_Planck)/H_0_Planck:.3f}%")
print()
print(f"  rho_c (Planck) = {rho_c_Planck:.6e} kg/m^3")
print(f"  rho_c (derived) = {rho_c:.6e} kg/m^3")
print(f"  Error: {100.0*abs(rho_c - rho_c_Planck)/rho_c_Planck:.3f}%")
print()

print("=" * 70)
print("Critical density derived from H_0 and G, both from epsilon_0, mu_0")
print("The cosmic density that determines spacetime curvature")
print("=" * 70)

input("Press Enter to exit...")
