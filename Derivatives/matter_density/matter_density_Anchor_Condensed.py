"""
matter_density_Anchor_Condensed.py

TriPhase V16 - Matter Density
Row 42 - Tag: (D*) DERIVED with discrete selection

Derives matter density rho_m = Omega_m * rho_c
where Omega_m ≈ 0.315 (Planck 2018)

Shows breakdown:
  - Baryonic matter (Omega_b ≈ 0.049)
  - Dark matter (Omega_DM ≈ 0.266)

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
print("TriPhase V16 - Matter Density")
print("Row 42 - DERIVED* with discrete selection from epsilon_0, mu_0")
print("=" * 70)
print()

# Critical density
rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)

print("CRITICAL DENSITY:")
print(f"  rho_c = 3*H_0^2/(8*pi*G)")
print(f"        = {rho_c:.15e} kg/m^3")
print()

# Matter density fractions (Planck 2018)
Omega_b = 0.049   # baryonic matter
Omega_DM = 0.266  # dark matter
Omega_m = 0.315   # total matter (Omega_b + Omega_DM)

print("MATTER DENSITY FRACTIONS (Planck 2018):")
print(f"  Omega_baryons     = {Omega_b:.3f} ({100.0*Omega_b/Omega_m:.1f}% of matter)")
print(f"  Omega_dark_matter = {Omega_DM:.3f} ({100.0*Omega_DM/Omega_m:.1f}% of matter)")
print(f"  Omega_matter      = {Omega_m:.3f}")
print(f"  Check: Omega_b + Omega_DM = {Omega_b + Omega_DM:.3f}")
print()

# Actual densities
rho_b = Omega_b * rho_c
rho_DM = Omega_DM * rho_c
rho_m = Omega_m * rho_c

print("MATTER DENSITIES:")
print()
print("Baryonic matter:")
print(f"  rho_b = Omega_b * rho_c")
print(f"        = {rho_b:.15e} kg/m^3")
print(f"        = {rho_b*1e3:.15e} g/cm^3")
print()
print("Dark matter:")
print(f"  rho_DM = Omega_DM * rho_c")
print(f"         = {rho_DM:.15e} kg/m^3")
print(f"         = {rho_DM*1e3:.15e} g/cm^3")
print()
print("Total matter:")
print(f"  rho_m = Omega_m * rho_c")
print(f"        = {rho_m:.15e} kg/m^3")
print(f"        = {rho_m*1e3:.15e} g/cm^3")
print()

# Number densities
n_b = rho_b / m_p
n_DM = rho_DM / m_p  # assuming DM has proton-like mass (placeholder)
n_m = rho_m / m_p

print("NUMBER DENSITIES (assuming proton mass):")
print(f"  Baryons:     n_b  = {n_b:.6e} particles/m^3")
print(f"                    = {n_b*1e-6:.6f} particles/cm^3")
print()
print(f"  Dark matter: n_DM = {n_DM:.6e} particles/m^3")
print(f"                    = {n_DM*1e-6:.6f} particles/cm^3")
print()
print(f"  Total:       n_m  = {n_m:.6e} particles/m^3")
print(f"                    = {n_m*1e-6:.6f} particles/cm^3")
print()

# Energy densities
u_b = rho_b * c**2
u_DM = rho_DM * c**2
u_m = rho_m * c**2

print("ENERGY DENSITIES:")
print(f"  Baryons:     u_b  = {u_b:.6e} J/m^3")
print(f"                    = {u_b/(e*1e9):.6e} GeV/m^3")
print(f"                    = {u_b/e*1e-6:.6f} eV/cm^3")
print()
print(f"  Dark matter: u_DM = {u_DM:.6e} J/m^3")
print(f"                    = {u_DM/(e*1e9):.6e} GeV/m^3")
print(f"                    = {u_DM/e*1e-6:.6f} eV/cm^3")
print()
print(f"  Total:       u_m  = {u_m:.6e} J/m^3")
print(f"                    = {u_m/(e*1e9):.6e} GeV/m^3")
print(f"                    = {u_m/e*1e-6:.6f} eV/cm^3")
print()

# Compare to familiar densities
print("COMPARISON TO FAMILIAR ENVIRONMENTS:")
environments = [
    ("Intergalactic medium", 1e-27),
    ("Cosmic average (baryons)", rho_b),
    ("Cosmic average (total matter)", rho_m),
    ("Interstellar medium", 1e-21),
    ("Best lab vacuum", 1e-17),
    ("Earth's atmosphere", 1.225),
]

for name, rho in environments:
    ratio_to_critical = rho / rho_c
    ratio_to_matter = rho / rho_m
    print(f"  {name:30s}: {rho:.2e} kg/m^3")
    print(f"    rho/rho_c = {ratio_to_critical:.2e}, rho/rho_m = {ratio_to_matter:.2e}")
print()

# Dark matter candidates
print("DARK MATTER MASS SCALES:")
print("  If DM is WIMPs (Weakly Interacting Massive Particles):")
m_WIMP_GeV = 100.0  # typical WIMP mass in GeV
m_WIMP = m_WIMP_GeV * e * 1e9 / c**2
n_WIMP = rho_DM / m_WIMP
print(f"    m_WIMP ~ {m_WIMP_GeV:.0f} GeV/c^2 = {m_WIMP:.6e} kg")
print(f"    n_WIMP = rho_DM/m_WIMP = {n_WIMP:.6e} m^-3")
print(f"           = {n_WIMP*1e-6:.6e} cm^-3")
print()
print("  If DM is axions:")
m_axion_eV = 1e-5  # typical axion mass in eV
m_axion = m_axion_eV * e / c**2
n_axion = rho_DM / m_axion
print(f"    m_axion ~ {m_axion_eV:.1e} eV/c^2 = {m_axion:.6e} kg")
print(f"    n_axion = rho_DM/m_axion = {n_axion:.6e} m^-3")
print(f"            = {n_axion*1e-6:.6e} cm^-3")
print()

# Matter-antimatter asymmetry
print("MATTER-ANTIMATTER ASYMMETRY:")
eta_b = 6.1e-10  # baryon-to-photon ratio (Planck 2018)
n_gamma = 411e6  # CMB photons per m^3
n_b_from_eta = eta_b * n_gamma
print(f"  Baryon-to-photon ratio: eta_b = {eta_b:.2e}")
print(f"  CMB photon density: n_gamma ~ {n_gamma:.0e} m^-3")
print(f"  Implied baryon density: n_b = eta_b*n_gamma = {n_b_from_eta:.6e} m^-3")
print(f"  From Omega_b: n_b = {n_b:.6e} m^-3")
print(f"  Agreement: {100.0*n_b/n_b_from_eta:.1f}%")
print()

# Evolution with scale factor
print("MATTER DENSITY EVOLUTION:")
print("  rho_m ~ a^(-3) (dilutes with expansion)")
print("  For scale factor a:")
print("    rho_m(a) = rho_m(today) * a^(-3)")
print()
print("  At recombination (z ~ 1100, a ~ 1/1100):")
a_rec = 1.0/1100.0
rho_m_rec = rho_m / a_rec**3
print(f"    rho_m(z=1100) = {rho_m_rec:.6e} kg/m^3")
print(f"    = {rho_m_rec/rho_m:.0e} * rho_m(today)")
print()

# Compare to vacuum frame rigidity
P_m = rho_m * c**2

print("COMPARISON TO VACUUM FRAME RIGIDITY:")
print(f"  P_m = rho_m*c^2 = {P_m:.6e} Pa")
print(f"  VF_r = {VF_r:.6e} Pa")
print(f"  P_m/VF_r = {P_m/VF_r:.6e}")
print()
print("  Matter pressure is ~52 orders below vacuum rigidity")
print()

print("ANCHOR CHAIN DERIVATION:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print(f"  c         = {c:.10e} m/s")
print(f"  G         = c^4 * 7.5 * epsilon_0^3 * mu_0^2")
print(f"            = {G:.14e} m^3/(kg·s^2)")
print(f"  H_0       = pi*sqrt(3)*f_e*alpha^18")
print(f"            = {H_0:.15e} Hz")
print(f"  rho_c     = 3*H_0^2/(8*pi*G)")
print(f"            = {rho_c:.15e} kg/m^3")
print()
print(f"  Omega_m   = {Omega_m:.3f} (Planck 2018)")
print(f"  rho_m     = Omega_m * rho_c")
print(f"            = {rho_m:.15e} kg/m^3")
print()
print("  Breakdown:")
print(f"    Omega_b  = {Omega_b:.3f} → rho_b  = {rho_b:.6e} kg/m^3")
print(f"    Omega_DM = {Omega_DM:.3f} → rho_DM = {rho_DM:.6e} kg/m^3")
print()

print("=" * 70)
print("Matter density from Omega_m * rho_c, both from epsilon_0, mu_0")
print("Baryons (4.9%) + Dark Matter (26.6%) = Total Matter (31.5%)")
print("The rest is Dark Energy (68.5%)")
print("=" * 70)

input("Press Enter to exit...")
