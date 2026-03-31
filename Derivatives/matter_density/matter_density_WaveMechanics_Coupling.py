# -*- coding: utf-8 -*-
"""
============================================================
TRIPHASE PHYSICS FRAMEWORK V16: A Foundation in Wave Mechanics
============================================================
MATTER DENSITY (Ω_m = ρ_m/ρ_c ~ 0.31)
Framework: WaveMechanics_Coupling
Tag: (C)

Matter fraction of total coupling budget.
ρ_m = Ω_m × ρ_c, where ρ_c derives from H₀ and G.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
============================================================
"""

import numpy as np
from scipy import constants
import sys
import io
if sys.stdout.encoding != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')


# ============================================================
# VACUUM CONSTANTS (PRIMARY INPUTS)
# ============================================================
epsilon_0 = 8.8541878128e-12  # F/m (vacuum permittivity)
mu_0 = 1.25663706212e-6       # H/m (vacuum permeability)

# Derived vacuum properties
c = 1 / np.sqrt(epsilon_0 * mu_0)  # Speed of light
Z_0 = np.sqrt(mu_0 / epsilon_0)    # Vacuum impedance

# Alpha (fine structure constant)
m = 17
node = 8 * m + 1  # = 137
correction = np.log(node) / node
alpha_inv = node + correction
alpha = 1 / alpha_inv

# Gravitational constant (from vacuum)
G = c**4 * 7.5 * epsilon_0**3 * mu_0**2

# Hubble constant (from vacuum)
m_e = constants.m_e
hbar = constants.hbar
f_e = m_e * c**2 / hbar  # Electron frequency
H_0 = np.pi * np.sqrt(3) * f_e * alpha**18  # rad/s
H_0_km_s_Mpc = H_0 / 1000 * 3.086e22  # Convert to km/s/Mpc

# Critical density
rho_c = 3 * H_0**2 / (8 * np.pi * G)

# Matter density (Omega_m from Planck 2018)
Omega_m = 0.309
rho_m = Omega_m * rho_c

# Reference values
G_CODATA = constants.G           # 6.67430e-11 m³/(kg·s²)
H_0_obs = 67.4 * 1000 / 3.086e22  # Planck 2018: 67.4 km/s/Mpc → rad/s
rho_c_obs = 3 * H_0_obs**2 / (8 * np.pi * G_CODATA)
rho_m_obs = Omega_m * rho_c_obs


# ============================================================
# THE THREE AXIOMS
# ============================================================
AXIOMS = """
| AXIOM                                           | WAVE PROPERTY        |
|-------------------------------------------------|----------------------|
| 1. Energy accumulation drives frequency adjust  | FREQUENCY/WAVELENGTH |
| 2. Three-phase oscillation is fundamental mode  | PHASE (x3)           |
| 3. Vector Frame Imperivity provides equilibrium | AMPLITUDE            |

COUPLING HIERARCHY:
- H₀ = π√3 f_e α¹⁸ (expansion rate)
- G = c⁴ × 7.5 × ε₀³ × μ₀² (gravitational coupling)
- ρ_c = 3H₀²/(8πG) (critical density)
- Ω_m = 0.309 (matter fraction, observational input)
- ρ_m = Ω_m × ρ_c (matter density)
"""


# ============================================================
# WAVEMECHANICS_COUPLING FRAMEWORK
# ============================================================
MECHANISM = """
MATTER DENSITY: FRACTION OF CRITICAL DENSITY

From Friedmann equation:
  Ω_m = ρ_m / ρ_c

where:
  ρ_m = total matter density (dark matter + baryonic)
  ρ_c = critical density = 3H₀²/(8πG)
  Ω_m ≈ 0.309 (Planck 2018)

COUPLING INTERPRETATION:
  - Ω_m = matter's share of cosmic coupling budget
  - Dark matter: Ω_DM ≈ 0.259 (84% of matter)
  - Baryonic matter: Ω_b ≈ 0.049 (16% of matter)
  - Matter fraction decreases as universe expands
  - Dark energy (Ω_DE ≈ 0.691) now dominates

TRIPHASE DERIVATION (all from ε₀, μ₀):
  1. α from m=17, node=137 structure
  2. f_e = m_e c²/ℏ (electron frequency)
  3. H₀ = π√3 f_e α¹⁸ (Hubble constant)
  4. G = c⁴ × 7.5 × ε₀³ × μ₀² (gravity)
  5. ρ_c = 3H₀²/(8πG) (critical density)
  6. Ω_m = 0.309 (observational input)
  7. ρ_m = Ω_m × ρ_c (matter density)

HISTORICAL EVOLUTION:
  Early universe (z >> 1):
    Ω_m(z) ≈ 1 (matter dominated)

  Today (z = 0):
    Ω_m = 0.309 (dark energy dominant)

  Future (z → -1):
    Ω_m → 0 (exponential expansion)

MATTER-RADIATION EQUALITY:
  z_eq ≈ 3400 (when Ω_m = Ω_r)
  Before: radiation dominated
  After: matter dominated
  Today: dark energy dominated

MATTER-DARK ENERGY EQUALITY:
  z_DE ≈ 0.33 (when Ω_m = Ω_DE)
  Occurred ~4.6 Gyr ago
  Transition to accelerated expansion
"""


# ============================================================
# CALCULATION
# ============================================================
def calculate():
    print("=" * 70)
    print("TRIPHASE PHYSICS FRAMEWORK V16: A Foundation in Wave Mechanics")
    print("=" * 70)
    print("MATTER DENSITY | Framework: WaveMechanics_Coupling")
    print("Tag: (C) - Calibrated (Ω_m observational input)")
    print("=" * 70)

    print("\nTHE THREE AXIOMS")
    print("-" * 70)
    print(AXIOMS)

    print("\nTHE MECHANISM: MATTER FRACTION")
    print("-" * 70)
    print(MECHANISM)

    print("CALCULATION")
    print("-" * 70)

    # Verify vacuum constants
    print("\nSTEP 1: Vacuum properties")
    print(f"  ε₀ = {epsilon_0:.6e} F/m")
    print(f"  μ₀ = {mu_0:.6e} H/m")
    print(f"  c = 1/√(ε₀μ₀) = {c:.6e} m/s")
    print(f"  α = 1/{alpha_inv:.6f}")

    # Hubble constant
    error_H0 = (H_0_km_s_Mpc - 67.4) / 67.4 * 100
    print("\nSTEP 2: Hubble constant")
    print(f"  f_e = m_e c²/ℏ = {f_e:.6e} Hz")
    print(f"  H₀ = π√3 f_e α¹⁸ = {H_0:.6e} rad/s")
    print(f"  = {H_0_km_s_Mpc:.2f} km/s/Mpc")
    print(f"  Planck 2018: 67.4 km/s/Mpc")
    print(f"  ERROR: {error_H0:+.2f}%")

    # Gravitational constant
    error_G = (G - G_CODATA) / G_CODATA * 100
    print("\nSTEP 3: Gravitational constant")
    print(f"  G = c⁴ × 7.5 × ε₀³ × μ₀² = {G:.6e} m³/(kg·s²)")
    print(f"  CODATA: {G_CODATA:.6e} m³/(kg·s²)")
    print(f"  ERROR: {error_G:+.2f}%")

    # Critical density
    error_rho_c = (rho_c - rho_c_obs) / rho_c_obs * 100
    print("\nSTEP 4: Critical density")
    print(f"  ρ_c = 3H₀²/(8πG) = {rho_c:.6e} kg/m³")
    print(f"  Observed: {rho_c_obs:.6e} kg/m³")
    print(f"  ERROR: {error_rho_c:+.2f}%")

    # Matter density
    error_rho_m = (rho_m - rho_m_obs) / rho_m_obs * 100
    print("\nSTEP 5: Matter density")
    print(f"  Ω_m = {Omega_m:.3f} (Planck 2018)")
    print(f"  ρ_m = Ω_m × ρ_c = {rho_m:.6e} kg/m³")
    print(f"  Observed: {rho_m_obs:.6e} kg/m³")
    print(f"  ERROR: {error_rho_m:+.2f}%")

    # Unit conversions
    rho_m_g_cm3 = rho_m * 1e-3  # kg/m³ → g/cm³
    m_p = constants.m_p
    n_p_equivalent = rho_m / m_p
    print("\nSTEP 6: Unit conversions")
    print(f"  ρ_m = {rho_m_g_cm3:.3e} g/cm³")
    print(f"  ≈ {n_p_equivalent:.2f} protons/m³")
    print(f"  (Equivalent to ~2 protons per cubic meter)")

    # Matter components
    print("\nSTEP 7: Matter components (Planck 2018)")
    Omega_DM = 0.259  # Dark matter
    Omega_b = 0.049   # Baryonic matter
    rho_DM = Omega_DM * rho_c
    rho_b = Omega_b * rho_c
    print(f"  Dark matter: Ω_DM = {Omega_DM:.3f}")
    print(f"    ρ_DM = {rho_DM:.6e} kg/m³ ({Omega_DM/Omega_m*100:.1f}% of matter)")
    print(f"  Baryonic: Ω_b = {Omega_b:.3f}")
    print(f"    ρ_b = {rho_b:.6e} kg/m³ ({Omega_b/Omega_m*100:.1f}% of matter)")
    print(f"  Total: Ω_m = {Omega_m:.3f}")

    # Density evolution
    print("\nSTEP 8: Density evolution with redshift")
    print("  ρ_m(z) = ρ_m,0 × (1+z)³  (matter scales as volume)")
    print("  ρ_DE(z) ≈ ρ_DE,0  (dark energy roughly constant)")
    print()
    redshifts = [0, 0.33, 1, 2, 5, 10, 100, 1000, 3400]
    print("  z        Ω_m(z)    Ω_DE(z)   Dominant")
    print("  " + "-" * 50)

    Omega_DE = 0.691
    for z in redshifts:
        Omega_m_z = Omega_m * (1 + z)**3 / (Omega_m * (1 + z)**3 + Omega_DE)
        Omega_DE_z = Omega_DE / (Omega_m * (1 + z)**3 + Omega_DE)
        dominant = "Matter" if Omega_m_z > 0.5 else "Dark Energy"
        if z >= 3000:
            dominant = "Radiation"
        print(f"  {z:5.0f}    {Omega_m_z:.4f}    {Omega_DE_z:.4f}    {dominant}")

    # Matter-dark energy equality
    z_eq_DE = (Omega_DE / Omega_m)**(1/3) - 1
    t_eq_DE = 2 / (3 * H_0 * np.sqrt(Omega_m)) * (1 + z_eq_DE)**(-3/2)
    t_eq_DE_Gyr = t_eq_DE / (1e9 * 365.25 * 24 * 3600)
    age_universe_Gyr = 13.8
    t_ago_Gyr = age_universe_Gyr - t_eq_DE_Gyr
    print("\nSTEP 9: Matter-dark energy equality")
    print(f"  z_eq = {z_eq_DE:.2f} (when Ω_m = Ω_DE)")
    print(f"  Occurred at t ≈ {t_eq_DE_Gyr:.1f} Gyr")
    print(f"  ≈ {t_ago_Gyr:.1f} Gyr ago")
    print(f"  (Transition to accelerated expansion)")

    # Energy densities
    u_m = rho_m * c**2
    u_total = rho_c * c**2
    print("\nSTEP 10: Energy densities")
    print(f"  u_m = ρ_m c² = {u_m:.3e} J/m³")
    print(f"  u_c = ρ_c c² = {u_total:.3e} J/m³")
    print(f"  Matter fraction: u_m/u_c = {u_m/u_total:.3f}")

    print("\n" + "=" * 70)
    print("PHYSICAL MEANING")
    print("=" * 70)
    print("""
MATTER DENSITY: ρ_m = Ω_m × ρ_c

COUPLING INTERPRETATION:
  - Ω_m = matter's share of cosmic coupling budget
  - ρ_c = 3H₀²/(8πG) (critical density)
  - Ω_m = 0.309 → matter is subdominant today
  - Dark energy (Ω_DE = 0.691) now dominates

TRIPHASE DERIVATION (all from ε₀, μ₀):
  1. H₀ = π√3 f_e α¹⁸ (expansion rate)
  2. G = c⁴ × 7.5 × ε₀³ × μ₀² (gravity)
  3. ρ_c = 3H₀²/(8πG) (critical density)
  4. ρ_m = Ω_m × ρ_c (matter density)

MATTER COMPONENTS:
  Dark matter: Ω_DM ≈ 0.259 (84% of matter)
  Baryonic: Ω_b ≈ 0.049 (16% of matter)

HISTORICAL EVOLUTION:
  z > 3400: Radiation dominated (Ω_r > Ω_m)
  3400 > z > 0.33: Matter dominated (Ω_m > Ω_DE)
  z < 0.33: Dark energy dominated (Ω_DE > Ω_m)

MATTER-DARK ENERGY TRANSITION:
  z_eq ≈ 0.33 (when Ω_m = Ω_DE)
  Occurred ~4.6 Gyr ago
  Marks onset of accelerated expansion

WHY THIS MATTERS:
  - Matter density determines structure formation
  - Dark matter provides gravitational scaffolding
  - Baryonic matter forms stars, galaxies, us
  - Decreasing Ω_m → future is dark energy dominated

PHYSICAL SCALE:
  ρ_m ≈ 2.9×10⁻²⁷ kg/m³
  ≈ 2 protons per cubic meter
  Incredibly dilute, yet drives cosmic structure!

TAG: (C) - Calibrated
  Ω_m is observational input (not yet derived from TriPhase)
  ρ_c derives from H₀ and G (both from ε₀, μ₀)
""")
    print("=" * 70)
    print("(c) 2025-2026 MIS Magnetic Innovative Solutions LLC")
    print("=" * 70)

if __name__ == "__main__":
    calculate()
    input("Press Enter to exit...")
