"""
TriPhase V16: Matter Density - QFT Framework
=============================================

QFT INTERPRETATION:
Matter density ρ_m ≈ 2.7×10⁻²⁷ kg/m³ represents the average mass density of all
matter (baryonic + dark) in the universe. This corresponds to Ω_m ≈ 0.315, or
about 31.5% of the critical density.

The matter density breaks down into:
  • Baryonic matter: Ω_b ≈ 0.049 (~5%, protons, neutrons, electrons)
  • Dark matter: Ω_DM ≈ 0.265 (~26%, unknown particle species)

In QFT, baryonic matter is well understood—protons and neutrons are bound states
of quarks and gluons governed by QCD. But dark matter remains a mystery. Candidates
include:
  • WIMPs (Weakly Interacting Massive Particles): χ with m_χ ~ 100 GeV
  • Axions: light bosons a with m_a ~ 10⁻⁵ eV from PQ symmetry breaking
  • Sterile neutrinos: heavy ν_s with m_s ~ keV, mixing with active neutrinos
  • Primordial black holes: M_PBH ~ 10²⁰-10³⁰ kg formed in early universe

Dark matter interacts gravitationally but not electromagnetically (no photon
coupling), making it invisible. Its presence is inferred from:
  • Galaxy rotation curves: v(r) = const, not v ~ 1/√r (Keplerian)
  • Gravitational lensing: light bending exceeds visible mass
  • CMB acoustic peaks: structure formation needs Ω_DM >> Ω_b
  • Large-scale structure: N-body simulations match observations with dark matter

The matter density evolves as ρ_m(a) = ρ_m,0 / a³ (dilutes with expansion),
unlike dark energy which remains constant: ρ_Λ(a) = const. This is why matter
dominated the early universe but dark energy dominates today—the densities
crossed at redshift z_eq ≈ 0.4 (about 4 billion years ago).

TriPhase computes ρ_m = ρ_crit × Ω_m where ρ_crit = 3H_0²/(8πG) and Ω_m ≈ 0.315.
Since H_0 emerges from α¹⁸ × f_e, the matter density connects to electromagnetic
vacuum structure through the Hubble expansion rate.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
TAG: (D) - Direct derivation from cosmological parameters
"""

import math

# ========== ANCHOR CHAIN (VERBATIM) ==========
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

# ========== QFT DERIVATION: MATTER DENSITY ==========
print("=" * 70)
print("  TRIPHASE V16: MATTER DENSITY (QFT Framework)")
print("=" * 70)
print()
print("QFT CONTEXT:")
print("  Matter density ρ_m includes all matter (baryonic + dark) in the")
print("  universe. Observations (Planck 2018) give Ω_m ≈ 0.315, meaning")
print("  matter comprises ~31% of critical density.")
print()
print("  Breakdown:")
print("    • Baryonic matter (p, n, e): Ω_b ≈ 0.049 (~5%)")
print("    • Dark matter (unknown):     Ω_DM ≈ 0.265 (~26%)")
print()
print("  Dark matter doesn't interact electromagnetically but dominates")
print("  gravitational clustering, enabling galaxy formation.")
print()

# Derivation
rho_crit = 3.0 * H_0**2 / (8.0 * math.pi * G)
Omega_m = 0.315  # Total matter density fraction (Planck 2018)
Omega_b = 0.049  # Baryonic matter fraction
Omega_DM = Omega_m - Omega_b  # Dark matter fraction

rho_m = rho_crit * Omega_m
rho_b = rho_crit * Omega_b
rho_DM = rho_crit * Omega_DM

n_baryons = rho_b / m_p  # Baryon number density

print("DERIVATION STEPS:")
print(f"  1. Critical density (from previous derivation):")
print(f"     ρ_crit = 3H_0²/(8πG)")
print(f"     ρ_crit = {rho_crit:.6e} kg/m³")
print()
print(f"  2. Matter density fractions (Planck 2018):")
print(f"     Ω_m (total matter)  = {Omega_m:.3f}")
print(f"     Ω_b (baryons)       = {Omega_b:.3f}")
print(f"     Ω_DM (dark matter)  = {Omega_DM:.3f}")
print()
print(f"  3. Matter density:")
print(f"     ρ_m = ρ_crit × Ω_m")
print(f"     = {rho_crit:.6e} kg/m³ × {Omega_m:.3f}")
print(f"     = {rho_m:.6e} kg/m³")
print()
print(f"  4. Component densities:")
print(f"     ρ_b (baryonic) = {rho_b:.6e} kg/m³")
print(f"     ρ_DM (dark)    = {rho_DM:.6e} kg/m³")
print()
print(f"  5. Baryon number density:")
print(f"     n_b = ρ_b / m_p")
print(f"     = {n_baryons:.3f} baryons/m³")
print(f"     ≈ {n_baryons * 1e-6:.3f} baryons/cm³")
print()

# Calibration
rho_m_Planck = 2.7e-27  # kg/m³ (approximate from Planck)
rho_galaxy = 1e-21  # kg/m³ (typical galaxy density)
rho_solar = 1400  # kg/m³ (Sun's average density)

print("=" * 70)
print("CALIBRATION CHECKPOINT:")
print(f"  TriPhase ρ_m:        {rho_m:.6e} kg/m³")
print(f"  Planck (approx):     {rho_m_Planck:.1e} kg/m³")
print()
print("  Comparison to other densities:")
print(f"    Cosmic matter avg:   {rho_m:.2e} kg/m³  (~0.3 protons/m³)")
print(f"    Typical galaxy:      {rho_galaxy:.1e} kg/m³")
print(f"    Sun (average):       {rho_solar:.1e} kg/m³")
print()
print(f"  ρ_galaxy / ρ_m ≈ {rho_galaxy / rho_m:.2e}  (galaxies are ~10⁶× denser)")
print("=" * 70)
print()

print("QFT INSIGHT:")
print("  Matter density ρ_m determines the universe's large-scale structure.")
print()
print("  1. STRUCTURE FORMATION:")
print("  In the early universe, quantum fluctuations from inflation seeded")
print("  density perturbations δρ/ρ ~ 10⁻⁵. These grew via gravitational")
print("  instability:")
print()
print("    δ̈ + 2Hδ̇ - 4πGρ_m δ = 0")
print()
print("  Matter-dominated era (a ~ t²/³): perturbations grow as δ ~ a ~ t²/³.")
print("  Dark energy era (a ~ e^(Ht)): perturbations freeze (δ = const).")
print()
print("  2. DARK MATTER HALOS:")
print("  N-body simulations show dark matter forms 'halos' with NFW profile:")
print()
print("    ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]")
print()
print("  Baryons fall into these halos, cooling and forming galaxies.")
print("  Without dark matter, structure formation would be too slow!")
print()
print("  3. DARK MATTER CANDIDATES:")
print()
print("  WIMPS (Weakly Interacting Massive Particles):")
print("    • Mass: m_χ ~ 100 GeV (electroweak scale)")
print("    • Interaction: weak force only (no EM, strong)")
print("    • Production: thermal freeze-out in early universe")
print("    • Detection: direct (xenon/germanium detectors), indirect (γ-rays)")
print()
print("  AXIONS:")
print("    • Mass: m_a ~ 10⁻⁵ eV (extremely light)")
print("    • Origin: Peccei-Quinn symmetry breaking (solves strong CP problem)")
print("    • Coupling: a→γγ (two-photon decay)")
print("    • Detection: resonant cavity experiments (ADMX)")
print()
print("  STERILE NEUTRINOS:")
print("    • Mass: m_s ~ keV (between active ν and WIMP)")
print("    • Production: oscillation mixing with active neutrinos")
print("    • Detection: X-ray line at ~3.5 keV (disputed)")
print()
print("  4. BARYONIC MATTER:")
print("  Only ~5% of cosmic matter is baryonic (protons, neutrons). Most")
print("  exists as:")
print("    • Ionized gas (WHIM): 40% (warm-hot intergalactic medium)")
print("    • Neutral hydrogen: 30% (21cm surveys)")
print("    • Stars: 10% (galaxies)")
print("    • Planets, dust: <1%")
print()
print("  5. MATTER-DARK ENERGY EQUALITY:")
print("  Matter density scales as ρ_m ~ a⁻³ (dilutes with expansion).")
print("  Dark energy ρ_Λ stays constant. They crossed at z_eq ≈ 0.4:")
print()
print("    a_eq = (Ω_Λ/Ω_m)^(1/3) ≈ 1.55  (about 4 Gyr ago)")
print()
print("  Before: matter-dominated → deceleration")
print("  After: dark energy-dominated → acceleration")
print("  We live just after the transition—a cosmic coincidence!")
print()
print("  TRIPHASE CONNECTION:")
print("  TriPhase derives ρ_m from H_0 ~ α¹⁸ × f_e, connecting matter density")
print("  to the fine structure constant:")
print()
print("    ρ_m ~ Ω_m × H_0² / G ~ α³⁶ × f_e² / G")
print()
print("  This suggests matter density is not an independent parameter but")
print("  emerges from the same electromagnetic vacuum structure (α, ε₀, μ₀)")
print("  that determines atomic physics. The 36th power of α ~ (1/137)³⁶")
print("  provides natural suppression from Compton frequency f_e ~ 10²⁰ Hz")
print("  to cosmic density scale ρ_m ~ 10⁻²⁷ kg/m³.")
print()
print("  Could dark matter itself be an electromagnetic vacuum effect—")
print("  perhaps the 'negative mass' component of EM field energy that")
print("  balances positive mass to give the observed Ω_m ≈ 0.315?")
print("=" * 70)

input("Press Enter to exit...")
