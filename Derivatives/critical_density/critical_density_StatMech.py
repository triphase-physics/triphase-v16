"""
TriPhase V16 — Critical Density (Statistical Mechanics Framework)
==================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
Critical density ρ_c = 3H_0²/(8πG) is the energy density required for the universe
to be spatially flat (k = 0 in the Friedmann equation). In statistical mechanics,
this represents the total energy density of all components in the cosmic partition
function: matter (baryons + dark matter), radiation (photons + neutrinos), and
dark energy (cosmological constant or quintessence). The density parameters Ω_i = ρ_i/ρ_c
satisfy the closure relation Σ Ω_i = 1 for a flat universe, which is confirmed
observationally to within 0.4% (Planck 2018). This flatness is not coincidental
but a prediction of cosmic inflation, which stretches any initial curvature to
negligible levels.

The grand canonical ensemble for cosmology assigns chemical potentials to conserved
charges (baryon number, lepton number) and temperature to each component. The
critical density emerges from the Friedmann equation as the energy density where
the kinetic energy of expansion exactly balances the gravitational potential energy,
analogous to the escape velocity problem in classical mechanics. For ρ > ρ_c, the
universe is closed and will eventually recollapse (Big Crunch); for ρ < ρ_c, it's
open and expands forever; ρ = ρ_c is the knife-edge case of flat geometry. TriPhase
connects ρ_c to fundamental constants via H_0 = π√3 f_e α^18, linking cosmological
density to atomic physics.

TAG: (D) — Direct TriPhase derivation from pure wave mechanics
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

# ========== STATISTICAL MECHANICS DERIVATION ==========
print("=" * 70)
print("TriPhase V16: Critical Density (Statistical Mechanics)")
print("=" * 70)
print()

print("STATISTICAL MECHANICS FRAMEWORK")
print("--------------------------------")
print("Grand canonical ensemble: All cosmic energy components")
print("Friedmann equation: H² = (8πG/3)ρ - kc²/a²")
print("Flatness condition: k = 0 requires ρ = ρ_c")
print("Closure relation: Σ Ω_i = 1 (Ω_i = ρ_i / ρ_c)")
print()

print("CRITICAL DENSITY CALCULATION")
print("----------------------------")
print(f"Hubble constant H_0 = {H_0:.6e} Hz")
print(f"Gravitational constant G = {G:.6e} m³/kg/s²")
print(f"Speed of light c = {c:.6e} m/s")
print()

# Critical density
rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)

print("Critical density ρ_c = 3H_0² / (8πG)")
print(f"  ρ_c = {rho_c:.6e} kg/m³")
print(f"  ρ_c = {rho_c * c**2:.6e} J/m³")
print()

# Convert to other units
# Number density of protons at critical density
n_c_protons = rho_c / m_p
print(f"Equivalent proton number density:")
print(f"  n_c = ρ_c / m_p = {n_c_protons:.6e} m^-3")
print(f"  n_c ≈ {n_c_protons:.2f} protons/m³ (about 6 per cubic meter)")
print()

# Energy density in eV/cm³
eV_per_J = 1.0 / 1.602176634e-19
cm3_per_m3 = 1e6
rho_c_eV_cm3 = (rho_c * c**2) * eV_per_J / cm3_per_m3
print(f"  ρ_c c² ≈ {rho_c_eV_cm3:.3f} eV/cm³")
print()

# COSMIC DENSITY BUDGET (Planck 2018 values)
print("COSMIC DENSITY BUDGET (Planck 2018)")
print("------------------------------------")
Omega_b = 0.0486      # Baryonic matter
Omega_cdm = 0.2589    # Cold dark matter
Omega_m = Omega_b + Omega_cdm  # Total matter
Omega_Lambda = 0.6847  # Dark energy
Omega_rad = 9.24e-5   # Radiation (photons + neutrinos)
Omega_k = 0.0007      # Curvature (residual, consistent with zero)
Omega_total = Omega_m + Omega_Lambda + Omega_rad

print(f"Ω_b (baryons):        {Omega_b:.4f}  ({Omega_b*100:.2f}%)")
print(f"Ω_cdm (dark matter):  {Omega_cdm:.4f}  ({Omega_cdm*100:.2f}%)")
print(f"Ω_m (total matter):   {Omega_m:.4f}  ({Omega_m*100:.2f}%)")
print(f"Ω_Λ (dark energy):    {Omega_Lambda:.4f}  ({Omega_Lambda*100:.2f}%)")
print(f"Ω_rad (radiation):    {Omega_rad:.5f}  ({Omega_rad*100:.4f}%)")
print(f"Ω_k (curvature):      {Omega_k:.4f}  ({Omega_k*100:.2f}%)")
print()
print(f"Ω_total = {Omega_total:.4f}")
print(f"Deviation from flatness: {abs(Omega_total - 1.0)*100:.2f}%")
print()

# Actual densities
rho_b = Omega_b * rho_c
rho_cdm = Omega_cdm * rho_c
rho_m = Omega_m * rho_c
rho_Lambda = Omega_Lambda * rho_c
rho_rad = Omega_rad * rho_c

print("Actual Densities:")
print(f"  ρ_b (baryons):      {rho_b:.6e} kg/m³")
print(f"  ρ_cdm (dark matter): {rho_cdm:.6e} kg/m³")
print(f"  ρ_Λ (dark energy):  {rho_Lambda:.6e} kg/m³")
print(f"  ρ_rad (radiation):  {rho_rad:.6e} kg/m³")
print()

# Baryon number density
n_b = rho_b / m_p
print(f"Baryon number density:")
print(f"  n_b = ρ_b / m_p = {n_b:.6e} m^-3")
print(f"  n_b ≈ {n_b:.3f} protons/m³")
print(f"  (This is the average density of hydrogen in the universe)")
print()

# ========== CALIBRATION CHECKPOINT ==========
print("CALIBRATION CHECKPOINT")
print("----------------------")
print("Critical density and flatness have been measured by multiple probes:")
print()
print("1. Planck CMB (2018):")
print("   Ω_total = 1.0007 ± 0.0019 (flat to 0.2%)")
print()
print("2. Baryon Acoustic Oscillations (BAO):")
print("   Ω_m = 0.315 ± 0.007")
print()
print("3. Type Ia Supernovae:")
print("   Ω_Λ = 0.684 ± 0.020")
print()
print("4. Weak Lensing + Galaxy Clustering:")
print("   Ω_m = 0.30 ± 0.02")
print()
print("All measurements converge on Ω_total ≈ 1, confirming flat geometry.")
print("This is a prediction of cosmic inflation (Guth 1981, Linde 1982).")
print()

# Hubble constant from TriPhase
Mpc = 3.0857e22  # meters
H_0_cosmo = H_0 * Mpc / 1000.0  # km/s/Mpc
print(f"TriPhase Hubble constant: H_0 = {H_0_cosmo:.4f} km/s/Mpc")
print(f"TriPhase critical density: ρ_c = {rho_c:.6e} kg/m³")
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT")
print("-----------------------------")
print("Critical density emerges from the Friedmann equation, which is itself a")
print("thermodynamic statement. Starting from Einstein's field equations:")
print()
print("  G_μν = (8πG/c⁴) T_μν")
print()
print("For a homogeneous, isotropic universe (FRW metric), the time-time component gives:")
print()
print("  H² + k/a² = (8πG/3) ρ")
print()
print("where H = ȧ/a is the Hubble parameter, k is spatial curvature, and ρ is the")
print("total energy density. For k = 0 (flat space), this simplifies to:")
print()
print("  ρ_c = 3H² / (8πG)")
print()
print("This is analogous to the virial theorem in statistical mechanics: kinetic energy")
print("(expansion) balances potential energy (gravity) when ρ = ρ_c. For ρ > ρ_c,")
print("gravity wins → closed universe. For ρ < ρ_c, expansion wins → open universe.")
print()
print("In the grand canonical ensemble for cosmology, the partition function includes:")
print()
print("  Z = Σ_states exp(-β[E_matter + E_radiation + E_Λ])")
print()
print("Each component contributes to the total density:")
print("  • Matter: ρ_m = (Ω_b + Ω_cdm) ρ_c ≈ 31% (particles, massive)")
print("  • Radiation: ρ_rad = Ω_rad ρ_c ≈ 0.01% (photons, neutrinos)")
print("  • Dark energy: ρ_Λ = Ω_Λ ρ_c ≈ 68% (vacuum energy)")
print()
print("The closure relation Ω_total = 1 is an empirical fact — the universe is flat.")
print("This was not guaranteed a priori; it's a consequence of cosmic inflation, which")
print("exponentially expands the universe by a factor e^60 or more, diluting any initial")
print("curvature to undetectable levels. Inflation is itself a statistical phenomenon:")
print("a scalar field (inflaton) in a metastable vacuum state undergoes a phase transition,")
print("releasing latent heat that drives exponential expansion.")
print()
print("TriPhase connects critical density to atomic physics through:")
print()
print("  H_0 = π√3 f_e α^18")
print()
print("This implies ρ_c ~ f_e² α^36 / G. The 36th power of α (≈ 10^-78) generates the")
print("enormous ratio between atomic and cosmic scales. This suggests the universe's")
print("total energy density may be encoded in electromagnetic vacuum structure — a")
print("profound hint that cosmology and quantum field theory are deeply connected.")
print()
print("The critical density is also the boundary between eternal expansion and recollapse.")
print("With Ω_Λ > 0, the universe accelerates, approaching de Sitter space asymptotically.")
print("The far-future partition function is that of a pure vacuum state at T_dS ~ 10^-30 K,")
print("the coldest possible temperature. All matter decays, stars die, black holes evaporate,")
print("leaving only quantum fluctuations in empty space — the ultimate heat death, governed")
print("by the statistical mechanics of the vacuum itself.")
print()
print("=" * 70)

input("Press Enter to exit...")
