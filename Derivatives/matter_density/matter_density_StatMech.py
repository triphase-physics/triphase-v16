"""
TriPhase V16 — Matter Density (Statistical Mechanics Framework)
================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
Matter density in cosmology refers to the combined density of baryonic matter
(protons, neutrons, electrons) and cold dark matter (CDM), contributing Ω_m ≈ 0.315
to the cosmic density budget. In statistical mechanics, this represents the non-
relativistic component of the cosmic partition function, with equation of state
w = P/ρc² ≈ 0 (pressure-less dust). The partition function for matter at cosmic
temperatures T << m_p c²/k_B includes both thermal kinetic energy and gravitational
potential energy from structure formation. Baryonic matter (Ω_b ≈ 0.049) is well
understood — it's ordinary atoms and molecules. Dark matter (Ω_cdm ≈ 0.266) remains
mysterious but is inferred from gravitational effects on galaxy rotation curves,
cluster dynamics, and gravitational lensing.

The matter density evolves as ρ_m(a) = ρ_m,0 / a³, where a(t) is the cosmic scale
factor. This simple dilution law follows from particle number conservation: as the
universe expands by a factor a, volume increases by a³, so density decreases by a³.
In the grand canonical ensemble, matter has a conserved charge (baryon number for
baryons, some unknown charge for dark matter), so its chemical potential μ evolves
to maintain N = const. The partition function must respect this constraint via
Lagrange multipliers. TriPhase connects matter density to critical density through
Ω_m = ρ_m / ρ_c, where ρ_c depends on H_0 derived from electromagnetic constants.

TAG: (C) — Calibrated to observational data (Ω_m = 0.315 from Planck 2018)
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
print("TriPhase V16: Matter Density (Statistical Mechanics)")
print("=" * 70)
print()

print("STATISTICAL MECHANICS FRAMEWORK")
print("--------------------------------")
print("Grand canonical ensemble: Conserved baryon/DM number")
print("Equation of state: w = P/(ρc²) ≈ 0 (pressure-less dust)")
print("Evolution: ρ_m(a) = ρ_m,0 / a³ (dilution from expansion)")
print("Observable: Ω_m = ρ_m / ρ_c = 0.315 ± 0.007 (Planck 2018)")
print()

print("CRITICAL DENSITY AND MATTER FRACTION")
print("-------------------------------------")
print(f"Hubble constant H_0 = {H_0:.6e} Hz")
print(f"Gravitational constant G = {G:.6e} m³/kg/s²")
print()

# Critical density
rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)
print(f"Critical density ρ_c = 3H_0²/(8πG)")
print(f"  ρ_c = {rho_c:.6e} kg/m³")
print()

# Matter density parameters (Planck 2018)
Omega_b = 0.0486      # Baryonic matter
Omega_cdm = 0.2589    # Cold dark matter
Omega_m = Omega_b + Omega_cdm  # Total matter

print("MATTER DENSITY COMPONENTS (Planck 2018)")
print("---------------------------------------")
print(f"Ω_b (baryons):       {Omega_b:.4f}  ({Omega_b*100:.2f}%)")
print(f"Ω_cdm (dark matter): {Omega_cdm:.4f}  ({Omega_cdm*100:.2f}%)")
print(f"Ω_m (total matter):  {Omega_m:.4f}  ({Omega_m*100:.2f}%)")
print()

# Actual densities
rho_b = Omega_b * rho_c
rho_cdm = Omega_cdm * rho_c
rho_m = Omega_m * rho_c

print("Actual Matter Densities (z = 0, today):")
print(f"  ρ_b (baryons):      {rho_b:.6e} kg/m³")
print(f"  ρ_cdm (dark matter): {rho_cdm:.6e} kg/m³")
print(f"  ρ_m (total matter):  {rho_m:.6e} kg/m³")
print()

# Baryon number density
n_b = rho_b / m_p
print(f"Baryon number density:")
print(f"  n_b = ρ_b / m_p = {n_b:.6e} m^-3")
print(f"  n_b ≈ {n_b:.4f} protons/m³")
print()
print(f"Average interparticle spacing:")
spacing_b = n_b**(-1.0/3.0)
print(f"  d = n_b^(-1/3) = {spacing_b:.3f} m")
print(f"  (Baryons are very sparse in intergalactic space)")
print()

# Dark matter to baryon ratio
DM_to_baryon = Omega_cdm / Omega_b
print(f"Dark matter to baryon ratio:")
print(f"  ρ_cdm / ρ_b = {DM_to_baryon:.2f}")
print(f"  (About 5.3 times more dark matter than baryons by mass)")
print()

# Matter-radiation equality
# ρ_m(a) = ρ_m,0 / a³, ρ_rad(a) = ρ_rad,0 / a⁴
# Equality at a_eq where ρ_m(a_eq) = ρ_rad(a_eq)
Omega_rad = 9.24e-5  # Radiation today (photons + neutrinos)
a_eq = Omega_rad / Omega_m
z_eq = 1.0 / a_eq - 1.0

print("MATTER-RADIATION EQUALITY")
print("-------------------------")
print(f"Ω_rad (today):        {Omega_rad:.5f}")
print(f"Scale factor a_eq:     {a_eq:.6e}")
print(f"Redshift z_eq = 1/a_eq - 1 ≈ {z_eq:.1f}")
print(f"  (Before z ~ 3400, radiation dominated; after, matter dominated)")
print()

# Matter density evolution with redshift
print("MATTER DENSITY EVOLUTION")
print("------------------------")
z_values = [0, 1, 2, 3, 10, 100, 1000]
print("Redshift z | a = 1/(1+z) | ρ_m(z) / ρ_m,0 = (1+z)³")
print("-" * 60)
for z in z_values:
    a = 1.0 / (1.0 + z)
    rho_m_z = rho_m * (1.0 + z)**3
    print(f"  {z:4d}     |  {a:.5f}     |  {(1+z)**3:10.2e}   ({rho_m_z:.3e} kg/m³)")
print()

# Structure formation: Jeans length
# λ_J = √(πc_s² / Gρ_m), where c_s is sound speed
T_IGM = 1e4  # K, typical IGM temperature at z ~ 3
k_B = 1.380649e-23  # J/K
c_s = math.sqrt(k_B * T_IGM / m_p)  # Sound speed
rho_m_z3 = rho_m * (1.0 + 3.0)**3
lambda_J = math.sqrt(math.pi * c_s**2 / (G * rho_m_z3))

print("STRUCTURE FORMATION (Jeans length)")
print("----------------------------------")
print(f"IGM temperature T ~ {T_IGM:.2e} K (z ~ 3)")
print(f"Sound speed c_s = √(k_BT/m_p) = {c_s:.3e} m/s")
print(f"Matter density ρ_m(z=3) = {rho_m_z3:.3e} kg/m³")
print(f"Jeans length λ_J = √(πc_s²/Gρ_m) = {lambda_J:.3e} m")
print(f"  λ_J ≈ {lambda_J / 3.086e22:.2f} Mpc")
print(f"  (Structures smaller than λ_J are pressure-supported)")
print()

# ========== CALIBRATION CHECKPOINT ==========
print("CALIBRATION CHECKPOINT")
print("----------------------")
print("Matter density Ω_m measured by multiple independent probes:")
print()
print("1. Planck CMB (2018):")
print("   Ω_m = 0.3153 ± 0.0073")
print()
print("2. Baryon Acoustic Oscillations (BAO):")
print("   Ω_m = 0.295 ± 0.012 (SDSS/BOSS)")
print()
print("3. Weak Gravitational Lensing:")
print("   Ω_m = 0.30 ± 0.02 (DES Year 3)")
print()
print("4. Galaxy Cluster Abundance:")
print("   Ω_m = 0.29 ± 0.03 (Planck SZ clusters)")
print()
print("5. Type Ia Supernovae + BAO:")
print("   Ω_m = 0.298 ± 0.022 (Pantheon+)")
print()
print(f"TriPhase uses Planck value: Ω_m = {Omega_m:.4f}")
print()
print("Dark matter evidence:")
print("  • Galaxy rotation curves (flat, not Keplerian)")
print("  • Cluster dynamics (virial theorem requires M > visible)")
print("  • Gravitational lensing (mass > luminous matter)")
print("  • CMB acoustic peaks (baryon-photon ratio)")
print("  • Large-scale structure (N-body simulations match obs.)")
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT")
print("-----------------------------")
print("Matter density in cosmology is the non-relativistic component of the cosmic")
print("partition function. At temperatures T << m_p c²/k_B ~ 10^13 K, matter particles")
print("have thermal velocities v << c, giving equation of state:")
print()
print("  P_m ~ ρ_m v² / c² ~ (k_BT / m_p c²) ρ_m c² << ρ_m c²")
print()
print("Thus w = P/(ρc²) ≈ 0 — effectively pressure-less 'dust'. The partition function")
print("for matter is:")
print()
print("  Z_m = ∫d³p exp(-β√(p²c² + m²c⁴))")
print()
print("For non-relativistic particles (pc << mc²), this reduces to:")
print()
print("  Z_m ~ (mk_BT/2πℏ²)^(3/2) exp(βmc²)")
print()
print("The exponential factor exp(βmc²) is enormous at cosmic temperatures (T ~ 10^4 K),")
print("so matter particles are Boltzmann-distributed, not Bose-Einstein or Fermi-Dirac.")
print()
print("Matter density evolution ρ_m(a) = ρ_m,0 / a³ follows from particle conservation:")
print()
print("  N_m = n_m V = const  →  n_m ∝ 1/V ∝ 1/a³  →  ρ_m = m n_m ∝ 1/a³")
print()
print("This contrasts with radiation ρ_rad ∝ 1/a⁴ (extra factor from redshift) and")
print("dark energy ρ_Λ = const. These different scalings determine the cosmic history:")
print()
print("  • z > 3400: Radiation domination (ρ_rad > ρ_m, ρ_Λ)")
print("  • 3400 > z > 0.4: Matter domination (ρ_m > ρ_rad, ρ_Λ)")
print("  • z < 0.4: Dark energy domination (ρ_Λ > ρ_m, ρ_rad)")
print()
print("Structure formation occurs during matter domination. Gravitational instability")
print("amplifies density perturbations δρ/ρ via the Jeans mechanism. The partition")
print("function for a self-gravitating gas includes both thermal and gravitational energy:")
print()
print("  Z = ∫Dρ(x) exp(-β[E_kinetic + E_gravitational])")
print()
print("Density fluctuations grow as δ ∝ a (matter-dominated) or δ ∝ a² (dark energy-")
print("dominated), seeding galaxies, clusters, and cosmic web. N-body simulations that")
print("sample this partition function reproduce observed large-scale structure.")
print()
print("Dark matter (Ω_cdm ≈ 0.26) dominates gravitational dynamics but doesn't interact")
print("electromagnetically. Candidates include WIMPs (weakly interacting massive particles),")
print("axions, sterile neutrinos, or more exotic options. The partition function for dark")
print("matter depends on its unknown microphysics — a major open problem in statistical")
print("mechanics and particle physics.")
print()
print("TriPhase connects matter density to cosmological expansion via:")
print()
print("  ρ_m = Ω_m ρ_c = Ω_m × 3H_0²/(8πG)")
print()
print("where H_0 = π√3 f_e α^18. This links matter density to atomic physics through the")
print("electron frequency f_e and fine structure constant α. If this connection is real,")
print("it suggests the cosmic matter budget may encode fundamental information about")
print("electromagnetic vacuum structure — a profound unification of quantum field theory")
print("and cosmology through statistical mechanics.")
print()
print("=" * 70)

input("Press Enter to exit...")
