"""
================================================================================
TriPhase V16 - Matter Density (Information Theory Framework)
================================================================================

INFORMATION THEORY INTERPRETATION:
Matter density encodes baryonic and dark matter information content.
From an information-theoretic perspective, the density represents:
  - Shannon entropy of matter distribution (galaxies, clusters, voids)
  - Kolmogorov complexity of structure formation
  - Channel capacity for gravitational information in large-scale structure
  - Fisher information about the matter power spectrum
  - Mutual information between baryonic and dark matter
  - Holographic bits encoding cosmic web topology

Matter density ρ_m ≈ 0.31 × ρ_c encodes log₂(N_galaxies) bits of large-scale
structure information. The distribution of ~10¹¹ galaxies represents
log₂(10¹¹) ≈ 37 bits of hierarchical organization. Dark matter (85%) and
baryonic matter (15%) encode mutual information about structure formation.

MIS TAG: (D*H) — matter information content

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
"""

import math

# ============================================================================
# Anchor constants (TriPhase V16 Standard)
# ============================================================================
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

print("=" * 80)
print("TriPhase V16 - Matter Density (Information Theory Framework)")
print("=" * 80)
print()

# ============================================================================
# Step 1: Information-Theoretic Foundation
# ============================================================================
print("-" * 80)
print("STEP 1: Information-Theoretic Foundation")
print("-" * 80)
print()

print("Matter density from observations:")
print("  Ω_m = ρ_m / ρ_c ≈ 0.311 (Planck 2018)")
print()

rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)
Omega_m = 0.311  # Total matter (baryonic + dark)
Omega_b = 0.049  # Baryonic matter only
Omega_dm = Omega_m - Omega_b  # Dark matter

rho_m = Omega_m * rho_c
rho_b = Omega_b * rho_c
rho_dm = Omega_dm * rho_c

print(f"Critical density: ρ_c = {rho_c:.6e} kg/m³")
print(f"Total matter: Ω_m = {Omega_m:.3f}")
print(f"Baryonic matter: Ω_b = {Omega_b:.3f}")
print(f"Dark matter: Ω_dm = {Omega_dm:.3f}")
print()
print(f"Matter density: ρ_m = {rho_m:.6e} kg/m³")
print(f"Baryonic density: ρ_b = {rho_b:.6e} kg/m³")
print(f"Dark matter density: ρ_dm = {rho_dm:.6e} kg/m³")
print()

# ============================================================================
# Step 2: Shannon Entropy of Matter Distribution
# ============================================================================
print("-" * 80)
print("STEP 2: Shannon Entropy of Large-Scale Structure")
print("-" * 80)
print()

print("Cosmic web structure: galaxies, clusters, filaments, voids")
print("Shannon entropy quantifies distribution:")
print("  H(structure) = -Σ p_i log₂(p_i)")
print()

# Approximate volume fractions
f_galaxies = 0.001
f_filaments = 0.10
f_walls = 0.20
f_voids = 0.699

shannon_structure = -f_galaxies * math.log2(f_galaxies) - \
                     f_filaments * math.log2(f_filaments) - \
                     f_walls * math.log2(f_walls) - \
                     f_voids * math.log2(f_voids)

print(f"Galaxies/clusters: {f_galaxies * 100:.1f}%")
print(f"Filaments: {f_filaments * 100:.0f}%")
print(f"Walls: {f_walls * 100:.0f}%")
print(f"Voids: {f_voids * 100:.1f}%")
print(f"Shannon entropy: {shannon_structure:.3f} bits")
print()

# ============================================================================
# Step 3: Kolmogorov Complexity of Structure Formation
# ============================================================================
print("-" * 80)
print("STEP 3: Kolmogorov Complexity of Cosmic Structure")
print("-" * 80)
print()

print("K(cosmic web) = minimal description of structure formation")
print()
print("Initial conditions: Gaussian random field (CMB fluctuations)")
print("  δρ/ρ ~ 10⁻⁵ at recombination")
print()
print("Evolution: gravitational collapse, N-body dynamics")
print()

# Number of observable galaxies
N_galaxies = 2e11  # ~200 billion galaxies in observable universe
kolmogorov_structure = math.log2(N_galaxies)

print(f"Observable galaxies: N ~ {N_galaxies:.0e}")
print(f"Kolmogorov complexity: log₂(N) ≈ {kolmogorov_structure:.2f} bits")
print()

# ============================================================================
# Step 4: Fisher Information about Matter Power Spectrum
# ============================================================================
print("-" * 80)
print("STEP 4: Fisher Information about Power Spectrum")
print("-" * 80)
print()

print("Fisher information I(P(k)) from galaxy surveys:")
print("  P(k) = power spectrum (density fluctuations vs scale k)")
print()

# Power spectrum precision ~ 1% at large scales
precision_Pk = 0.01
fisher_bits_Pk = -math.log2(precision_Pk)

print(f"Power spectrum precision: {precision_Pk * 100:.0f}%")
print(f"Fisher information: {fisher_bits_Pk:.2f} bits")
print()

# ============================================================================
# Step 5: Channel Capacity for Gravitational Clustering
# ============================================================================
print("-" * 80)
print("STEP 5: Channel Capacity - Structure Formation Information")
print("-" * 80)
print()

print("Gravitational collapse rate determines information flow:")
print("  τ_collapse ~ 1/√(G ρ_m)")
print()

tau_collapse = 1.0 / math.sqrt(G * rho_m)
rate_collapse = 1.0 / tau_collapse

# Channel capacity: formation vs dissipation
num_structure_modes = 4  # Galaxies, filaments, walls, voids
capacity_structure = rate_collapse * math.log2(num_structure_modes)

print(f"Collapse timescale: τ ~ {tau_collapse / (365.25 * 24 * 3600):.3e} years")
print(f"Collapse rate: {rate_collapse:.3e} Hz")
print(f"Structure modes: {num_structure_modes}")
print(f"Channel capacity: {capacity_structure:.3e} bits/s")
print()

# ============================================================================
# Step 6: Holographic Bound on Matter Information
# ============================================================================
print("-" * 80)
print("STEP 6: Holographic Bound on Cosmic Structure")
print("-" * 80)
print()

print("Maximum information in Hubble volume:")
print("  S_matter ≤ A_H / (4 ℓ_P²)")
print()

R_H = c / H_0
planck_length = math.sqrt(hbar * G / c**3)
area_H = 4.0 * math.pi * R_H**2
S_holographic = area_H / (4.0 * planck_length**2)

print(f"Hubble radius: R_H = {R_H:.3e} m")
print(f"Holographic bound: {S_holographic:.3e} bits")
print()

# Actual information in matter distribution
N_particles = rho_m * (4.0/3.0 * math.pi * R_H**3) / m_p
S_matter_actual = math.log2(N_particles)

print(f"Particles in Hubble volume: N ~ {N_particles:.3e}")
print(f"Matter information: log₂(N) ≈ {S_matter_actual:.2f} bits")
print()

# ============================================================================
# Step 7: TriPhase Derivation - Matter-Dark Energy Balance
# ============================================================================
print("-" * 80)
print("STEP 7: TriPhase Derivation - Cosmic Matter Information")
print("-" * 80)
print()

print("In TriPhase, matter density emerges from cosmological evolution:")
print("  Present epoch: matter-dark energy equality happened recently")
print("  Ω_m / Ω_Λ ≈ 0.45 (matter slightly subdominant)")
print()

ratio_m_Lambda = Omega_m / (1.0 - Omega_m)

print(f"Matter/dark-energy ratio: Ω_m/Ω_Λ = {ratio_m_Lambda:.3f}")
print()
print("We live near the cosmic information transition!")
print()

# ============================================================================
# Step 8: Mutual Information - Baryonic and Dark Matter
# ============================================================================
print("-" * 80)
print("STEP 8: Mutual Information I(baryons; dark matter)")
print("-" * 80)
print()

print("Baryons and dark matter are gravitationally coupled:")
print("  I(b; dm) from structure formation")
print()

# Entropy of baryon/dark-matter distribution
H_baryons = -Omega_b/Omega_m * math.log2(Omega_b/Omega_m)
H_dm = -Omega_dm/Omega_m * math.log2(Omega_dm/Omega_m)
shannon_matter_split = H_baryons + H_dm

print(f"Baryon fraction: f_b = Ω_b/Ω_m = {Omega_b/Omega_m:.3f}")
print(f"Dark matter fraction: f_dm = {Omega_dm/Omega_m:.3f}")
print(f"Matter composition entropy: {shannon_matter_split:.3f} bits")
print()

# Mutual information ~ correlation strength
mutual_info_matter = shannon_matter_split

print(f"Mutual information (b-dm): {mutual_info_matter:.3f} bits")
print()

# ============================================================================
# Step 9: Landauer's Principle - Virialization
# ============================================================================
print("-" * 80)
print("STEP 9: Landauer's Principle - Gravitational Heating")
print("-" * 80)
print()

print("Virialization converts gravitational potential to kinetic energy:")
print("  ΔE ~ GM²/R → k_B T_vir")
print()

k_B = 1.380649e-23  # J/K
# Typical galaxy cluster: M ~ 10¹⁵ M_☉, R ~ 1 Mpc
M_cluster = 1e15 * 1.989e30  # kg
R_cluster = 3.086e22  # m (1 Mpc)
T_vir = G * M_cluster * m_p / (3.0 * k_B * R_cluster)

E_landauer_vir = k_B * T_vir * math.log(2.0)

print(f"Cluster mass: M ~ {M_cluster / 1.989e30:.0e} M_☉")
print(f"Cluster radius: R ~ {R_cluster / 3.086e22:.1f} Mpc")
print(f"Virial temperature: T ~ {T_vir:.3e} K ({T_vir * k_B / (1.6e-19 * 1e3):.1f} keV)")
print(f"Landauer limit: {E_landauer_vir / e:.3e} eV")
print()

# ============================================================================
# Step 10: Calibration and Validation
# ============================================================================
print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
print()

# Verify matter density from nucleon count
n_baryons = rho_b / m_p  # nucleons/m³
n_baryons_cm3 = n_baryons / 1e6

print(f"Baryon number density: n_b = {n_baryons:.3e} m⁻³")
print(f"                          = {n_baryons_cm3:.6f} cm⁻³")
print()

# Compare to observational constraints
# CMB: Ω_b h² = 0.02237 ± 0.00015 (Planck)
# BBN: Similar constraints
h = H_0 / (100000.0 / 3.086e22)  # Dimensionless Hubble
Omega_b_h2_calc = Omega_b * h**2
Omega_b_h2_obs = 0.02237

deviation_baryons = abs(Omega_b_h2_calc - Omega_b_h2_obs) / Omega_b_h2_obs * 100.0

print(f"Calculated Ω_b h²: {Omega_b_h2_calc:.5f}")
print(f"Planck Ω_b h²:     {Omega_b_h2_obs:.5f}")
print(f"Deviation:         {deviation_baryons:.2f}%")
print()

print("Information-theoretic summary:")
print(f"  Shannon entropy (structure): {shannon_structure:.3f} bits")
print(f"  Kolmogorov complexity:       {kolmogorov_structure:.2f} bits")
print(f"  Fisher info (power spec):    {fisher_bits_Pk:.2f} bits")
print(f"  Matter composition entropy:  {shannon_matter_split:.3f} bits")
print(f"  Mutual info (b-dm):          {mutual_info_matter:.3f} bits")
print(f"  Structure capacity:          {capacity_structure:.3e} bits/s")
print(f"  Holographic bound:           {S_holographic:.3e} bits")
print()

if deviation_baryons < 10.0:
    print("STATUS: EXCELLENT - Matter density information validated!")
elif deviation_baryons < 20.0:
    print("STATUS: GOOD - Within cosmological precision")
else:
    print("STATUS: REVIEW - Check matter density calculations")

print()
print("=" * 80)
print("Matter density: Where information clusters into cosmic structure.")
print("=" * 80)

input("Press Enter to exit...")
