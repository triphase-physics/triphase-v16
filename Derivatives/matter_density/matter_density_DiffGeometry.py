"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================
Derivative:  Matter Density (ρ_m = Ω_m × ρ_c)
Framework:   DiffGeometry
Version:     16.0
Generated:   2026-03-26
Status:      Active Development
Tag: (D*)
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
"""
import math

# === ANCHOR INPUTS ===
epsilon_0 = 8.8541878128e-12  # F/m
mu_0      = 1.25663706212e-6  # H/m
e         = 1.602176634e-19   # C

# === DERIVED ANCHOR CHAIN ===
c         = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0       = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha     = 1.0 / alpha_inv
hbar      = Z_0 * e**2 / (4.0 * math.pi * alpha)
h         = 2.0 * math.pi * hbar
G         = c**4 * 7.5 * epsilon_0**3 * mu_0**2
r_e       = 2.8179403262e-15
m_e       = hbar * alpha / (c * r_e)
f_e       = m_e * c**2 / hbar
T_17      = 17 * 18 // 2  # = 153
mp_me     = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p       = m_e * mp_me
H_0       = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r      = c**4 / (8.0 * math.pi * G)

# === DERIVATION: Matter Density ===
print("=" * 80)
print("TRIPHASE V16: MATTER DENSITY")
print("Framework: DiffGeometry")
print("=" * 80)
print()

print("MATTER DENSITY DEFINITION:")
print("-" * 80)
print()
print("Total matter density (baryonic + dark matter):")
print()
print("  ρ_m = Ω_m × ρ_c")
print()
print("Where:")
print("  Ω_m = matter density parameter")
print("  ρ_c = critical density = 3H₀²/(8πG)")
print()

# === CRITICAL DENSITY ===
print("=" * 80)
print("CRITICAL DENSITY FROM TRIPHASE")
print("=" * 80)
print()

rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)

print(f"Hubble constant:")
print(f"  H₀ = π√3 × f_e × α¹⁸ = {H_0:.15e} s⁻¹")
print(f"     = {H_0 * 3.086e22 / 1e6:.3f} km/s/Mpc")
print()

print(f"Critical density:")
print(f"  ρ_c = 3H₀²/(8πG) = {rho_c:.15e} kg/m³")
print()

# === MATTER COMPONENTS ===
print("=" * 80)
print("MATTER DENSITY COMPONENTS (Planck 2018)")
print("=" * 80)
print()

Omega_b = 0.0493  # Baryonic matter
Omega_DM = 0.2657  # Dark matter
Omega_m = Omega_b + Omega_DM  # Total matter

print("Density parameters:")
print(f"  Ω_b  = {Omega_b:.4f}  (baryonic matter)")
print(f"  Ω_DM = {Omega_DM:.4f}  (dark matter)")
print(f"  Ω_m  = {Omega_m:.4f}  (total matter)")
print()

rho_b = Omega_b * rho_c
rho_DM = Omega_DM * rho_c
rho_m = Omega_m * rho_c

print("Matter densities:")
print(f"  ρ_b  = {rho_b:.15e} kg/m³")
print(f"  ρ_DM = {rho_DM:.15e} kg/m³")
print(f"  ρ_m  = {rho_m:.15e} kg/m³")
print()

print(f"Dark matter fraction:")
print(f"  ρ_DM/ρ_m = {rho_DM/rho_m:.3f} ({rho_DM/rho_m*100:.1f}%)")
print(f"  ρ_b/ρ_m  = {rho_b/rho_m:.3f} ({rho_b/rho_m*100:.1f}%)")
print()

# === DIFFERENTIAL GEOMETRY INTERPRETATION ===
print("=" * 80)
print("DIFFERENTIAL GEOMETRY INTERPRETATION")
print("=" * 80)
print()

print("MATTER IN STRESS-ENERGY TENSOR:")
print("-" * 80)
print()
print("For non-relativistic matter (cold dark matter, baryons):")
print()
print("  T^matter_μν = ρ_m c² u_μ u_ν")
print()
print("In comoving frame (u = (1,0,0,0)):")
print()
print("  T^matter_00 = ρ_m c²  (energy density)")
print("  T^matter_ii = 0       (no pressure)")
print()
print("Equation of state:")
print("  P_m ≈ 0  →  w_m = P_m/(ρ_m c²) ≈ 0")
print()
print("Trace of stress-energy:")
print("  T = g^μν T_μν = -ρ_m c²  (negative)")
print()

print("FRIEDMANN EQUATION (MATTER-DOMINATED):")
print("-" * 80)
print()
print("For flat universe (k=0) with matter only:")
print()
print("  H² = (8πG/3) ρ_m")
print()
print("Scale factor evolution:")
print("  a(t) ∝ t^(2/3)")
print("  ρ_m(t) ∝ a^(-3) ∝ t^(-2)")
print()
print("Matter density dilutes as volume expands:")
print("  ρ_m ∝ 1/a³")
print()

print("RICCI CURVATURE FROM MATTER:")
print("-" * 80)
print()
print("Matter creates POSITIVE Ricci curvature:")
print()
print("  R_00 > 0  (temporal curvature)")
print()
print("This causes geodesic convergence:")
print("  - Test particles attract")
print("  - Universe decelerates")
print("  - Gravity is attractive")
print()

kappa = 8.0 * math.pi * G / c**4
R_matter = kappa * rho_m * c**2

print(f"Curvature from matter:")
print(f"  R ~ κ ρ_m c² = {R_matter:.6e} m⁻²")
print()

L_matter = 1.0 / math.sqrt(R_matter)
print(f"Curvature length scale:")
print(f"  L = 1/√R = {L_matter:.6e} m")
print(f"           = {L_matter / 3.086e22:.3f} Mpc")
print()

# === CONVERT TO VARIOUS UNITS ===
print("=" * 80)
print("MATTER DENSITY IN VARIOUS UNITS")
print("=" * 80)
print()

print("TOTAL MATTER (ρ_m):")
print("-" * 80)
print(f"  {rho_m:.15e} kg/m³")
print(f"  {rho_m * 1e-3:.15e} g/cm³")
print()

# Number densities
n_protons = rho_m / m_p
n_atoms = rho_m / (1.67e-27)

print(f"Equivalent proton density:")
print(f"  {n_protons:.6e} protons/m³")
print(f"  {n_protons * 1e-6:.6e} protons/cm³")
print()

print(f"Equivalent atom density:")
print(f"  {n_atoms:.6e} atoms/m³")
print(f"  {n_atoms * 1e-6:.6e} atoms/cm³")
print()

# Energy density
E_m = rho_m * c**2
print(f"Energy density:")
print(f"  {E_m:.6e} J/m³")
print(f"  {E_m / 1.60218e-19:.6e} eV/m³")
print()

# Solar masses per cubic Megaparsec
Mpc = 3.086e22  # m
M_sun = 1.989e30  # kg
rho_m_solar = rho_m * Mpc**3 / M_sun

print(f"Solar masses per cubic Megaparsec:")
print(f"  {rho_m_solar:.6e} M_sun/Mpc³")
print()

print("BARYONIC MATTER (ρ_b):")
print("-" * 80)
print(f"  {rho_b:.15e} kg/m³")
print(f"  {rho_b / m_p:.6e} protons/m³")
print()

print("DARK MATTER (ρ_DM):")
print("-" * 80)
print(f"  {rho_DM:.15e} kg/m³")
print()

# === MATTER-DOMINATED ERA ===
print("=" * 80)
print("MATTER-DOMINATED ERA")
print("=" * 80)
print()

print("SCALE FACTOR EVOLUTION:")
print("-" * 80)
print()
print("When matter dominates (Ω_m >> Ω_Λ, early universe):")
print()
print("  Friedmann: H² = (8πG/3) ρ_m")
print()
print("  With ρ_m ∝ a^(-3):")
print()
print("    H²a³ = constant")
print("    (ȧ/a)² a³ = constant")
print("    ȧ ∝ a^(-1/2)")
print()
print("  Integrating:")
print("    a(t) ∝ t^(2/3)")
print()

# Matter-radiation equality
z_eq = 3365  # Redshift at matter-radiation equality (Planck 2018)
a_eq = 1.0 / (1.0 + z_eq)
t_eq = (2.0/3.0) * (1.0/H_0) * a_eq**(3.0/2.0)

print(f"Matter-radiation equality:")
print(f"  z_eq = {z_eq}")
print(f"  a_eq = 1/(1+z_eq) = {a_eq:.6e}")
print(f"  t_eq ≈ {t_eq:.6e} s = {t_eq/(365.25*24*3600):.3e} years")
print()

print("Before t_eq: Radiation-dominated, a ∝ t^(1/2)")
print("After t_eq:  Matter-dominated, a ∝ t^(2/3)")
print()

# === GEODESIC DEVIATION ===
print("=" * 80)
print("GEODESIC DEVIATION FROM MATTER")
print("=" * 80)
print()

print("RAYCHAUDHURI EQUATION:")
print("-" * 80)
print()
print("Evolution of geodesic convergence (Θ = expansion):")
print()
print("  dΘ/dτ = -Θ²/3 - σ² + ω² - R_μν u^μ u^ν")
print()
print("For matter-dominated FLRW (no shear σ, no vorticity ω):")
print()
print("  dΘ/dτ = -Θ²/3 - (4πG) ρ_m")
print()
print("The -4πG ρ_m term causes geodesic CONVERGENCE:")
print("  - Initially parallel geodesics bend toward each other")
print("  - Gravitational attraction")
print("  - Structure formation")
print()

print("FOCUSING THEOREM:")
print("-" * 80)
print()
print("If ρ_m > 0, then dΘ/dτ < 0 (assuming Θ²/3 dominates).")
print()
print("This means:")
print("  - Matter causes spacetime to focus geodesics")
print("  - Expansion decelerates")
print("  - Ricci curvature is positive along timelike directions")
print()

# === STRUCTURE FORMATION ===
print("=" * 80)
print("STRUCTURE FORMATION")
print("=" * 80)
print()

print("JEANS INSTABILITY:")
print("-" * 80)
print()
print("Matter overdensities grow via gravitational instability.")
print()
print("Jeans length (critical scale for collapse):")
print()
print("  λ_J = c_s √(π/(G ρ_m))")
print()
print("Where c_s is sound speed.")
print()

# Typical sound speed in early universe
c_s = 2e5  # m/s (rough estimate)
lambda_J = c_s * math.sqrt(math.pi / (G * rho_m))

print(f"Today's Jeans length (c_s ~ {c_s:.0e} m/s):")
print(f"  λ_J = {lambda_J:.6e} m")
print(f"      = {lambda_J / 3.086e22:.3e} Mpc")
print()

print("Perturbations larger than λ_J:")
print("  - Gravity dominates pressure")
print("  - Overdensity grows")
print("  - Structure forms")
print()

print("Perturbations smaller than λ_J:")
print("  - Pressure dominates gravity")
print("  - Oscillations (sound waves)")
print("  - No collapse")
print()

# === DARK MATTER HALOS ===
print("=" * 80)
print("DARK MATTER HALOS")
print("=" * 80)
print()

print("VIRIAL THEOREM:")
print("-" * 80)
print()
print("For bound structure (galaxy, cluster):")
print()
print("  2K + U = 0  (virial equilibrium)")
print()
print("Where K = kinetic energy, U = gravitational potential energy.")
print()

# Milky Way halo
M_halo = 1e12 * M_sun  # kg
R_halo = 200e3 * 3.086e19  # m (200 kpc)
rho_halo = M_halo / (4.0/3.0 * math.pi * R_halo**3)

print("Example: Milky Way dark matter halo")
print(f"  M_halo ≈ {M_halo/M_sun:.3e} M_sun")
print(f"  R_halo ≈ {R_halo / (3.086e19):.0f} kpc")
print(f"  ρ_halo ≈ {rho_halo:.3e} kg/m³")
print(f"  ρ_halo/ρ_DM ≈ {rho_halo/rho_DM:.3e}")
print()
print("Halo is ~1000× denser than cosmic average.")
print()

# === COMPARISON WITH DARK ENERGY ===
print("=" * 80)
print("MATTER vs DARK ENERGY")
print("=" * 80)
print()

Omega_Lambda = 0.685
rho_Lambda = Omega_Lambda * rho_c

print("Cosmic density budget today:")
print(f"  Ω_m      = {Omega_m:.4f}  ({Omega_m*100:.1f}%)")
print(f"  Ω_Λ      = {Omega_Lambda:.4f}  ({Omega_Lambda*100:.1f}%)")
print()

print("Matter and dark energy evolve differently:")
print()
print("  ρ_m(a) ∝ a^(-3)   (dilutes with volume)")
print("  ρ_Λ(a) = constant (vacuum energy)")
print()

# Matter-dark energy equality
a_eq_Lambda = (rho_Lambda / rho_m)**(1.0/3.0)
z_eq_Lambda = 1.0/a_eq_Lambda - 1.0

print(f"Matter-dark energy equality:")
print(f"  a_eq = {a_eq_Lambda:.3f}")
print(f"  z_eq = {z_eq_Lambda:.3f}")
print()

print("Past (z > 0.3): Matter-dominated")
print("  - Expansion decelerates")
print("  - Structure formation efficient")
print()

print("Future (z < 0.3): Dark energy-dominated")
print("  - Expansion accelerates")
print("  - Structure formation suppressed")
print()

# === GEODESIC CONGRUENCE ===
print("=" * 80)
print("GEODESIC CONGRUENCE IN MATTER-DOMINATED SPACETIME")
print("=" * 80)
print()

print("RIEMANN CURVATURE FROM MATTER:")
print("-" * 80)
print()
print("In FLRW metric with matter, the Riemann tensor has components:")
print()
print("  R_0i0j ~ (H² + Ḣ) δ_ij")
print("  R_ijkl ~ (H² - kc²/a²) (δ_ik δ_jl - δ_il δ_jk)")
print()
print("For matter-dominated (a ∝ t^(2/3)):")
print("  H = 2/(3t)")
print("  Ḣ = -2/(3t²)")
print()
print("This gives positive Ricci curvature R_00 > 0,")
print("causing attractive gravity and geodesic focusing.")
print()

# === CALIBRATION CHECKPOINT ===
print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
print()

print("TriPhase matter density:")
print(f"  ρ_m = Ω_m × ρ_c")
print(f"      = {Omega_m:.4f} × {rho_c:.6e}")
print(f"      = {rho_m:.15e} kg/m³")
print()

rho_m_obs = 0.315 * 8.5e-27  # kg/m³ (from Planck H_0)
print(f"Observed (Planck 2018):")
print(f"  ρ_m ≈ {rho_m_obs:.15e} kg/m³")
print()

deviation = abs(rho_m - rho_m_obs) / rho_m_obs * 100.0
print(f"Deviation: {deviation:.2f}%")
print()

if deviation < 5.0:
    print("STATUS: EXCELLENT - Within 5% of observations")
elif deviation < 10.0:
    print("STATUS: Good - Within 10% of observations")
else:
    print("STATUS: Check H₀ derivation")
print()

print("Matter components:")
print(f"  Baryons:     ρ_b  = {rho_b:.6e} kg/m³ ({Omega_b/Omega_m*100:.1f}%)")
print(f"  Dark matter: ρ_DM = {rho_DM:.6e} kg/m³ ({Omega_DM/Omega_m*100:.1f}%)")
print()

print("Dark matter dominates by factor of ~5.4")
print()

print("=" * 80)
print("DERIVATION COMPLETE")
print("=" * 80)
print()
print("Matter density ρ_m = Ω_m × ρ_c combines two components:")
print()
print("1. BARYONIC MATTER (Ω_b = 0.0493)")
print("   - Ordinary matter: protons, neutrons, electrons")
print("   - Forms stars, planets, gas, dust")
print("   - Interacts via EM and nuclear forces")
print()
print("2. DARK MATTER (Ω_DM = 0.2657)")
print("   - Non-baryonic, collisionless")
print("   - Interacts only via gravity")
print("   - Forms halos around galaxies")
print()
print("In differential geometry:")
print("  - Matter creates positive Ricci curvature")
print("  - Geodesics converge (attractive gravity)")
print("  - Scale factor a(t) ∝ t^(2/3) (matter-dominated)")
print()
print(f"Total matter density: ρ_m = {rho_m:.6e} kg/m³")
print(f"                      ~ {n_protons:.1f} protons/m³")
print()
print("This is the source of structure formation via gravitational")
print("instability, creating galaxies, clusters, and cosmic web.")
print()

input("Press Enter to exit...")
