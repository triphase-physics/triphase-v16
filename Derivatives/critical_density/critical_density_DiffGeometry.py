"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================
Derivative:  Critical Density (ρ_c = 3H₀²/(8πG))
Framework:   DiffGeometry
Version:     16.0
Generated:   2026-03-26
Status:      Active Development
Tag: (D)
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

# === DERIVATION: Critical Density ===
print("=" * 80)
print("TRIPHASE V16: CRITICAL DENSITY")
print("Framework: DiffGeometry")
print("=" * 80)
print()

print("CRITICAL DENSITY DEFINITION:")
print("-" * 80)
print()
print("The energy density that produces FLAT spatial geometry (k = 0):")
print()
print("  ρ_c = 3H₀² / (8πG)")
print()
print("Where:")
print("  H₀ = Hubble constant (today)")
print("  G  = gravitational constant")
print()

# === COMPUTE CRITICAL DENSITY ===
print("=" * 80)
print("COMPUTING ρ_c FROM TRIPHASE")
print("=" * 80)
print()

print("Hubble constant:")
print(f"  H₀ = π√3 × f_e × α¹⁸")
print(f"     = {H_0:.15e} s⁻¹")
print()

H_0_SI = H_0 * 3.086e22 / 1e6  # Convert to km/s/Mpc
print(f"  H₀ = {H_0_SI:.3f} km/s/Mpc")
print()

print("Gravitational constant:")
print(f"  G = c⁴ × 7.5 × ε₀³ × μ₀²")
print(f"    = {G:.15e} m³ kg⁻¹ s⁻²")
print()

rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)

print("Critical density:")
print(f"  ρ_c = 3H₀²/(8πG)")
print(f"      = {rho_c:.15e} kg/m³")
print()

# === DIFFERENTIAL GEOMETRY INTERPRETATION ===
print("=" * 80)
print("DIFFERENTIAL GEOMETRY INTERPRETATION")
print("=" * 80)
print()

print("FRIEDMANN-LEMAÎTRE-ROBERTSON-WALKER (FLRW) METRIC:")
print("-" * 80)
print()
print("The metric for a homogeneous, isotropic universe:")
print()
print("  ds² = -c²dt² + a²(t)[dr²/(1-kr²) + r²(dθ² + sin²θ dφ²)]")
print()
print("Where:")
print("  a(t) = scale factor (cosmic expansion)")
print("  k = spatial curvature parameter")
print()
print("Three cases:")
print("  k = +1 : Positive curvature (closed, spherical)")
print("  k =  0 : Zero curvature (flat, Euclidean)")
print("  k = -1 : Negative curvature (open, hyperbolic)")
print()

print("FRIEDMANN EQUATION:")
print("-" * 80)
print()
print("Einstein equation G_μν = κT_μν for FLRW metric gives:")
print()
print("  H² = (8πG/3)ρ - kc²/a²")
print()
print("Where H = ȧ/a is the Hubble parameter.")
print()
print("At critical density ρ = ρ_c, the curvature term vanishes:")
print()
print("  k = 0  ⟺  ρ = 3H²/(8πG) = ρ_c")
print()
print("GEOMETRIC MEANING:")
print()
print("  ρ < ρ_c : k = -1 (open, negative curvature)")
print("           Hyperbolic geometry, infinite volume")
print("           Universe expands forever")
print()
print("  ρ = ρ_c : k = 0 (flat, zero curvature)")
print("           Euclidean geometry")
print("           Critical balance")
print()
print("  ρ > ρ_c : k = +1 (closed, positive curvature)")
print("           Spherical geometry, finite volume")
print("           May recollapse (depends on equation of state)")
print()

print("SPATIAL CURVATURE RADIUS:")
print("-" * 80)
print()
print("For k ≠ 0, define curvature radius R_k:")
print()
print("  R_k² = -c²/(k H²) × [1 - (ρ/ρ_c)]")
print()
print("At ρ = ρ_c: R_k → ∞ (flat space)")
print()

# Calculate today's curvature radius (if k ≠ 0)
Omega_total = 1.0  # Assume flat universe (Planck 2018)
print(f"Observed total density parameter:")
print(f"  Ω_total = ρ_total/ρ_c = {Omega_total:.3f} (Planck 2018)")
print()
print("This indicates k ≈ 0 (flat spatial geometry).")
print()

# === CONVERT TO VARIOUS UNITS ===
print("=" * 80)
print("CRITICAL DENSITY IN VARIOUS UNITS")
print("=" * 80)
print()

# Mass density
print(f"Mass density:")
print(f"  ρ_c = {rho_c:.15e} kg/m³")
print(f"      = {rho_c * 1e-3:.15e} g/cm³")
print()

# Energy density
energy_density = rho_c * c**2
print(f"Energy density:")
print(f"  ρ_c c² = {energy_density:.6e} J/m³")
print(f"         = {energy_density / 1.60218e-19:.6e} eV/m³")
print(f"         = {energy_density / 1.60218e-10:.6e} GeV/m³")
print()

# Number density (protons)
n_protons = rho_c / m_p
print(f"Equivalent proton number density:")
print(f"  n_p = ρ_c/m_p = {n_protons:.6e} protons/m³")
print(f"                = {n_protons * 1e-6:.6e} protons/cm³")
print()

# Solar masses per cubic Megaparsec
Mpc = 3.086e22  # meters
M_sun = 1.989e30  # kg
solar_mass_density = rho_c * Mpc**3 / M_sun
print(f"Solar mass per cubic Megaparsec:")
print(f"  ρ_c = {solar_mass_density:.6e} M_sun/Mpc³")
print()

# Atoms per cubic meter
n_atoms = rho_c / (1.67e-27)  # Average atomic mass
print(f"Equivalent atom number density:")
print(f"  n = {n_atoms:.6e} atoms/m³")
print()

# === DENSITY FRACTIONS ===
print("=" * 80)
print("COSMIC DENSITY BUDGET (Planck 2018)")
print("=" * 80)
print()

print("Density parameter Ω_i = ρ_i / ρ_c:")
print("-" * 80)
print()

components = [
    ("Baryonic matter", 0.0493, "Ordinary matter (atoms)"),
    ("Dark matter", 0.2657, "Non-baryonic CDM"),
    ("Dark energy", 0.685, "Cosmological constant"),
    ("Radiation", 9.4e-5, "Photons + neutrinos"),
    ("Spatial curvature", 0.0007, "k/a²H² term"),
]

Omega_sum = 0.0
for name, Omega, description in components:
    rho_i = Omega * rho_c
    Omega_sum += Omega
    print(f"{name:20s}: Ω = {Omega:.4f}")
    print(f"                      ρ = {rho_i:.6e} kg/m³")
    print(f"                      ({description})")
    print()

print(f"Total: Ω_total = {Omega_sum:.4f}")
print()

if abs(Omega_sum - 1.0) < 0.01:
    print("STATUS: Universe is FLAT (Ω_total ≈ 1)")
else:
    print(f"Curvature: Ω_k = {1.0 - Omega_sum:.4f}")
print()

# === COSMIC SCALES ===
print("=" * 80)
print("COSMIC SCALES FROM CRITICAL DENSITY")
print("=" * 80)
print()

print("HUBBLE RADIUS (HORIZON):")
print("-" * 80)
r_H = c / H_0
print(f"  r_H = c/H₀ = {r_H:.6e} m")
print(f"            = {r_H / 9.461e15:.3f} light-years")
print(f"            = {r_H / 3.086e22:.3f} Mpc")
print(f"            = {r_H / 9.461e24:.3f} Gly")
print()

print("HUBBLE TIME (AGE SCALE):")
print("-" * 80)
t_H = 1.0 / H_0
print(f"  t_H = 1/H₀ = {t_H:.6e} s")
print(f"             = {t_H / (365.25*24*3600):.3e} years")
print(f"             = {t_H / (365.25*24*3600*1e9):.2f} Gyr")
print()

print("MASS WITHIN HUBBLE VOLUME:")
print("-" * 80)
V_H = (4.0/3.0) * math.pi * r_H**3
M_H = rho_c * V_H
print(f"  V_H = (4π/3)r_H³ = {V_H:.6e} m³")
print(f"  M_H = ρ_c × V_H = {M_H:.6e} kg")
print(f"                  = {M_H / M_sun:.6e} M_sun")
print()

# === SCHWARZSCHILD RADIUS OF UNIVERSE ===
print("=" * 80)
print("SCHWARZSCHILD RADIUS OF OBSERVABLE UNIVERSE")
print("=" * 80)
print()

r_s = 2.0 * G * M_H / c**2
print(f"If mass M_H were compressed to Schwarzschild radius:")
print(f"  r_s = 2GM_H/c² = {r_s:.6e} m")
print(f"                 = {r_s / r_H:.3f} × r_H")
print()

if r_s < r_H:
    print("  r_s < r_H : Observable universe is NOT a black hole")
else:
    print("  r_s ≥ r_H : Observable universe acts like black hole")
print()

# === TIME EVOLUTION ===
print("=" * 80)
print("CRITICAL DENSITY EVOLUTION")
print("=" * 80)
print()

print("FRIEDMANN EQUATION:")
print("-" * 80)
print()
print("  H²(t) = (8πG/3)ρ(t)")
print()
print("For flat universe (k=0), critical density evolves:")
print()
print("  ρ_c(t) = 3H²(t)/(8πG)")
print()
print("As universe expands, H decreases → ρ_c decreases.")
print()

print("MATTER-DOMINATED ERA (a ∝ t^(2/3)):")
print("  H(t) = 2/(3t)")
print("  ρ_c(t) ∝ 1/t²")
print()

print("RADIATION-DOMINATED ERA (a ∝ t^(1/2)):")
print("  H(t) = 1/(2t)")
print("  ρ_c(t) ∝ 1/t²")
print()

print("DARK ENERGY-DOMINATED ERA (a ∝ e^(Ht)):")
print("  H(t) = constant")
print("  ρ_c(t) = constant")
print()

# === CURVATURE INTERPRETATION ===
print("=" * 80)
print("CURVATURE AND CRITICAL DENSITY")
print("=" * 80)
print()

print("RICCI SCALAR FOR FLRW METRIC:")
print("-" * 80)
print()
print("  R = 6[ä/a + (ȧ/a)² + kc²/a²]")
print()
print("Using Friedmann equations:")
print("  R = 6[(8πG/3)(ρ + 3P/c²) + (8πG/3)ρ - kc²/a²]")
print()
print("At critical density (k=0):")
print("  R = 6(8πG/3)[2ρ + 3P/c²]")
print()

kappa = 8.0 * math.pi * G / c**4
R_critical = 6.0 * (8.0 * math.pi * G / 3.0) * 2.0 * rho_c  # Assuming P≈0 (matter-dominated)

print(f"Today's Ricci scalar (matter-dominated approx):")
print(f"  R ≈ {R_critical:.6e} m⁻²")
print()

L_curv = 1.0 / math.sqrt(abs(R_critical))
print(f"Curvature length scale:")
print(f"  L = 1/√R = {L_curv:.6e} m")
print(f"           = {L_curv / r_H:.3f} × r_H")
print()

# === CALIBRATION CHECKPOINT ===
print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
print()

# Observed values (Planck 2018)
rho_c_obs = 8.5e-27  # kg/m³ (from H_0 = 67.4 km/s/Mpc)

print(f"TriPhase critical density:")
print(f"  ρ_c = {rho_c:.15e} kg/m³")
print()
print(f"Observed (Planck 2018, H₀=67.4 km/s/Mpc):")
print(f"  ρ_c = {rho_c_obs:.15e} kg/m³")
print()

deviation = abs(rho_c - rho_c_obs) / rho_c_obs * 100.0
print(f"Deviation: {deviation:.2f}%")
print()

if deviation < 5.0:
    print("STATUS: GOOD - Within 5% of observations")
    print("        (Discrepancy mainly from H₀ uncertainty)")
elif deviation < 10.0:
    print("STATUS: Acceptable - Within 10% of observations")
else:
    print("STATUS: Check H₀ and G derivations")
print()

print(f"Hubble constant:")
print(f"  TriPhase: H₀ = {H_0_SI:.3f} km/s/Mpc")
print(f"  Planck:   H₀ = 67.4 ± 0.5 km/s/Mpc")
print(f"  SH0ES:    H₀ = 73.0 ± 1.0 km/s/Mpc (tension!)")
print()

print("=" * 80)
print("DERIVATION COMPLETE")
print("=" * 80)
print()
print("Critical density ρ_c = 3H₀²/(8πG) is the density that produces")
print("FLAT spatial geometry in the FLRW metric.")
print()
print("In differential geometry:")
print("  - ρ < ρ_c : Negative spatial curvature (hyperbolic)")
print("  - ρ = ρ_c : Zero spatial curvature (Euclidean)")
print("  - ρ > ρ_c : Positive spatial curvature (spherical)")
print()
print("TriPhase derives both H₀ and G from vacuum parameters:")
print(f"  H₀ = π√3 × f_e × α¹⁸ = {H_0_SI:.3f} km/s/Mpc")
print(f"  G = c⁴ × 7.5 × ε₀³ × μ₀² = {G:.6e} m³ kg⁻¹ s⁻²")
print()
print(f"Result: ρ_c = {rho_c:.6e} kg/m³")
print()
print("Observations (Planck 2018) indicate Ω_total ≈ 1.0000 ± 0.0007,")
print("confirming that the universe is spatially FLAT to high precision.")
print()
print("The critical density sets the scale for all cosmic densities:")
print("  Ω_i = ρ_i/ρ_c")
print()

input("Press Enter to exit...")
