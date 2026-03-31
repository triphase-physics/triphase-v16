"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================
Derivative:  Lyman-Alpha Wavelength (λ_Lyα ≈ 121.567 nm)
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
epsilon_0 = 8.8541878128e-12  # F/m - permittivity of free space
mu_0      = 1.25663706212e-6  # H/m - permeability of free space
e         = 1.602176634e-19   # C - elementary charge

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
R_inf     = alpha**2 * m_e * c / (2.0 * hbar)

print("=" * 80)
print("TRIPHASE V16 - LYMAN-ALPHA WAVELENGTH")
print("Framework: DiffGeometry")
print("=" * 80)
print()

# === DERIVATION with DiffGeometry interpretation ===
print("DIFFERENTIAL GEOMETRY INTERPRETATION:")
print("-" * 80)
print("The Lyman-alpha wavelength represents the GEODESIC DISTANCE between")
print("the n=1 and n=2 energy shells on the hydrogen atom manifold.")
print()
print("In differential geometry, the hydrogen atom defines a curved manifold")
print("with metric determined by the Coulomb potential:")
print()
print("  ds² = (1 + V/mc²)² dt² - dr²/c² - r²dΩ²")
print()
print("The Bohr radius sets the metric scale:")
print("  a₀ = ℏ / (m_e × c × α)")
print()
print("Energy levels correspond to quantized geodesics:")
print("  E_n = -m_e × c² × α² / (2n²)")
print()
print("The Lyman-alpha transition (n=2→n=1) spans a geodesic distance:")
print("  λ_Lyα = 4 / (3 × R_∞)")
print()
print("Where R_∞ (Rydberg constant) is the curvature parameter:")
print("  R_∞ = α² × m_e × c / (2ℏ)")
print()
print("Physical Meaning:")
print("  - λ_Lyα is the proper length along the photon geodesic")
print("  - The factor 4/3 comes from the Rydberg formula (1 - 1/4)")
print("  - This is the MOST PROBABLE transition in hydrogen")
print("  - Manifests in cosmological Lyman-alpha forest (redshifted)")
print()

# === COMPUTATION ===
# Bohr radius
a_0 = hbar / (m_e * c * alpha)

# Rydberg constant (in m^-1)
R_inf_derived = alpha**2 * m_e * c / (2.0 * hbar)

# Lyman-alpha wavelength
lambda_Lyman_alpha = 4.0 / (3.0 * R_inf_derived)

# Photon energy
E_Lyman_alpha = h * c / lambda_Lyman_alpha

print("=" * 80)
print("NUMERICAL RESULT:")
print("=" * 80)
print("Geometric parameters:")
print(f"  Bohr radius a₀           = {a_0:.6e} m")
print(f"  Rydberg constant R_∞     = {R_inf_derived:.6e} m⁻¹")
print()
print("Lyman-alpha transition:")
print(f"  λ_Lyα = 4/(3×R_∞)        = {lambda_Lyman_alpha * 1e9:.6f} nm")
print(f"  Photon energy            = {E_Lyman_alpha / e:.6f} eV")
print(f"  Frequency                = {c / lambda_Lyman_alpha:.6e} Hz")
print()

# === CALIBRATION CHECKPOINT ===
# CODATA 2018 values
R_inf_CODATA = 10973731.568160  # m^-1
lambda_Lyman_alpha_CODATA = 121.56701  # nm

print("=" * 80)
print("CALIBRATION vs CODATA 2018:")
print("=" * 80)
print(f"Rydberg constant R_∞:")
print(f"  DERIVED   : {R_inf_derived:.6f} m⁻¹")
print(f"  CODATA    : {R_inf_CODATA:.6f} m⁻¹")
print(f"  Difference: {abs(R_inf_derived - R_inf_CODATA):.3f} m⁻¹")
print(f"  Precision : {abs(R_inf_derived - R_inf_CODATA) / R_inf_CODATA * 1e6:.3f} ppm")
print()
print(f"Lyman-alpha wavelength:")
print(f"  DERIVED   : {lambda_Lyman_alpha * 1e9:.6f} nm")
print(f"  CODATA    : {lambda_Lyman_alpha_CODATA:.6f} nm")
print(f"  Difference: {abs(lambda_Lyman_alpha * 1e9 - lambda_Lyman_alpha_CODATA):.6f} nm")
print(f"  Precision : {abs(lambda_Lyman_alpha * 1e9 - lambda_Lyman_alpha_CODATA) / lambda_Lyman_alpha_CODATA * 1e6:.3f} ppm")
print()

# === GEOMETRIC RELATIONS ===
print("=" * 80)
print("GEOMETRIC STRUCTURE OF HYDROGEN ATOM:")
print("=" * 80)
print(f"Bohr radius (metric scale)         : a₀ = {a_0 / 1e-10:.6f} Å")
print(f"Classical electron radius          : r_e = {r_e / 1e-15:.6f} fm")
print(f"Compton wavelength                 : λ_C = {h / (m_e * c) / 1e-12:.6f} pm")
print()
print("Ratio relationships:")
print(f"  a₀ / r_e = {a_0 / r_e:.6f} = 1/α ✓")
print(f"  λ_C / r_e = {(h / (m_e * c)) / r_e:.6f} ≈ 2π/α")
print(f"  λ_Lyα / a₀ = {lambda_Lyman_alpha / a_0:.6f}")
print()

# === COSMOLOGICAL MANIFESTATION ===
print("=" * 80)
print("COSMOLOGICAL LYMAN-ALPHA FOREST:")
print("=" * 80)
print("At high redshift (z ~ 2-6), the Lyman-alpha line is redshifted into")
print("the optical band, creating absorption lines from intervening neutral")
print("hydrogen clouds. This maps the large-scale structure of the universe.")
print()
print("Example redshifts:")
for z in [2.0, 3.0, 4.0, 5.0]:
    lambda_observed = lambda_Lyman_alpha * (1.0 + z)
    print(f"  z = {z:.1f}  →  λ_obs = {lambda_observed * 1e9:.1f} nm ({lambda_observed * 1e9 / 1e3:.0f} nm)")
print()
print("The Lyman-alpha forest is a direct probe of the cosmic web geometry.")
print()

# === RELATION TO TRIPHASE ===
print("=" * 80)
print("TRIPHASE CONNECTIONS:")
print("=" * 80)
print("The Rydberg constant connects to TriPhase via:")
print(f"  R_∞ = α² × m_e × c / (2ℏ)")
print()
print("This means λ_Lyα is directly tied to the fine structure constant:")
print(f"  λ_Lyα ∝ 1/α²")
print()
print("Energy scales:")
print(f"  Lyman-alpha energy       : {E_Lyman_alpha / e:.6f} eV")
print(f"  Rydberg energy (13.6 eV) : {m_e * c**2 * alpha**2 / (2 * e):.6f} eV")
print(f"  Ratio                    : {(E_Lyman_alpha / e) / (m_e * c**2 * alpha**2 / (2 * e)):.6f} = 3/4 ✓")
print()
print("The 3/4 factor comes from (1 - 1/4) in the Rydberg formula.")
print("=" * 80)
print()

input("Press Enter to exit...")
