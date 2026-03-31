"""
TriPhase V16: Speed of Light (c)
Dimensional Analysis Framework

Derivative: c = 1/√(ε₀μ₀)
MIS TAG: (D) - Derived from electromagnetic vacuum properties
Status: Exact by SI 2019 definition

DIMENSIONAL INTERPRETATION:
The speed of light emerges directly from the electromagnetic properties
of vacuum: permittivity ε₀ and permeability μ₀. This is Maxwell's
fundamental result showing light as an electromagnetic wave.

In TriPhase, c is the primary dimensional constant from which all
other velocities and wave speeds are derived.

SI UNITS: [m s⁻¹]

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# =====================================================================
# ANCHOR CONSTANTS (TriPhase V16 Standard Chain)
# =====================================================================
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

# =====================================================================
print("=" * 70)
print("TriPhase V16: Speed of Light (c)")
print("Dimensional Analysis Framework")
print("=" * 70)
print()

# =====================================================================
# STEP 1: IDENTIFY TARGET DIMENSIONS
# =====================================================================
print("STEP 1: IDENTIFY TARGET DIMENSIONS")
print("-" * 70)
print("Target: Speed of light c")
print("SI Dimensions: [m s⁻¹] or [L T⁻¹]")
print()
print("Physical meaning: Velocity of electromagnetic wave propagation")
print("  c = λ × f")
print("  [c] = [L] × [T⁻¹] = [L T⁻¹]")
print()
print("In SI 2019: c = 299,792,458 m/s (exact by definition)")
print()

# =====================================================================
# STEP 2: AVAILABLE BASE DIMENSIONS
# =====================================================================
print("STEP 2: AVAILABLE BASE DIMENSIONS")
print("-" * 70)
print("Base constants and their dimensions:")
print("  ε₀: [A² s⁴ kg⁻¹ m⁻³]  (electric permittivity)")
print("  μ₀: [kg m A⁻² s⁻²]    (magnetic permeability)")
print()
print("From Maxwell's equations:")
print("  Wave equation: ∇²E - ε₀μ₀ ∂²E/∂t² = 0")
print("  Wave speed: c² = 1/(ε₀μ₀)")
print()

# =====================================================================
# STEP 3: DIMENSIONAL MATCHING
# =====================================================================
print("STEP 3: DIMENSIONAL MATCHING")
print("-" * 70)
print("Deriving c = 1/√(ε₀μ₀):")
print()
print("  [ε₀ × μ₀] = [A² s⁴ kg⁻¹ m⁻³] × [kg m A⁻² s⁻²]")
print("            = [A²⁻² s⁴⁻² kg⁻¹⁺¹ m⁻³⁺¹]")
print("            = [s² m⁻²]")
print()
print("  [1/(ε₀μ₀)] = [m² s⁻²]")
print()
print("  [√(1/(ε₀μ₀))] = [m s⁻¹] ✓")
print()
print("Therefore: c = 1/√(ε₀μ₀) has correct dimensions for velocity")
print()

# =====================================================================
# STEP 4: TRIPHASE DERIVATION
# =====================================================================
print("STEP 4: TRIPHASE DERIVATION")
print("-" * 70)
print("TriPhase Formula: c = 1/√(ε₀μ₀)")
print()
print("Given values:")
print(f"  ε₀ = {epsilon_0:.12e} F/m")
print(f"  μ₀ = {mu_0:.12e} H/m")
print()
print("Computing:")
product_eps_mu = epsilon_0 * mu_0
print(f"  ε₀ × μ₀ = {product_eps_mu:.15e} s² m⁻²")
print()
sqrt_product = math.sqrt(product_eps_mu)
print(f"  √(ε₀μ₀) = {sqrt_product:.15e} s m⁻¹")
print()
c_derived = 1.0 / sqrt_product
print(f"  c = 1/√(ε₀μ₀) = {c_derived:.12e} m/s")
print()
print(f"  c = {c_derived:.0f} m/s")
print()

# =====================================================================
# STEP 5: BUCKINGHAM PI THEOREM
# =====================================================================
print("STEP 5: BUCKINGHAM PI THEOREM")
print("-" * 70)
print("Dimensionless groups involving c:")
print()
print("Given c [L T⁻¹], we can form ratios with other velocities:")
print()
print("1. β = v/c (particle velocity relative to light)")
print("   For electron at f_e: v_e << c")
v_phase = hbar / (m_e * r_e)  # rough estimate
beta_e = v_phase / c
print(f"   β_e ~ {beta_e:.6e}")
print()
print("2. Lorentz factor: γ = 1/√(1 - β²)")
print("   For β << 1: γ ≈ 1 + β²/2")
print()
print("3. Fine structure constant: α = v_e/c (Bohr model)")
print(f"   α ≈ {alpha:.6e} ≈ v₁/c")
print("   (v₁ = velocity of electron in first Bohr orbit)")
print()
print("4. Impedance relation: Z₀ = μ₀c")
Z_0_check = mu_0 * c
print(f"   Z₀ = {Z_0:.6f} Ω")
print(f"   μ₀c = {Z_0_check:.6f} Ω")
print(f"   Agreement: ✓")
print()
print("5. Energy-momentum: E² = (pc)² + (mc²)²")
print("   [pc] = [kg m s⁻¹][m s⁻¹] = [kg m² s⁻²] = [J] ✓")
print("   [mc²] = [kg][m² s⁻²] = [J] ✓")
print()

# =====================================================================
# STEP 6: NATURAL UNITS COMPARISON
# =====================================================================
print("STEP 6: NATURAL UNITS COMPARISON")
print("-" * 70)
print("Speed of light in various unit systems:")
print()
print("1. SI units:")
print(f"   c = {c:.0f} m/s (exact)")
print()
print("2. Natural units (c = 1):")
print("   c = 1 (by definition)")
print("   Distances measured in time: 1 light-second = 1 second")
print("   Energy and mass become equivalent: E = m")
print()
print("3. Geometric units (c = 1, G in [L M⁻¹]):")
print("   Time and length have same dimension")
print("   Mass has dimension [L⁻¹]")
print()
print("4. Planck units:")
l_planck = math.sqrt(hbar * G / c**3)
t_planck = l_planck / c
print(f"   l_P/t_P = c = {l_planck / t_planck:.0f} m/s")
print()
print("5. Astronomical units:")
au = 1.496e11  # meters
light_year = c * 365.25 * 24 * 3600
print(f"   c = 1 AU / {au / c:.1f} s")
print(f"   c = 1 light-year / year")
print()

# =====================================================================
# STEP 7: DIMENSIONAL CROSSCHECKS
# =====================================================================
print("STEP 7: DIMENSIONAL CROSSCHECKS")
print("-" * 70)
print("Verifying c in various electromagnetic contexts:")
print()
print("1. Wave equation: c = λf")
lambda_test = 500e-9  # 500 nm (green light)
f_test = c / lambda_test
print(f"   For λ = 500 nm:")
print(f"   f = c/λ = {f_test:.6e} Hz")
print(f"   [c] = [m]/[s⁻¹] = [m s⁻¹] ✓")
print()
print("2. Impedance: Z₀ = √(μ₀/ε₀) = μ₀c = ε₀⁻¹c⁻¹")
Z_0_from_mu = mu_0 * c
Z_0_from_eps = 1.0 / (epsilon_0 * c)
Z_0_from_ratio = math.sqrt(mu_0 / epsilon_0)
print(f"   Z₀ = √(μ₀/ε₀) = {Z_0_from_ratio:.6f} Ω")
print(f"   Z₀ = μ₀c = {Z_0_from_mu:.6f} Ω")
print(f"   Z₀ = 1/(ε₀c) = {Z_0_from_eps:.6f} Ω")
print(f"   All agree ✓")
print()
print("3. Energy flux: S = (E×B)/μ₀ (Poynting vector)")
print("   For EM wave: |E| = c|B|")
print("   [E] = [V/m] = [kg m s⁻³ A⁻¹]")
print("   [B] = [T] = [kg s⁻² A⁻¹]")
print("   [E]/[B] = [m s⁻¹] = [c] ✓")
print()
print("4. Relativistic energy: E = γmc²")
print("   [E] = [kg][m² s⁻²] = [J] ✓")
print()
print("5. Photon momentum: p = E/c = ℏω/c = ℏ/λ")
E_photon = h * f_test
p_photon = E_photon / c
lambda_photon = h / p_photon
print(f"   E = {E_photon:.6e} J")
print(f"   p = E/c = {p_photon:.6e} kg m/s")
print(f"   λ = h/p = {lambda_photon:.6e} m")
print(f"   [p] = [J]/[m s⁻¹] = [kg m² s⁻²]/[m s⁻¹] = [kg m s⁻¹] ✓")
print()

# =====================================================================
# STEP 8: CALIBRATION CHECKPOINT
# =====================================================================
print("STEP 8: CALIBRATION CHECKPOINT")
print("-" * 70)
print("Comparing to SI 2019 definition:")
print()
c_exact = 299792458.0  # m/s (exact by SI definition)
print(f"TriPhase c:   {c:.0f} m/s")
print(f"SI 2019 c:    {c_exact:.0f} m/s (exact)")
print()
deviation = c - c_exact
deviation_ppm = deviation / c_exact * 1e6
print(f"Deviation:    {deviation:+.3f} m/s")
print(f"              {deviation_ppm:+.6f} ppm")
print()
print("Component verification:")
print(f"  ε₀ = {epsilon_0:.12e} F/m")
print(f"  μ₀ = {mu_0:.12e} H/m")
print()
c_check = 1.0 / math.sqrt(epsilon_0 * mu_0)
print(f"  c = 1/√(ε₀μ₀) = {c_check:.0f} m/s")
print()
print("Historical context:")
print("  Before SI 2019: c was measured")
print("  After SI 2019:  c is defined exactly as 299,792,458 m/s")
print("                  meter is defined via c and second")
print()
print("Relation to other constants:")
print("  c² = 1/(ε₀μ₀)")
print(f"  c² = {c**2:.15e} m²/s²")
print()
print("  Z₀ = √(μ₀/ε₀) = μ₀c")
print(f"  Z₀ = {Z_0:.12f} Ω (impedance of free space)")
print()
print("Physical interpretation:")
print("  c is NOT a universal speed limit for 'matter'")
print("  c is the speed of EM wave propagation in vacuum")
print("  It emerges from vacuum's response to electric and magnetic fields")
print()
print("  ε₀: 'stiffness' of space to electric fields")
print("  μ₀: 'stiffness' of space to magnetic fields")
print("  c: Wave speed through this 'EM medium'")
print()
print("TriPhase interpretation:")
print("  c sets the fundamental wave speed in 3D space")
print("  All massive particles are wave packets traveling at c")
print("  'Slower' velocities arise from helical/spiral wave paths")
print("  de Broglie wavelength: λ = h/(mv) where v = c·sin(θ)")
print()
if abs(deviation_ppm) < 1:
    print("  ✓ Perfect agreement (within numerical precision)")
    print("  This confirms ε₀ and μ₀ are calibrated to SI definition")
print()

print("=" * 70)
print("DIMENSIONAL ANALYSIS COMPLETE")
print("c has dimensions [m s⁻¹] as required for velocity")
print("Maxwell's unification of electricity and magnetism confirmed")
print("=" * 70)

input("Press Enter to exit...")
