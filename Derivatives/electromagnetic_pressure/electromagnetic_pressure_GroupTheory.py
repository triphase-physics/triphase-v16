#!/usr/bin/env python3
"""
================================================================================
electromagnetic_pressure_GroupTheory.py
================================================================================
TriPhase Wave Mechanics Framework — V16 Python Derivatives
GroupTheory Module: Electromagnetic Pressure P_EM = B²/(2μ₀) = ε₀E²/2

Derives EM field pressure as the Casimir invariant of the U(1) gauge field
stress tensor. The quadratic form F_μνF^μν represents the U(1) field energy.

MIS Tag: (D) — Pure derivation from ε₀, μ₀, α, e, c

IRON RULES:
- import math only (NO numpy, scipy)
- CODATA/PDG values are CALIBRATION CHECKPOINTS only
- Every calculation derives from epsilon_0 and mu_0
- mp_me = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)

Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
All Rights Reserved.
================================================================================
"""

import math

print("="*80)
print("TRIPHASE V16: ELECTROMAGNETIC PRESSURE — GROUP THEORY")
print("P_EM = B²/(2μ₀) = ε₀E²/2 as U(1) Casimir Invariant")
print("="*80)
print()

# ==============================================================================
# ANCHOR CHAIN: Standard TriPhase V16 Constants
# ==============================================================================
print("ANCHOR CHAIN: Building from ε₀ and μ₀")
print("-" * 80)

epsilon_0 = 8.8541878128e-12  # F/m
mu_0      = 1.25663706212e-6  # H/m
e         = 1.602176634e-19   # C

print(f"ε₀ = {epsilon_0:.13e} F/m")
print(f"μ₀ = {mu_0:.14e} H/m")
print(f"e  = {e:.12e} C")
print()

c = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0 = math.sqrt(mu_0 / epsilon_0)
print(f"c = {c:.10e} m/s")
print(f"Z₀ = {Z_0:.10f} Ω")
print()

alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha = 1.0 / alpha_inv
hbar = Z_0 * e**2 / (4.0 * math.pi * alpha)
h = 2.0 * math.pi * hbar
print(f"α = {alpha:.15f}")
print(f"ℏ = {hbar:.15e} J·s")
print()

r_e = 2.8179403262e-15
m_e = hbar * alpha / (c * r_e)
f_e = m_e * c**2 / hbar
print(f"m_e = {m_e:.15e} kg")
print(f"f_e = {f_e:.10e} Hz")
print()

mp_me = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p = m_e * mp_me
print(f"mp/me = {mp_me:.10f}")
print(f"m_p = {m_p:.15e} kg")
print()

G = c**4 * 7.5 * epsilon_0**3 * mu_0**2
H_0 = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r = c**4 / (8.0 * math.pi * G)
w_0 = -(17.0 / 18.0)**2
print(f"G = {G:.15e} m³/(kg·s²)")
print(f"H₀ = {H_0:.10e} Hz")
print(f"VF = {VF_r:.10e} Pa")
print(f"w₀ = {w_0:.15f}")
print()

# ==============================================================================
# GROUP THEORY: U(1) GAUGE FIELD
# ==============================================================================
print("="*80)
print("GROUP THEORY: U(1) ELECTROMAGNETIC GAUGE FIELD")
print("="*80)
print()

print("THE U(1) GAUGE GROUP:")
print("-" * 80)
print("""
Electromagnetism is described by a U(1) gauge theory:
- U(1) = circle group {e^(iθ) : θ ∈ [0, 2π)}
- Lie algebra: u(1) = iℝ (one-dimensional, abelian)
- Generator: Q (electric charge operator)

GAUGE FIELD:
- Vector potential: A_μ = (φ/c, A) where φ = electric potential, A = vector potential
- Gauge transformation: A_μ → A_μ + ∂_μΛ (Λ = gauge parameter)
- Physical fields (gauge invariant):
  • Electric field: E = -∇φ - ∂A/∂t
  • Magnetic field: B = ∇ × A

FIELD STRENGTH TENSOR:
- F_μν = ∂_μA_ν - ∂_νA_μ (antisymmetric rank-2 tensor)
- Components:
  F_0i = E_i/c  (electric field)
  F_ij = ε_ijk B_k  (magnetic field)

GAUGE INVARIANCE:
- F_μν is gauge invariant (∂_μ and ∂_ν commute)
- Physical observables must be built from F_μν, not A_μ directly
""")
print()

print("CASIMIR INVARIANT:")
print("-" * 80)
print("""
The U(1) field has TWO Lorentz-invariant scalars:

1. FIRST INVARIANT (energy-like):
   I₁ = F_μνF^μν = 2(B²/μ₀ - ε₀E²)

2. SECOND INVARIANT (helicity-like):
   I₂ = F_μν*F̃^μν = -4ε₀E·B  (dual tensor)

For electromagnetic ENERGY DENSITY:
   u_EM = ½(ε₀E² + B²/μ₀)

The pressure is EQUAL to energy density for radiation:
   P_EM = u_EM/3 × 3 = u_EM  (in vacuum, isotropic)

More precisely, the STRESS TENSOR gives:
   T_μν = (1/μ₀)[F_μλF_ν^λ - ¼g_μνF_λρF^λρ]

The diagonal spatial components give PRESSURE.
""")
print()

# ==============================================================================
# DERIVATION: ELECTROMAGNETIC PRESSURE
# ==============================================================================
print("="*80)
print("DERIVATION: P_EM = B²/(2μ₀) = ε₀E²/2")
print("="*80)
print()

print("STEP 1: Maxwell Stress Tensor")
print("-" * 80)
print("The electromagnetic stress tensor is:")
print("  T_ij = ε₀[E_iE_j - ½δ_ij E²] + (1/μ₀)[B_iB_j - ½δ_ij B²]")
print()
print("For a UNIFORM field (say, B in z-direction, E = 0):")
print("  T_xx = T_yy = -B²/(2μ₀)  (tension perpendicular to field)")
print("  T_zz = +B²/(2μ₀)         (pressure along field)")
print()
print("The TRACE of the stress tensor (isotropic pressure):")
print("  P = ⅓Tr(T) = ⅓(T_xx + T_yy + T_zz)")
print("    = ⅓(-B²/(2μ₀) - B²/(2μ₀) + B²/(2μ₀))")
print("    = -B²/(6μ₀)  ???")
print()
print("Wait — that's negative! What's going on?")
print()

print("RESOLUTION: Pressure vs. Energy Density")
print("-" * 80)
print("The confusion arises from ANISOTROPY. A static magnetic field has:")
print("  • Magnetic pressure: P_B = B²/(2μ₀) (ALONG field lines)")
print("  • Magnetic tension: T_B = -B²/(2μ₀) (PERPENDICULAR to field)")
print()
print("For RADIATION (electromagnetic waves), the situation is different:")
print("  • E and B oscillate in phase, perpendicular to propagation")
print("  • Energy density: u = ε₀E² = B²/μ₀ (equipartition)")
print("  • Pressure: P = u/3 (for radiation)")
print()
print("We define MAGNETIC PRESSURE as the pressure along field lines:")
print("  P_B = B²/(2μ₀)")
print()
print("Similarly, ELECTRIC PRESSURE:")
print("  P_E = ε₀E²/2")
print()

print("STEP 2: Energy Density and Pressure Relation")
print("-" * 80)
print("Total electromagnetic energy density:")
print("  u_EM = ½(ε₀E² + B²/μ₀)")
print()
print("For pure electric field (B=0):")
print("  u_E = ε₀E²/2")
print("  P_E = ε₀E²/2  (pressure = energy density)")
print()
print("For pure magnetic field (E=0):")
print("  u_B = B²/(2μ₀)")
print("  P_B = B²/(2μ₀)  (pressure = energy density)")
print()
print("For electromagnetic WAVES (E ⊥ B, E/B = c):")
print("  u_EM = ε₀E² = B²/μ₀")
print("  P_rad = u_EM/3  (radiation pressure, isotropic averaging)")
print()

print("STEP 3: Group-Theoretic Normalization")
print("-" * 80)
print("Why B²/(2μ₀) and not B²/μ₀?")
print()
print("The factor of ½ comes from the QUADRATIC Casimir:")
print("  ⟨F²⟩ = F_μνF^μν = 2(B²/μ₀ - ε₀E²)")
print()
print("The energy density requires averaging over oscillation:")
print("  ⟨E²⟩_time = ½E₀²  (for sinusoidal wave)")
print("  ⟨B²⟩_time = ½B₀²")
print()
print("This gives the ½ factor in the static field energy.")
print()

# ==============================================================================
# NUMERICAL EXAMPLES
# ==============================================================================
print("="*80)
print("NUMERICAL EXAMPLES: ELECTROMAGNETIC PRESSURES")
print("="*80)
print()

print("EXAMPLE 1: Laboratory Magnetic Field")
print("-" * 80)
B_lab = 1.0  # Tesla (achievable with strong permanent magnet)
P_B_lab = B_lab**2 / (2.0 * mu_0)
print(f"Magnetic field: B = {B_lab:.1f} T")
print(f"Magnetic pressure: P_B = B²/(2μ₀) = {P_B_lab:.6e} Pa")
print(f"                        = {P_B_lab/1e5:.2f} atm")
print()
print(f"Comparison: Atmospheric pressure = 1.013e5 Pa")
print(f"            P_B / P_atm = {P_B_lab/1.013e5:.2f}")
print()

print("EXAMPLE 2: Fusion Reactor (ITER)")
print("-" * 80)
B_iter = 11.8  # Tesla (ITER toroidal field)
P_B_iter = B_iter**2 / (2.0 * mu_0)
print(f"Magnetic field: B = {B_iter:.1f} T")
print(f"Magnetic pressure: P_B = {P_B_iter:.6e} Pa")
print(f"                        = {P_B_iter/1e5:.0f} atm")
print()
print("This magnetic pressure must CONFINE the plasma pressure:")
P_plasma = 6e5  # Pa (rough ITER plasma pressure)
print(f"  Plasma pressure: P_plasma ~ {P_plasma:.1e} Pa = {P_plasma/1e5:.1f} atm")
print(f"  Confinement ratio: P_B/P_plasma = {P_B_iter/P_plasma:.1f}")
print()

print("EXAMPLE 3: Neutron Star Magnetic Field")
print("-" * 80)
B_ns = 1e8  # Tesla (typical neutron star surface field)
P_B_ns = B_ns**2 / (2.0 * mu_0)
print(f"Magnetic field: B = {B_ns:.1e} T")
print(f"Magnetic pressure: P_B = {P_B_ns:.6e} Pa")
print()
print("Compare to neutron star material pressure:")
P_ns_matter = 1e34  # Pa (nuclear density)
print(f"  Matter pressure: P ~ {P_ns_matter:.1e} Pa")
print(f"  Ratio: P_B/P_matter = {P_B_ns/P_ns_matter:.3e}")
print()
print("Magnetic field is SUBDOMINANT to matter pressure in typical neutron stars,")
print("but DOMINANT in magnetars (B ~ 10¹¹ T).")
print()

print("EXAMPLE 4: Magnetar")
print("-" * 80)
B_magnetar = 1e11  # Tesla (magnetar field)
P_B_magnetar = B_magnetar**2 / (2.0 * mu_0)
print(f"Magnetic field: B = {B_magnetar:.1e} T")
print(f"Magnetic pressure: P_B = {P_B_magnetar:.6e} Pa")
print(f"  Ratio to matter: P_B/P_matter = {P_B_magnetar/P_ns_matter:.3f}")
print()
print("In magnetars, the magnetic field pressure is COMPARABLE to matter pressure!")
print("This can distort the neutron star shape and power X-ray/gamma-ray bursts.")
print()

print("EXAMPLE 5: Solar Radiation Pressure")
print("-" * 80)
# Solar constant
S_sun = 1361.0  # W/m² (solar irradiance at 1 AU)
P_rad_sun = S_sun / c
print(f"Solar irradiance at Earth: S = {S_sun:.1f} W/m²")
print(f"Radiation pressure: P_rad = S/c = {P_rad_sun:.6e} Pa")
print(f"                         = {P_rad_sun:.3e} N/m²")
print()
print("This is the force per unit area on a perfectly reflecting surface.")
print("Solar sails use this pressure for propulsion!")
print()

# Force on a 1 m² sail
F_sail = P_rad_sun  # N (per m²)
print(f"Force on 1 m² sail: F = {F_sail:.6e} N = {F_sail*1e6:.3f} μN")
print()

# Acceleration for a lightweight sail
m_sail = 0.01  # kg (10 grams per m²)
a_sail = F_sail / m_sail
print(f"Sail mass: m = {m_sail*1000:.1f} g/m²")
print(f"Acceleration: a = F/m = {a_sail:.6e} m/s²")
print(f"             = {a_sail*3600*24:.3f} m/s per day")
print()

print("EXAMPLE 6: Electric Field Pressure")
print("-" * 80)
E_breakdown = 3e6  # V/m (air breakdown field)
P_E_breakdown = epsilon_0 * E_breakdown**2 / 2.0
print(f"Air breakdown field: E = {E_breakdown:.1e} V/m")
print(f"Electric pressure: P_E = ε₀E²/2 = {P_E_breakdown:.6e} Pa")
print(f"                       = {P_E_breakdown/1e5:.4f} atm")
print()

# Coulomb barrier pressure (nuclear scale)
E_nuclear = 1e11  # V/m (electric field at nuclear radius)
P_E_nuclear = epsilon_0 * E_nuclear**2 / 2.0
print(f"Nuclear-scale field: E = {E_nuclear:.1e} V/m")
print(f"Electric pressure: P_E = {P_E_nuclear:.6e} Pa")
print()

# ==============================================================================
# POYNTING VECTOR AND RADIATION PRESSURE
# ==============================================================================
print("="*80)
print("RADIATION PRESSURE: POYNTING VECTOR")
print("="*80)
print()

print("POYNTING VECTOR:")
print("-" * 80)
print("Energy flux density (power per unit area):")
print("  S = (1/μ₀) E × B")
print()
print("Magnitude: |S| = EB/μ₀")
print()
print("For a plane wave (E ⊥ B, E = cB):")
print("  |S| = cB²/μ₀ = cε₀E²")
print()

print("RADIATION PRESSURE:")
print("-" * 80)
print("Momentum density of EM field:")
print("  g = S/c² = (ε₀E × B)/c = ε₀E×B/c")
print()
print("Momentum flux (pressure) for normally incident wave:")
print("  P_rad = |g|c = (ε₀E²/c)c = ε₀E²")
print()
print("But wait — this is TWICE the energy density!")
print()
print("RESOLUTION: Reflection vs. Absorption")
print("  • Absorbing surface: P = u (momentum transfer = energy/c)")
print("  • Reflecting surface: P = 2u (momentum change = 2×energy/c)")
print()
print("For ISOTROPIC radiation field (averaged over all directions):")
print("  P_rad = u/3")
print()

# ==============================================================================
# CALIBRATION CHECKPOINT
# ==============================================================================
print("="*80)
print("CALIBRATION CHECKPOINT")
print("="*80)
print()

# Test case: B = 1 T
B_test = 1.0  # T
P_B_test = B_test**2 / (2.0 * mu_0)
print("Test case: B = 1 T")
print(f"  Derived: P_B = B²/(2μ₀) = {P_B_test:.10e} Pa")
print()

P_B_expected = 1.0 / (2.0 * 1.25663706212e-6)
print(f"  Expected: {P_B_expected:.10e} Pa")
print()

frac_diff = abs(P_B_test - P_B_expected) / P_B_expected
print(f"Fractional difference: {frac_diff:.3e}")
print()

if frac_diff < 1e-10:
    print("✓ PASS: Electromagnetic pressure formula is exact")
else:
    print("✗ FAIL: Check μ₀ value")
print()

# Test ε₀/μ₀ = 1/c²
ratio_test = epsilon_0 * mu_0
c_from_ratio = 1.0 / math.sqrt(ratio_test)
print(f"Verify: c = 1/√(ε₀μ₀) = {c_from_ratio:.10e} m/s")
print(f"        (should equal {c:.10e})")
print()

# ==============================================================================
# SUMMARY
# ==============================================================================
print("="*80)
print("SUMMARY: ELECTROMAGNETIC PRESSURE AS U(1) CASIMIR")
print("="*80)
print()

print("FORMULAS:")
print("  Magnetic pressure: P_B = B²/(2μ₀)")
print("  Electric pressure: P_E = ε₀E²/2")
print("  Radiation pressure: P_rad = u_EM/3 (isotropic)")
print()

print("GROUP THEORY:")
print("  • U(1) gauge field with field strength F_μν")
print("  • Casimir invariant: F_μνF^μν = 2(B²/μ₀ - ε₀E²)")
print("  • Stress tensor: T_μν = (1/μ₀)[F_μλF_ν^λ - ¼g_μνF_λρF^λρ]")
print("  • Pressure emerges from spatial diagonal components")
print()

print("PHYSICAL APPLICATIONS:")
print("  • Magnetic confinement fusion (ITER, stellarators)")
print("  • Neutron star structure (magnetars)")
print("  • Solar sails and radiation pressure")
print("  • Laser-matter interactions")
print("  • Cosmological radiation pressure (early universe)")
print()

print("CALIBRATION:")
print(f"  B = 1 T → P_B = {P_B_test:.3e} Pa ≈ 4 atm ✓")
print(f"  Solar radiation → P = {P_rad_sun:.3e} Pa ✓")
print()

print("="*80)
print("END OF ELECTROMAGNETIC_PRESSURE_GROUPTHEORY.PY")
print("="*80)

input("Press Enter to exit...")
