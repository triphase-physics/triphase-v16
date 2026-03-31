"""
================================================================================
TriPhase V16 - Z Boson Mass Derivative
Framework: THERMODYNAMICS
Tag: (D*H) - Derived with hypothetical components
================================================================================

THERMODYNAMICS FRAMEWORK:
Interprets each physical quantity through statistical mechanics / thermodynamic
concepts: partition functions Z and free energy F = -k_BT ln(Z), entropy
S = -∂F/∂T, equipartition theorem (½k_BT per degree of freedom), Boltzmann
distributions, Maxwell-Boltzmann statistics, phase transitions, order parameters,
critical phenomena, equations of state, thermodynamic potentials (U, H, F, G),
degrees of freedom counting, mode counting, heat capacity, specific heat,
adiabatic processes, black-body radiation, Planck distribution, thermodynamic
stability conditions, Stefan-Boltzmann law, Wien displacement, chemical potential,
Gibbs free energy.

PHYSICAL DERIVATION:
The Z boson mass emerges from the same electroweak phase transition as the W,
but with a different mass due to the mixing of neutral gauge bosons. The Z is
a superposition of the W³ (SU(2)) and B (U(1)) gauge bosons:

    Z_μ = cos(θ_W) W³_μ - sin(θ_W) B_μ

The orthogonal combination gives the massless photon. The Z mass is related
to the W mass through the Weinberg angle:

    M_Z = M_W / cos(θ_W)

This is a thermodynamic constraint — the same Higgs VEV v generates both masses,
but the mixing angle determines their ratio.

The formula:
    M_Z = M_W / cos(θ_W)

is exact at tree level and receives small radiative corrections.

Thermodynamic interpretation: The Z mass, like the W mass, vanishes above T_EW
and emerges below T_EW through the same order parameter ⟨φ⟩ = v. The M_W/M_Z
ratio is a thermodynamic prediction independent of temperature.

Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# ============================================================================
# STANDARD ANCHOR CHAIN
# ============================================================================
epsilon_0 = 8.8541878128e-12   # F/m
mu_0      = 1.25663706212e-6   # H/m
e         = 1.602176634e-19    # C (exact SI)
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
T_17      = 17 * 18 // 2   # = 153
mp_me     = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p       = m_e * mp_me
H_0       = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r      = c**4 / (8.0 * math.pi * G)

# ============================================================================
# THERMODYNAMIC DERIVATION - Z BOSON MASS
# ============================================================================

print("=" * 80)
print("Z BOSON MASS - THERMODYNAMICS FRAMEWORK")
print("=" * 80)
print()
print("THERMODYNAMIC INTERPRETATION:")
print("The Z boson mass emerges from the same EW phase transition as M_W.")
print("The mass ratio M_Z/M_W = 1/cos(θ_W) is a thermodynamic constraint.")
print()

# Weinberg angle
sin2_theta_W = 0.23122
cos_theta_W = math.sqrt(1.0 - sin2_theta_W)
theta_W_deg = math.acos(cos_theta_W) * 180.0 / math.pi

# W boson mass (from previous derivation)
M_W = m_p * mp_me * alpha / (2.0 * sin2_theta_W)
M_W_GeV = M_W * c**2 / (1.602176634e-19 * 1e9)

# Z boson mass
M_Z = M_W / cos_theta_W
M_Z_GeV = M_Z * c**2 / (1.602176634e-19 * 1e9)

print("THERMODYNAMIC DERIVATION:")
print(f"  W boson mass M_W            = {M_W_GeV:.4f} GeV/c²")
print(f"  Weinberg angle θ_W          = {theta_W_deg:.2f}°")
print(f"  sin²θ_W                     = {sin2_theta_W:.5f}")
print(f"  cos θ_W                     = {cos_theta_W:.5f}")
print()
print(f"  M_Z = M_W / cos(θ_W)")
print(f"      = {M_W_GeV:.4f} / {cos_theta_W:.5f}")
print(f"      = {M_Z_GeV:.4f} GeV/c²")
print()

# Gauge boson mixing
print("GAUGE BOSON MIXING:")
print(f"  The neutral electroweak sector has 2 gauge bosons before SSB:")
print(f"    W³_μ from SU(2)_L (weak isospin)")
print(f"    B_μ from U(1)_Y (weak hypercharge)")
print()
print(f"  After symmetry breaking, they mix:")
print(f"    Z_μ = cos(θ_W) W³_μ - sin(θ_W) B_μ  (massive)")
print(f"    A_μ = sin(θ_W) W³_μ + cos(θ_W) B_μ  (massless photon)")
print()
print(f"  The mixing angle θ_W is determined by gauge couplings:")
print(f"    tan(θ_W) = g'/g where g' is U(1) coupling, g is SU(2) coupling")
print()

# Mass generation from Higgs VEV
v_GeV = 246.0
g_weak = 2.0 * M_W_GeV / v_GeV
g_prime = g_weak * math.tan(math.acos(cos_theta_W))

print("MASS GENERATION MECHANISM:")
print(f"  Higgs VEV v                 = {v_GeV:.1f} GeV")
print(f"  SU(2) coupling g            = {g_weak:.4f}")
print(f"  U(1) coupling g'            = {g_prime:.4f}")
print()
print(f"  M_W = (g/2) × v             = {M_W_GeV:.2f} GeV")
print(f"  M_Z = (√(g² + g'²)/2) × v   = {M_Z_GeV:.2f} GeV")
print(f"  M_γ = 0                     (U(1)_EM unbroken)")
print()

# Electroweak phase transition
k_B = 1.380649e-23  # J/K
T_EW_GeV = 160.0
T_EW_K = T_EW_GeV * 1e9 * 1.602176634e-19 / k_B

print("ELECTROWEAK PHASE TRANSITION:")
print(f"  Critical temperature T_EW   = {T_EW_GeV:.1f} GeV = {T_EW_K:.3e} K")
print()
print(f"  T > T_EW: M_W = M_Z = 0 (SU(2)×U(1) symmetric)")
print(f"  T < T_EW: M_W = {M_W_GeV:.1f} GeV, M_Z = {M_Z_GeV:.1f} GeV (broken to U(1)_EM)")
print()

# Free energy contribution
print("FREE ENERGY ANALYSIS:")
print(f"  The Z boson contributes to the finite-T effective potential:")
print(f"    V_eff ∋ (T⁴/π²) × n_B × J_B(M_Z²(φ)/T²)")
print(f"  where n_B = 3 (Z polarizations) and J_B is the bosonic integral")
print()
print(f"  At high T: Thermal mass dominates, M_Z,thermal ∝ gT")
print(f"  At low T: Vacuum mass dominates, M_Z = gv/cos(θ_W)")
print()

# Partition function
print("PARTITION FUNCTION:")
print(f"  Z boson thermal partition function:")
print(f"    Z_Z = ∏_k [1 - exp(-E_k/T)]^(-3)")
print(f"  where E_k = √(k² + M_Z²(T))")
print()
print(f"  The free energy F = -T ln(Z_Z) depends on M_Z(T)")
print(f"  Below T_EW, M_Z tracks the Higgs VEV: M_Z(T) ∝ v(T)")
print()

# Temperature-dependent mass
T_below_GeV = 100.0  # Example: T < T_EW
v_T = v_GeV * math.sqrt(1.0 - (T_below_GeV/T_EW_GeV)**2) if T_below_GeV < T_EW_GeV else 0
M_Z_T_GeV = (math.sqrt(g_weak**2 + g_prime**2) / 2.0) * v_T

print("TEMPERATURE-DEPENDENT MASS:")
print(f"  Order parameter v(T) ∝ √(1 - (T/T_EW)²) for T < T_EW")
print(f"  At T = {T_below_GeV:.0f} GeV:")
print(f"    v(T) = {v_T:.1f} GeV")
print(f"    M_Z(T) = {M_Z_T_GeV:.1f} GeV")
print(f"  As T → 0: M_Z(T) → M_Z = {M_Z_GeV:.1f} GeV")
print()

# Rho parameter
rho = (M_W_GeV / (M_Z_GeV * cos_theta_W))**2
print("RHO PARAMETER (CUSTODIAL SYM):")
print(f"  ρ ≡ M_W² / (M_Z² cos²θ_W)   = {rho:.6f}")
print(f"  ρ = 1 at tree level (custodial SU(2) symmetry)")
print(f"  Deviations from 1 signal new physics (e.g., Δρ ∝ m_t²)")
print()

# ============================================================================
# CALIBRATION COMPARISON
# ============================================================================

print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)

# PDG Z boson mass
M_Z_PDG = 91.1876  # GeV/c²
M_Z_PDG_unc = 0.0021  # GeV/c²

deviation = ((M_Z_GeV - M_Z_PDG) / M_Z_PDG) * 100

# Also check the mass ratio
ratio_derived = M_Z_GeV / M_W_GeV
ratio_expected = 1.0 / cos_theta_W
ratio_PDG = M_Z_PDG / 80.377

print()
print(f"TriPhase derived M_Z        = {M_Z_GeV:.4f} GeV/c²")
print(f"PDG measured M_Z            = {M_Z_PDG:.4f} ± {M_Z_PDG_unc:.4f} GeV/c²")
print()
print(f"Deviation from PDG          = {deviation:+.2f}%")
print()
print(f"Mass ratio M_Z/M_W:")
print(f"  TriPhase derived            = {ratio_derived:.5f}")
print(f"  Expected 1/cos(θ_W)         = {ratio_expected:.5f}")
print(f"  PDG measured                = {ratio_PDG:.5f}")
print()
print("NOTE: Z mass is known to 0.002% precision from LEP Z-pole scans.")
print("      M_Z/M_W ratio tests the electroweak theory structure.")
print()

print("THERMODYNAMIC SIGNIFICANCE:")
print("  The Z boson mass is a thermodynamic order parameter indicator:")
print("    - Vanishes in the symmetric phase (T > T_EW)")
print("    - Emerges through the same SSB as M_W")
print("    - The ratio M_Z/M_W is a thermodynamic constraint")
print("  Precision Z measurements test the EW phase transition mechanism.")
print()

print("=" * 80)
print("END Z BOSON MASS THERMODYNAMICS DERIVATIVE")
print("=" * 80)

input("Press Enter to exit...")
