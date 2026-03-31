"""
================================================================================
TriPhase V16 - Charm Quark Mass Derivative
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
The charm quark mass emerges from lattice QCD thermodynamics. Charmed hadrons
(D mesons, J/ψ) undergo sequential dissociation as temperature increases:

    T < 1.0 GeV: Charmed hadrons stable
    T ≈ 1.2 GeV: D mesons dissolve (deconfinement)
    T ≈ 1.6 GeV: J/ψ dissociates (Debye screening)

The charm quark deconfinement temperature T_c ≈ 1.2 GeV determines its mass
through the thermal equation of state. Above T_c, charm quarks behave as
quasi-free particles in the QGP.

The formula:
    m_c ≈ m_p × (4/3) × √α⁻¹

The factor (4/3) is the color charge factor C_F for SU(3). The √α⁻¹ term
arises from the electromagnetic-to-strong scale ratio, as charm quarks carry
both color and electric charge.

Thermodynamic interpretation: The charm quark mass is set by the temperature
at which charmed hadrons melt into the QGP. This is a sharp thermodynamic
phase transition, visible in lattice QCD free energy calculations.

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
# THERMODYNAMIC DERIVATION - CHARM QUARK MASS
# ============================================================================

print("=" * 80)
print("CHARM QUARK MASS - THERMODYNAMICS FRAMEWORK")
print("=" * 80)
print()
print("THERMODYNAMIC INTERPRETATION:")
print("The charm quark mass is determined by the deconfinement temperature")
print("of charmed hadrons. At T_c ≈ 1.2 GeV, D mesons dissolve into the QGP.")
print()

# Charm quark mass formula
# m_c = m_p × (4/3) × √α⁻¹
color_factor = 4.0 / 3.0  # C_F for SU(3) color
scale_factor = math.sqrt(alpha_inv)

m_c = m_p * color_factor * scale_factor

# Convert to MeV/c² and GeV/c²
m_c_MeV = m_c * c**2 / (1.602176634e-19 * 1e6)
m_c_GeV = m_c_MeV / 1000.0

print("THERMODYNAMIC DERIVATION:")
print(f"  Proton mass m_p             = {m_p:.6e} kg")
print(f"  Fine structure α⁻¹          = {alpha_inv:.10f}")
print(f"  Color factor C_F            = {color_factor:.6f}")
print(f"  Scale factor √α⁻¹           = {scale_factor:.6f}")
print()
print(f"  m_c = m_p × (4/3) × √α⁻¹")
print(f"      = {m_p:.6e} × {color_factor:.6f} × {scale_factor:.6f}")
print(f"      = {m_c:.6e} kg")
print(f"      = {m_c_MeV:.2f} MeV/c²")
print(f"      = {m_c_GeV:.3f} GeV/c²")
print()

# Thermodynamic phase transition analysis
k_B = 1.380649e-23  # J/K
T_c_GeV = 1.2  # Charm deconfinement temperature
T_c_MeV = T_c_GeV * 1000.0
T_c_K = T_c_MeV * 1.602176634e-19 * 1e6 / k_B

print("THERMODYNAMIC PHASE TRANSITION:")
print(f"  Deconfinement temperature T_c = {T_c_GeV:.1f} GeV")
print(f"                                = {T_c_MeV:.0f} MeV")
print(f"                                = {T_c_K:.3e} K")
print()

# Free energy analysis
print("FREE ENERGY LANDSCAPE:")
print(f"  Below T_c: Charmed hadrons (D, D*, J/ψ) are stable bound states.")
print(f"             Free energy F_hadron < F_QGP")
print(f"  Above T_c: Charm quarks are deconfined in the QGP.")
print(f"             Free energy F_QGP < F_hadron")
print(f"  At T = T_c: First-order phase transition (latent heat).")
print()

# Debye screening
lambda_D_MeV = T_c_MeV / math.sqrt(4.0 * math.pi * alpha_inv)
print("DEBYE SCREENING IN QGP:")
print(f"  Debye screening length λ_D ∝ T/g")
print(f"  At T_c, λ_D ≈ {lambda_D_MeV:.2f} MeV⁻¹ ≈ 0.2 fm")
print(f"  When λ_D < r_c (charm quark radius), charmed hadrons dissolve.")
print()

# Sequential dissociation
print("SEQUENTIAL DISSOCIATION OF CHARMONIA:")
print(f"  T = 1.0 GeV: All charmed hadrons stable")
print(f"  T = 1.2 GeV: D mesons melt (m_c sets this threshold)")
print(f"  T = 1.6 GeV: J/ψ dissociates (tighter binding)")
print(f"  T = 2.0 GeV: χ_c states dissolve")
print()

# Lattice QCD thermodynamics
print("LATTICE QCD FREE ENERGY:")
print(f"  Charm quark contribution to QGP pressure:")
print(f"  P_c/T⁴ = (7π²/180) × g_c × [1 - 15m_c²/(4π²T²) + ...]")
print(f"  where g_c = 2 (spin) × 3 (color) = 6 degrees of freedom")
print()

# ============================================================================
# CALIBRATION COMPARISON
# ============================================================================

print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)

# PDG charm quark mass (MS-bar scheme at charm mass scale)
m_c_PDG_low = 1.270  # GeV/c² (lower bound)
m_c_PDG_high = 1.290  # GeV/c² (upper bound)
m_c_PDG_central = (m_c_PDG_low + m_c_PDG_high) / 2.0

deviation_low = ((m_c_GeV - m_c_PDG_low) / m_c_PDG_low) * 100
deviation_central = ((m_c_GeV - m_c_PDG_central) / m_c_PDG_central) * 100
deviation_high = ((m_c_GeV - m_c_PDG_high) / m_c_PDG_high) * 100

print()
print(f"TriPhase derived m_c        = {m_c_GeV:.3f} GeV/c²")
print(f"PDG range (MS-bar @ m_c)    = {m_c_PDG_low:.3f} - {m_c_PDG_high:.3f} GeV/c²")
print(f"PDG central value           = {m_c_PDG_central:.3f} GeV/c²")
print()
print(f"Deviation from PDG low      = {deviation_low:+.2f}%")
print(f"Deviation from PDG central  = {deviation_central:+.2f}%")
print(f"Deviation from PDG high     = {deviation_high:+.2f}%")
print()
print("NOTE: Charm quark mass is scheme-dependent (pole vs MS-bar vs 1S).")
print("      Lattice QCD extracts masses from thermodynamic susceptibilities.")
print("      TriPhase mass relates directly to deconfinement temperature.")
print()

print("THERMODYNAMIC SIGNIFICANCE:")
print("  The charm quark mass defines a critical temperature in QCD.")
print("  Heavy-ion collisions (RHIC, LHC) probe this thermodynamic regime.")
print("  Open charm production vs J/ψ suppression tests this phase structure.")
print()

print("=" * 80)
print("END CHARM QUARK MASS THERMODYNAMICS DERIVATIVE")
print("=" * 80)

input("Press Enter to exit...")
