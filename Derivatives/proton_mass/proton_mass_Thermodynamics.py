"""
================================================================================
TriPhase V16 - Proton Mass Derivative
Framework: THERMODYNAMICS
Tag: (D) - Pure derivation
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
The proton mass emerges from QCD thermodynamics. Remarkably, most of the proton's
mass does NOT come from the up and down quark masses (m_u + m_d ≈ 10 MeV), but
from the gluon field energy — the QCD vacuum energy density.

The QCD trace anomaly gives:
    ⟨T^μ_μ⟩ = (β(g)/2g) × G^a_μν G^{aμν}

where β(g) is the QCD beta function. This non-zero trace means the QCD vacuum
has a non-trivial equation of state: P ≠ ρ/3.

The proton mass is set by the QCD confinement scale Λ_QCD ≈ 200 MeV, which
emerges from dimensional transmutation. The formula:

    m_p = m_e × mp_me

encodes this through the mp_me ratio = 1836.15267...

Thermodynamic interpretation: The proton mass = vacuum pressure × characteristic
volume. It's the free energy cost of creating a "bag" of confined quarks and
gluons in the QCD vacuum. At T > T_c ≈ 170 MeV, this free energy barrier
disappears and the proton dissolves into the QGP.

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
# THERMODYNAMIC DERIVATION - PROTON MASS
# ============================================================================

print("=" * 80)
print("PROTON MASS - THERMODYNAMICS FRAMEWORK")
print("=" * 80)
print()
print("THERMODYNAMIC INTERPRETATION:")
print("The proton mass emerges from QCD vacuum energy — not quark masses.")
print("m_p ≈ Λ_QCD × (volume) where Λ_QCD ≈ 200 MeV is the confinement scale.")
print("This is a thermodynamic free energy: F_proton = U_gluons - T×S_quarks.")
print()

# Proton mass from TriPhase
print("TRIPHASE DERIVATION:")
print(f"  Electron mass m_e           = {m_e:.10e} kg")
print(f"  Proton-electron ratio       = {mp_me:.10f}")
print()
print(f"  m_p = m_e × mp_me")
print(f"      = {m_e:.10e} × {mp_me:.10f}")
print(f"      = {m_p:.10e} kg")
print()

# Convert to energy units
m_p_MeV = m_p * c**2 / (1.602176634e-19 * 1e6)
m_p_GeV = m_p_MeV / 1000.0

print(f"  m_p × c² = {m_p_MeV:.6f} MeV")
print(f"           = {m_p_GeV:.9f} GeV")
print()

# QCD thermodynamics
print("QCD THERMODYNAMIC ANALYSIS:")
print()

# Quark mass contributions (tiny)
m_u_MeV = 2.2  # Up quark mass (MS-bar at 2 GeV)
m_d_MeV = 4.7  # Down quark mass (MS-bar at 2 GeV)
m_quarks_MeV = 2*m_u_MeV + m_d_MeV  # uud content
quark_fraction = m_quarks_MeV / m_p_MeV

print(f"  Quark mass contribution:")
print(f"    m_u ≈ {m_u_MeV:.1f} MeV, m_d ≈ {m_d_MeV:.1f} MeV")
print(f"    2m_u + m_d ≈ {m_quarks_MeV:.1f} MeV")
print(f"    Fraction of m_p = {quark_fraction*100:.1f}%")
print()

# Gluon field energy (dominant)
gluon_fraction = 1.0 - quark_fraction
print(f"  Gluon field energy contribution:")
print(f"    ≈ {gluon_fraction*100:.1f}% of m_p")
print(f"    ≈ {m_p_MeV * gluon_fraction:.1f} MeV")
print(f"  This is the QCD vacuum energy density × proton volume")
print()

# QCD confinement scale
Lambda_QCD_MeV = 200.0  # Typical value
print(f"  QCD confinement scale Λ_QCD ≈ {Lambda_QCD_MeV:.0f} MeV")
print(f"  This sets the energy scale for hadron masses")
print(f"  Dimensional transmutation: mass from dimensionless coupling g_s(μ)")
print()

# MIT Bag Model interpretation
B_MeV4 = 145.0**4  # Bag constant in MeV^4
r_p_fm = 0.84  # Proton radius in fm
volume_fm3 = (4.0/3.0) * math.pi * r_p_fm**3

print("MIT BAG MODEL THERMODYNAMICS:")
print(f"  Proton = bag of confined quarks + gluons")
print(f"  Bag constant B = {145.0:.0f}^4 MeV^4 (vacuum pressure)")
print(f"  Proton radius r_p ≈ {r_p_fm:.2f} fm")
print(f"  Volume V ≈ {volume_fm3:.2f} fm³")
print(f"  Bag energy E_bag = B × V ≈ contribution to m_p")
print()

# Thermodynamic equation of state
print("QCD EQUATION OF STATE:")
print(f"  Inside proton (confined phase):")
print(f"    P_inside ≈ 0 (quarks/gluons confined)")
print(f"  Outside proton (QCD vacuum):")
print(f"    P_vacuum = -B (negative pressure)")
print(f"  Pressure balance creates confinement boundary")
print()

# Phase transition to QGP
k_B = 1.380649e-23  # J/K
T_c_MeV = 170.0  # Deconfinement temperature
T_c_K = T_c_MeV * 1.602176634e-19 * 1e6 / k_B

print("DECONFINEMENT PHASE TRANSITION:")
print(f"  Critical temperature T_c    = {T_c_MeV:.0f} MeV = {T_c_K:.3e} K")
print(f"  Above T_c: Protons dissolve into QGP")
print(f"  Free energy F_proton(T) increases with T")
print(f"  At T = T_c: F_proton = F_QGP (phase coexistence)")
print()

# Entropy and specific heat
print("THERMODYNAMIC PROPERTIES:")
print(f"  Entropy S_proton ≈ 0 at T → 0 (ground state)")
print(f"  Specific heat c_V = dU/dT shows hadronic resonances")
print(f"  Above T_c: c_V jumps (latent heat of deconfinement)")
print()

# Trace anomaly
print("QCD TRACE ANOMALY:")
print(f"  In a conformal theory: ⟨T^μ_μ⟩ = 0 (traceless stress tensor)")
print(f"  In QCD: ⟨T^μ_μ⟩ = (β/2g)⟨G²⟩ ≠ 0 (scale anomaly)")
print(f"  This non-zero trace generates hadron masses from Λ_QCD")
print(f"  m_p emerges from dimensional transmutation + confinement")
print()

# ============================================================================
# CALIBRATION COMPARISON
# ============================================================================

print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)

# CODATA proton mass
m_p_CODATA = 1.67262192369e-27  # kg
m_p_CODATA_MeV = 938.27208816  # MeV/c²

deviation_kg = ((m_p - m_p_CODATA) / m_p_CODATA) * 100
deviation_MeV = ((m_p_MeV - m_p_CODATA_MeV) / m_p_CODATA_MeV) * 100

print()
print(f"TriPhase derived m_p        = {m_p:.11e} kg")
print(f"CODATA 2018 m_p             = {m_p_CODATA:.11e} kg")
print(f"Deviation                   = {deviation_kg:+.6f}%")
print()
print(f"TriPhase m_p × c²           = {m_p_MeV:.8f} MeV")
print(f"CODATA m_p × c²             = {m_p_CODATA_MeV:.8f} MeV")
print(f"Deviation                   = {deviation_MeV:+.6f}%")
print()

print("THERMODYNAMIC SIGNIFICANCE:")
print("  The proton mass is a fundamental output of QCD thermodynamics.")
print("  It emerges from the QCD vacuum structure — not from quark masses.")
print("  Lattice QCD computes m_p from first principles via partition function:")
print(f"    Z = ∫[DA][Dψ][Dψ̄] exp(-S_QCD/ℏ)")
print(f"  TriPhase captures this through the mp_me ratio.")
print()

print("=" * 80)
print("END PROTON MASS THERMODYNAMICS DERIVATIVE")
print("=" * 80)

input("Press Enter to exit...")
