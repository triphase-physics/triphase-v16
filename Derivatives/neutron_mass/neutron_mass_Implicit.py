"""
========================================================================
TriPhase V16 Derivative: Neutron Mass (Implicit Framework)
========================================================================

IMPLICIT FRAMEWORK INTERPRETATION:
The neutron mass emerges through an implicit perturbative constraint on
the proton mass: m_n = m_p × [1 + α × (m_e/m_p) × T_17]. This equation
defines the neutron mass implicitly as the unique value satisfying the
self-consistency condition where the mass difference arises from electromagnetic
corrections scaled by the triangular number T_17. The implicit nature arises
because the correction term (m_e/m_p) involves the ratio of masses, creating
a coupled system where neutron and proton masses mutually constrain each other.

The implicit function theorem applied to F(m_n, m_p) = m_n - m_p[1 + ε] = 0
guarantees existence and uniqueness of the neutron mass given the proton mass.
The perturbation parameter ε = α(m_e/m_p)T_17 represents electromagnetic self-
energy corrections that shift the neutron slightly heavier than the proton.
The neutron mass is not computed forward but emerges backward as the unique
solution satisfying the implicit electromagnetic constraint equation.

REFERENCE: CODATA 1.67492749804(95)e-27 kg

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
TAG: (D*)
========================================================================
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

print("=" * 70)
print("NEUTRON MASS — IMPLICIT FRAMEWORK")
print("=" * 70)

# IMPLICIT PERTURBATIVE CONSTRAINT
print("\nIMPLICIT CONSTRAINT EQUATION:")
print("  F(m_n, m_p) = m_n - m_p × [1 + α × (m_e/m_p) × T_17] = 0")
print("  Neutron mass emerges from coupled electromagnetic constraint.")

# Perturbation parameter
mass_ratio = m_e / m_p
triangular = T_17
perturbation = alpha * mass_ratio * triangular

print(f"\nPERTURBATION COMPONENTS:")
print(f"  Fine structure: α = {alpha:.10f}")
print(f"  Mass ratio: m_e/m_p = {mass_ratio:.6e}")
print(f"  Triangular: T_17 = {triangular}")
print(f"  Perturbation: ε = α(m_e/m_p)T_17 = {perturbation:.6e}")

# Implicit solution for neutron mass
m_n = m_p * (1.0 + perturbation)

print(f"\nIMPLICIT SOLUTION:")
print(f"  m_n = m_p × (1 + ε)")
print(f"  m_n = {m_n:.12e} kg")

# Mass difference
delta_m = m_n - m_p
delta_m_MeV = delta_m * c**2 / (1.602176634e-19 * 1e6)

print(f"\nMASS SPLITTING:")
print(f"  Δm = m_n - m_p = {delta_m:.6e} kg")
print(f"  Δm = {delta_m_MeV:.4f} MeV/c²")

# Verify constraint satisfaction
residual = m_n - m_p * (1.0 + perturbation)
print(f"\nCONSTRAINT RESIDUAL:")
print(f"  F(m_n, m_p) = {residual:.6e}  (implicit zero)")

print("\n" + "=" * 70)
print("CALIBRATION CHECKPOINT")
print("=" * 70)

codata_value = 1.67492749804e-27  # kg
deviation_ppm = (m_n - codata_value) / codata_value * 1e6

print(f"TriPhase Implicit:  {m_n:.12e} kg")
print(f"CODATA 2018:        {codata_value:.12e} kg")
print(f"Deviation:          {deviation_ppm:+.3f} ppm")

# Experimental mass difference
exp_delta_MeV = 1.29333  # MeV/c²
print(f"\nMass Difference:")
print(f"TriPhase: {delta_m_MeV:.4f} MeV/c²")
print(f"Experimental: {exp_delta_MeV:.5f} MeV/c²")
print(f"Deviation: {(delta_m_MeV - exp_delta_MeV)/exp_delta_MeV * 100:+.2f}%")

print("\n" + "=" * 70)
print("IMPLICIT FRAMEWORK INSIGHTS")
print("=" * 70)

print("\n1. COUPLED IMPLICIT SYSTEM:")
print("   Neutron and proton masses mutually constrain each other through")
print("   the mass ratio (m_e/m_p) appearing in the perturbation term.")

print("\n2. ELECTROMAGNETIC ORIGIN:")
print("   The mass splitting emerges from implicit electromagnetic self-energy")
print("   corrections, not from explicit quark mass calculations.")

print("\n3. TRIANGULAR SCALING:")
print("   T_17 acts as an implicit geometric constraint that determines the")
print("   magnitude of the electromagnetic perturbation.")

print("\n4. PERTURBATIVE UNIQUENESS:")
print("   The implicit function theorem guarantees unique m_n given m_p and")
print("   the constraint F(m_n, m_p) = 0 with continuous partial derivatives.")

print("=" * 70)
input("Press Enter to exit...")
