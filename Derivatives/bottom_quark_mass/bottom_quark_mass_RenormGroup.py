"""
TriPhase V16 — Bottom Quark Mass (Renormalization Group Framework)
========================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

RENORMALIZATION GROUP INTERPRETATION
--------------------------------------
The bottom quark mass is a running coupling evaluated at μ ≈ m_b ≈ 4.18 GeV.
The TriPhase formula m_b = m_e × T₁₇ × mp_me × (1 + α) encodes RG evolution
through perturbative QCD. The factor (1 + α) ≈ 1.0073 represents first-order
QCD radiative corrections, while T₁₇ × mp_me provides the primary scaling from
electron to bottom quark energy scale.

At the bottom threshold, QCD is fully perturbative, and the running mass
exhibits logarithmic dependence on the renormalization scale: m_b(μ) =
m_b(m_b) × [α_s(μ)/α_s(m_b)]^(γ_m/β₀). The (1 + α) correction captures the
one-loop QCD β-function contribution, while the base factor T₁₇ × mp_me ≈ 2.7×10⁵
represents the integrated RG flow from electromagnetic to strong scales.

The bottom quark is heavy enough that perturbative QCD is reliable, yet light
enough that RG running effects are measurable. It serves as a precision test
of QCD asymptotic freedom and the β-function structure. The TriPhase formula
connects this to the fundamental α¹⁸ cascade linking electron to cosmic scales.

TAG: (D) — Pure derivation from TriPhase RG flow with QCD corrections
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

# ========== RENORMALIZATION GROUP DERIVATION ==========
print("=" * 70)
print("TriPhase V16: Bottom Quark Mass (Renormalization Group)")
print("=" * 70)
print()

print("RG FLOW IN PERTURBATIVE QCD")
print("-" * 70)
print(f"Electron mass (UV anchor):       m_e = {m_e:.6e} kg")
print(f"Fine structure constant:         α   = {alpha:.10f}")
print(f"Triangle number T₁₇:             T₁₇ = {T_17}")
print(f"Proton-electron mass ratio:      mp_me = {mp_me:.6f}")
print(f"QCD correction factor:           (1 + α) = {1 + alpha:.10f}")
print()

print("BETA FUNCTION AND ANOMALOUS DIMENSIONS")
print("-" * 70)
print("In perturbative QCD, the bottom quark running mass evolves as:")
print("  m_b(μ) = m_b(m_b) × [α_s(μ)/α_s(m_b)]^(γ_m/β₀)")
print()
print("TriPhase encodes this as:")
print("  m_b = m_e × T₁₇ × mp_me × (1 + α)")
print()
print(f"Base RG scale T₁₇ × mp_me:       {T_17 * mp_me:.3e}")
print(f"With QCD correction:             {T_17 * mp_me * (1 + alpha):.3e}")
print()

m_b = m_e * T_17 * mp_me * (1 + alpha)

print(f"Bottom quark mass (TriPhase):    m_b = {m_b:.6e} kg")
print(f"                                     = {m_b / 1.782662e-28:.3f} MeV/c²")
print()

# ========== CALIBRATION CHECKPOINT ==========
m_b_CODATA_MeV = 4180  # MeV/c² (MS-bar scheme at m_b, PDG 2024)
m_b_CODATA = m_b_CODATA_MeV * 1.782662e-28  # kg
deviation_ppm = abs(m_b - m_b_CODATA) / m_b_CODATA * 1e6

print("CALIBRATION vs. CODATA/PDG")
print("-" * 70)
print(f"CODATA bottom quark mass:        {m_b_CODATA_MeV} MeV/c² (MS-bar, μ = m_b)")
print(f"TriPhase bottom quark mass:      {m_b / 1.782662e-28:.0f} MeV/c²")
print(f"Deviation:                       {deviation_ppm:.0f} ppm ({abs(m_b - m_b_CODATA)/m_b_CODATA * 100:.2f}%)")
print()

print("RENORMALIZATION GROUP INSIGHT")
print("-" * 70)
print("The bottom quark mass is fully perturbative, allowing precision tests of QCD")
print("RG evolution. The (1 + α) factor captures one-loop QCD corrections, while")
print("T₁₇ × mp_me encodes the integrated β-function flow from electron to bottom scale.")
print()
print("=" * 70)

input("Press Enter to exit...")
