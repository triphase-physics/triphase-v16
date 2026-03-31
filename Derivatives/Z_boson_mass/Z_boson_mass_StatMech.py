"""
TriPhase V16 — Z Boson Mass (Statistical Mechanics Framework)
==============================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
The Z boson mass (~91.1876 GeV) emerges from the same electroweak symmetry breaking
mechanism as the W boson, but with a crucial difference: the Z is a mixture of the
neutral SU(2)_L gauge boson (W⁰) and the U(1)_Y hypercharge boson (B). The mass
relationship m_Z = m_W / cos(θ_W) connects the Z and W masses through the weak mixing
angle θ_W, also known as the Weinberg angle. In statistical mechanics terms, the
partition function for the electroweak sector must respect gauge invariance, which
constrains the relative masses of W and Z bosons. The mixing angle θ_W ≈ 28.7° is
itself a statistical observable — it parameterizes the relative coupling strengths
g_2 (SU(2)_L) and g_1 (U(1)_Y) through tan(θ_W) = g_1 / g_2.

At finite temperature above T_EW ~ 159 GeV, both W and Z bosons become massless,
and the partition function exhibits restored SU(2)_L × U(1)_Y symmetry. Below T_EW,
the Higgs condensate ⟨φ⟩ = v/√2 gives masses m_W = g_2 v / 2 and m_Z = v√(g_1² + g_2²) / 2.
The Z mass is larger because it couples to both gauge groups. The precise measurement
of m_Z at LEP (σ_m_Z ~ 2.1 MeV) makes it one of the best-determined fundamental
constants, providing a stringent test of electroweak theory and constraining the
Higgs mass through radiative corrections in the partition function.

TAG: (D*) — TriPhase derivation requiring phenomenological coefficient
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

# ========== STATISTICAL MECHANICS DERIVATION ==========
print("=" * 70)
print("TriPhase V16: Z Boson Mass (Statistical Mechanics)")
print("=" * 70)
print()

print("STATISTICAL MECHANICS FRAMEWORK")
print("--------------------------------")
print("Ensemble: Canonical (electroweak broken phase)")
print("Gauge mixing: Z = W⁰ cos(θ_W) - B sin(θ_W)")
print("Mass relation: m_Z = m_W / cos(θ_W)")
print("Observable: Weinberg angle θ_W from coupling ratio")
print()

print("ELECTROWEAK GAUGE BOSON SPECTRUM")
print("---------------------------------")
print(f"Proton mass m_p = {m_p:.6e} kg")
print(f"Fine structure α = {alpha:.10f}")
print()

# TriPhase formula for Z boson mass
# First compute W mass, then apply m_Z = m_W / cos(θ_W)
coeff_W = 85.7  # W-to-proton ratio from previous derivation
m_W_tph = m_p * coeff_W

# Weak mixing angle: sin²(θ_W) ≈ 0.2312 (PDG)
sin2_theta_W = 0.23122
cos_theta_W = math.sqrt(1.0 - sin2_theta_W)
theta_W_deg = math.asin(math.sqrt(sin2_theta_W)) * 180.0 / math.pi

m_Z_tph = m_W_tph / cos_theta_W

print(f"W boson mass m_W = {m_W_tph:.6e} kg")
print(f"Weak mixing angle θ_W = {theta_W_deg:.4f}°")
print(f"sin²(θ_W) = {sin2_theta_W:.5f}")
print(f"cos(θ_W) = {cos_theta_W:.5f}")
print()
print(f"m_Z = m_W / cos(θ_W)")
print(f"m_Z (TriPhase) = {m_Z_tph:.6e} kg")
print()

# Convert to GeV/c²
m_Z_GeV = m_Z_tph * c**2 / (1.602176634e-10)
m_W_GeV = m_W_tph * c**2 / (1.602176634e-10)

print(f"m_W (TriPhase) = {m_W_GeV:.6f} GeV/c²")
print(f"m_Z (TriPhase) = {m_Z_GeV:.6f} GeV/c²")
print(f"m_Z / m_W ratio = {m_Z_GeV / m_W_GeV:.6f}")
print()

# ========== CALIBRATION CHECKPOINT ==========
print("CALIBRATION CHECKPOINT")
print("----------------------")
# PDG 2022: m_Z = 91.1876 ± 0.0021 GeV (one of best-known constants)
m_Z_pdg = 91.1876  # GeV/c²
m_W_pdg = 80.377   # GeV/c²
ratio_pdg = m_Z_pdg / m_W_pdg

deviation_Z = (m_Z_GeV - m_Z_pdg) / m_Z_pdg * 1e6
deviation_ratio = (m_Z_GeV / m_W_GeV - ratio_pdg) / ratio_pdg * 1e6

print(f"PDG 2022 m_Z:           {m_Z_pdg:.6f} GeV/c²")
print(f"TriPhase m_Z:           {m_Z_GeV:.6f} GeV/c²")
print(f"Deviation:              {deviation_Z:.0f} ppm")
print()
print(f"PDG m_Z/m_W ratio:      {ratio_pdg:.6f}")
print(f"TriPhase m_Z/m_W:       {m_Z_GeV / m_W_GeV:.6f}")
print(f"Ratio deviation:        {deviation_ratio:.0f} ppm")
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT")
print("-----------------------------")
print("The Z boson mass exemplifies how gauge symmetry constraints propagate through")
print("the partition function. When SU(2)_L × U(1)_Y breaks to U(1)_EM, three degrees")
print("of freedom are 'eaten' by W⁺, W⁻, and Z bosons (longitudinal polarizations),")
print("while the photon remains massless. The mass eigenstates (γ, Z) are mixtures")
print("of the interaction eigenstates (W⁰, B), with mixing angle θ_W determined by")
print("the ratio of coupling constants: tan(θ_W) = g_1 / g_2.")
print()
print("In statistical mechanics, this mixing arises from diagonalizing the mass matrix")
print("in the partition function. The Higgs VEV ⟨φ⟩ couples to both W⁰ and B through")
print("the covariant derivative, generating a 2×2 mass matrix. Diagonalization yields:")
print("  • Photon A: massless (m_γ = 0)")
print("  • Z boson: massive (m_Z = v√(g_1² + g_2²) / 2)")
print()
print("The relation m_Z = m_W / cos(θ_W) is a consequence of gauge invariance — it's")
print("not a free parameter but a prediction. The precise LEP measurement of m_Z to")
print("2.1 MeV (23 ppm) tests this prediction with extraordinary precision. Radiative")
print("corrections from virtual top quarks and Higgs bosons shift m_Z through loop")
print("diagrams, encoding quantum fluctuations in the partition function. These")
print("corrections allowed theorists to predict m_t ≈ 175 GeV and m_H ≈ 125 GeV")
print("before their direct discovery — a triumph of statistical field theory.")
print()
print("TriPhase connects Z mass to hadronic scale: m_Z/m_p ≈ 97.2, encoding the")
print("hierarchy between QCD confinement (Λ_QCD ~ 200 MeV) and electroweak breaking")
print("(v ~ 246 GeV). Both scales emerge from dimensional transmutation, yet their")
print("ratio remains unexplained in the Standard Model — a hint of deeper structure.")
print()
print("=" * 70)

input("Press Enter to exit...")
