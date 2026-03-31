"""
TriPhase V16 — Z Boson Mass (Renormalization Group Framework)
========================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

RENORMALIZATION GROUP INTERPRETATION
--------------------------------------
The Z boson mass is related to the W mass via the Weinberg angle θ_W, which
encodes the mixing between U(1)_Y and SU(2)_L gauge groups at electroweak
symmetry breaking. The TriPhase formula M_Z = M_W/cos(θ_W) with sin²(θ_W) = απ
connects the Z mass to gauge coupling unification. The Weinberg angle is NOT
a free parameter but arises from the RG running of g₁ and g₂ to a common scale.

In the electroweak RG framework, the Weinberg angle runs with energy scale:
sin²θ_W(μ) = g'²/(g² + g'²), where g and g' are the SU(2)_L and U(1)_Y couplings.
At the Z pole, sin²θ_W ≈ 0.231, close to the TriPhase prediction απ ≈ 0.229.
This near-equality suggests that the fine structure constant α and the geometric
factor π are fundamental to gauge coupling unification.

The Z mass, together with the W mass and Higgs mass, determines the electroweak
vacuum stability. The RG flow of these masses toward the Planck scale reveals
whether the Higgs potential develops instabilities. TriPhase predicts M_Z through
the α¹⁸ cascade, suggesting that electroweak physics is embedded in the same
geometric structure that generates the cosmic Hubble scale.

TAG: (D) — Pure derivation from Weinberg angle and gauge unification
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
print("TriPhase V16: Z Boson Mass (Renormalization Group)")
print("=" * 70)
print()

print("GAUGE COUPLING UNIFICATION AND WEINBERG ANGLE")
print("-" * 70)
print(f"Fine structure constant:         α   = {alpha:.10f}")
print(f"Geometric factor:                π   = {math.pi:.10f}")
print(f"Weinberg angle prediction:       sin²(θ_W) = απ = {alpha * math.pi:.10f}")
print(f"Weinberg angle:                  θ_W = {math.degrees(math.asin(math.sqrt(alpha * math.pi))):.3f}°")
print()

print("W BOSON MASS (FROM PREVIOUS DERIVATION)")
print("-" * 70)
M_W = m_p * T_17 / (2 * alpha)
print(f"W boson mass:                    M_W = {M_W / 1.782662e-28 / 1000:.3f} GeV/c²")
print()

print("Z BOSON MASS FROM ELECTROWEAK MIXING")
print("-" * 70)
print("The Z boson mass is related to M_W via:")
print("  M_Z = M_W / cos(θ_W)")
print()
print("where the Weinberg angle θ_W encodes U(1)_Y × SU(2)_L mixing.")
print()
print("TriPhase predicts sin²(θ_W) = απ, connecting electroweak unification")
print("to the fundamental constants α and π.")
print()

sin2_theta_W = alpha * math.pi
cos_theta_W = math.sqrt(1 - sin2_theta_W)
M_Z = M_W / cos_theta_W

print(f"cos(θ_W):                        {cos_theta_W:.10f}")
print(f"Z boson mass (TriPhase):         M_Z = {M_Z:.6e} kg")
print(f"                                     = {M_Z / 1.782662e-28:.3f} MeV/c²")
print(f"                                     = {M_Z / 1.782662e-28 / 1000:.3f} GeV/c²")
print()

# ========== CALIBRATION CHECKPOINT ==========
M_Z_CODATA_GeV = 91.1876  # GeV/c² (PDG 2024, Z pole mass)
M_Z_CODATA = M_Z_CODATA_GeV * 1.782662e-25  # kg
sin2_theta_W_CODATA = 0.23121  # On-shell scheme (PDG 2024)
deviation_ppm = abs(M_Z - M_Z_CODATA) / M_Z_CODATA * 1e6
sin2_deviation_pct = abs(sin2_theta_W - sin2_theta_W_CODATA) / sin2_theta_W_CODATA * 100

print("CALIBRATION vs. CODATA/PDG")
print("-" * 70)
print(f"CODATA Z boson mass:             {M_Z_CODATA_GeV:.4f} GeV/c² (Z pole)")
print(f"TriPhase Z boson mass:           {M_Z / 1.782662e-28 / 1000:.3f} GeV/c²")
print(f"Deviation:                       {deviation_ppm:.0f} ppm ({abs(M_Z - M_Z_CODATA)/M_Z_CODATA * 100:.2f}%)")
print()
print(f"CODATA sin²(θ_W):                {sin2_theta_W_CODATA:.5f} (on-shell)")
print(f"TriPhase sin²(θ_W) = απ:         {sin2_theta_W:.5f}")
print(f"Weinberg angle deviation:        {sin2_deviation_pct:.2f}%")
print()

print("RENORMALIZATION GROUP INSIGHT")
print("-" * 70)
print("The Z mass emerges from electroweak gauge coupling unification. The")
print("Weinberg angle θ_W runs with energy scale, approaching the TriPhase value")
print("sin²(θ_W) = απ at high energies. This suggests that α and π are fundamental")
print("to gauge unification, not merely electromagnetic coupling constants. The RG")
print("flow M_Z(μ) connects electroweak breaking to the α¹⁸ cosmic cascade.")
print()
print("=" * 70)

input("Press Enter to exit...")
