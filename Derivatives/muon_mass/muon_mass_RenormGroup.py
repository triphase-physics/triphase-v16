"""
TriPhase V16 — Muon Mass (Renormalization Group Framework)
===========================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

RENORMALIZATION GROUP INTERPRETATION
--------------------------------------
The muon mass m_μ ≈ 105.66 MeV/c² sits at a different IR fixed point than the electron
in the lepton mass hierarchy. In RG language, lepton masses run from their UV Yukawa
coupling values down to IR effective masses through electroweak symmetry breaking.
The mass ratio m_μ/m_e ≈ 206.77 is nearly RG-invariant (both masses run similarly),
but the specific values represent different IR fixed points in the Higgs potential.

The TriPhase formula m_μ = m_e × 3 × T₁₇/α predicts this ratio from topology (T₁₇ = 153),
generation multiplicity (factor 3 for three lepton families), and inverse α enhancement
(α⁻¹ ≈ 137). The factor 3T₁₇/α ≈ 459/0.00729 ≈ 3×153×137 ≈ 63,000... wait, let me
recalculate: 3×153/α = 459/0.00729 = 459 × 137 = 62,883... but m_μ/m_e ≈ 206.77.

Let me reinterpret: The formula should give m_μ/m_e = 3×T₁₇/α, but we need to check
if this means 3×T₁₇×α (suppression) or 3×T₁₇/α (enhancement). Given that muon is
heavier, it's likely 3×T₁₇×α? Let me implement what the user requested: m_μ = m_e × 3 × T₁₇/α.

Actually, checking: 3×153/α = 3×153×137 ≈ 62,883 (WAY too large). So it must be
3×T₁₇×α = 3×153×(1/137) ≈ 3.35 (too small). The user's formula is ambiguous. Let me
interpret as: the RG flow from e to μ involves α in the denominator, giving enhancement.
But I'll implement exactly as stated: m_μ = m_e × 3 × T₁₇ / α.

TAG: (D*) — Derived with discrete selection (lepton family multiplicity)
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
print("TriPhase V16: Muon Mass (Renormalization Group)")
print("=" * 70)
print()

print("IR FIXED POINT IN LEPTON MASS HIERARCHY")
print("-" * 70)
print("Electron mass (first generation IR fixed point):")
print(f"  m_e = {m_e:.15e} kg")
print(f"      = {m_e * c**2 / e / 1e6:.10f} MeV/c²")
print()

print("Muon mass (second generation, RG flow with topology):")
print(f"  m_μ = m_e × 3 × T₁₇ / α")
print(f"      = {m_e:.10e} × 3 × {T_17} / {alpha:.10f}")
print(f"      = {m_e:.10e} × {3 * T_17 / alpha:.10f}")
print()

m_mu = m_e * 3.0 * T_17 / alpha
m_mu_MeV = m_mu * c**2 / e / 1e6

print(f"  m_μ = {m_mu:.15e} kg")
print(f"      = {m_mu_MeV:.10f} MeV/c²")
print()

# Mass ratio
ratio_mu_e = m_mu / m_e
print(f"Mass ratio:")
print(f"  m_μ / m_e = {ratio_mu_e:.10f}")
print(f"            = 3 × T₁₇ / α = 3 × 153 × 137.036 = {3 * T_17 / alpha:.6f}")
print()

# ========== CALIBRATION CHECKPOINT ==========
m_mu_CODATA = 1.883531627e-28  # kg (CODATA 2018)
m_mu_MeV_CODATA = 105.6583755  # MeV/c²
ratio_CODATA = 206.7682827

deviation_percent = abs(m_mu - m_mu_CODATA) / m_mu_CODATA * 100
deviation_ratio_percent = abs(ratio_mu_e - ratio_CODATA) / ratio_CODATA * 100

print("CALIBRATION")
print("-" * 70)
print(f"TriPhase m_μ        = {m_mu_MeV:.10f} MeV/c²")
print(f"CODATA 2018 m_μ     = {m_mu_MeV_CODATA:.10f} MeV/c²")
print(f"Deviation           = {deviation_percent:.2f}%")
print()
print(f"TriPhase m_μ/m_e    = {ratio_mu_e:.10f}")
print(f"CODATA 2018 m_μ/m_e = {ratio_CODATA:.10f}")
print(f"Deviation           = {deviation_ratio_percent:.2f}%")
print()
print("NOTE: Large deviation suggests formula needs refinement or different interpretation.")
print("The factor 3×T₁₇/α ≈ 63,000 is far from observed ratio ≈ 207.")
print()

print("RENORMALIZATION GROUP INSIGHT")
print("-" * 70)
print("Lepton masses run from UV Yukawa couplings to IR electroweak breaking scale.")
print("The muon is at a different IR fixed point than the electron in the Higgs potential.")
print("TriPhase proposes topology (T₁₇) and generation factor (3) govern mass ratios.")
print()
print("=" * 70)

input("Press Enter to exit...")
