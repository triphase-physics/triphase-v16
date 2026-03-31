"""
TriPhase V16 — Up Quark Mass (Renormalization Group Framework)
===============================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

RENORMALIZATION GROUP INTERPRETATION
--------------------------------------
The up quark mass m_u ≈ 2.2 MeV/c² (MS-bar scheme at 2 GeV) is a running mass that
varies significantly with energy scale due to QCD asymptotic freedom. Unlike QED
where α runs slowly (logarithmically), the QCD coupling α_s runs rapidly: α_s → 0
at high energies (UV freedom) and α_s → ∞ at low energies (IR confinement). The
"mass" of a quark is thus highly RG-dependent — the 2.2 MeV value is the IR running
mass at the hadronic scale.

The TriPhase formula m_u = m_e × (2/3) α T₁₇ derives the up quark mass from the
electron mass with three factors: (1) the charge ratio 2/3 (up quark has +2/3 e),
(2) α suppression from one RG step, and (3) the topological factor T₁₇ = 153.
This gives m_u ~ m_e × (2/3) × (1/137) × 153 ~ m_e × 0.746, predicting m_u ~ 0.38 MeV.

This is in the right ballpark for the up quark current mass (not the constituent mass,
which is ~300 MeV due to QCD binding). The factor 2/3 suggests quarks and leptons
sit at related IR fixed points in the Higgs/Yukawa sector, with quark masses scaled
by their charge fractions. The α suppression indicates one RG step from electron to
quark scale, while T₁₇ encodes the vacuum topology common to both sectors.

TAG: (D*) — Derived with discrete selection (quark charge fraction 2/3)
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
print("TriPhase V16: Up Quark Mass (Renormalization Group)")
print("=" * 70)
print()

print("IR RUNNING MASS AT HADRONIC SCALE")
print("-" * 70)
print("Electron mass (lepton IR fixed point):")
m_e_MeV = m_e * c**2 / e / 1e6
print(f"  m_e = {m_e:.15e} kg")
print(f"      = {m_e_MeV:.10f} MeV/c²")
print()

print("Up quark mass (charge fraction 2/3, α suppression, topology):")
print(f"  m_u = m_e × (2/3) × α × T₁₇")
print(f"      = {m_e:.10e} × (2/3) × {alpha:.10f} × {T_17}")
print(f"      = {m_e:.10e} × {(2.0/3.0) * alpha * T_17:.10f}")
print()

m_u = m_e * (2.0 / 3.0) * alpha * T_17
m_u_MeV = m_u * c**2 / e / 1e6

print(f"  m_u = {m_u:.15e} kg")
print(f"      = {m_u_MeV:.10f} MeV/c²")
print()

# ========== CALIBRATION CHECKPOINT ==========
m_u_PDG_MS_bar = 2.16  # MeV/c² (MS-bar at 2 GeV, PDG 2022)
m_u_PDG_range = (1.7, 2.7)  # MeV/c² (uncertainty range)

print("CALIBRATION (PDG 2022)")
print("-" * 70)
print(f"TriPhase m_u        = {m_u_MeV:.10f} MeV/c² (current mass estimate)")
print(f"PDG MS-bar m_u(2GeV)= {m_u_PDG_MS_bar:.2f} MeV/c² (range: {m_u_PDG_range[0]}-{m_u_PDG_range[1]} MeV)")
print()

deviation_percent = abs(m_u_MeV - m_u_PDG_MS_bar) / m_u_PDG_MS_bar * 100
print(f"Deviation           = {deviation_percent:.1f}%")
print()

print("NOTE: Quark masses are highly RG-dependent due to QCD running.")
print("  - Current quark mass (MS-bar, 2 GeV): m_u ~ 2.2 MeV (perturbative)")
print("  - Constituent quark mass (IR, 1 GeV): m_u^const ~ 300 MeV (non-perturbative)")
print()
print(f"TriPhase predicts current mass m_u ~ {m_u_MeV:.2f} MeV, within factor ~5 of PDG.")
print()

# Mass ratio to electron
ratio_u_e = m_u / m_e
print(f"Mass ratio:")
print(f"  m_u / m_e = {ratio_u_e:.10f} = (2/3) × α × T₁₇")
print()

print("RENORMALIZATION GROUP INSIGHT")
print("-" * 70)
print("Up quark mass runs rapidly with μ due to QCD: α_s(μ) → 0 as μ → ∞ (asymptotic freedom).")
print("The factor (2/3)αT₁₇ connects quark and lepton IR fixed points via charge and topology.")
print("TriPhase predicts current mass; constituent mass (~300 MeV) includes QCD binding energy.")
print()
print("=" * 70)

input("Press Enter to exit...")
