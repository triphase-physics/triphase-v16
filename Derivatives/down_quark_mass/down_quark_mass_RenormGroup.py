"""
TriPhase V16 — Down Quark Mass (Renormalization Group Framework)
=================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

RENORMALIZATION GROUP INTERPRETATION
--------------------------------------
The down quark mass m_d ≈ 4.7 MeV/c² (MS-bar scheme at 2 GeV) is roughly twice the
up quark mass, reflecting the isospin structure of the light quark sector. Like the
up quark, the down quark mass is a running mass that varies with energy scale due to
QCD's asymptotic freedom and infrared confinement. The ratio m_d/m_u ≈ 2.2 is nearly
RG-invariant, suggesting both quarks flow to related IR fixed points in the QCD/Higgs
coupling space.

The TriPhase formula m_d = m_e × (4/3) α T₁₇ derives the down quark mass from the
electron mass with the charge magnitude factor 4/3 (down quark has charge -1/3 e,
and the factor 4/3 = 2 × 2/3 encodes the isospin doubling). This gives m_d/m_u = 2,
precisely the observed ratio to leading order. The α suppression and T₁₇ topology
factor are the same as for the up quark, showing that u and d quarks sit at the
same RG scale but different isospin subspaces.

The fact that m_d > m_u (despite having smaller charge magnitude) is a puzzle in the
Standard Model, resolved here by the isospin structure: the factor 4/3 vs 2/3 arises
from SU(2)_L weak isospin, where up and down quarks form a doublet. The RG flow
preserves this isospin structure, with m_d/m_u ≈ 2 emerging as a fixed point ratio.

TAG: (D*) — Derived with discrete selection (quark charge fraction -1/3, isospin factor)
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

# Up quark mass (for comparison)
m_u = m_e * (2.0 / 3.0) * alpha * T_17

# ========== RENORMALIZATION GROUP DERIVATION ==========
print("=" * 70)
print("TriPhase V16: Down Quark Mass (Renormalization Group)")
print("=" * 70)
print()

print("IR RUNNING MASS AT HADRONIC SCALE (ISOSPIN DOUBLET)")
print("-" * 70)
print("Electron mass (lepton IR fixed point):")
m_e_MeV = m_e * c**2 / e / 1e6
print(f"  m_e = {m_e:.15e} kg")
print(f"      = {m_e_MeV:.10f} MeV/c²")
print()

print("Down quark mass (charge -1/3, isospin factor 4/3, topology):")
print(f"  m_d = m_e × (4/3) × α × T₁₇")
print(f"      = {m_e:.10e} × (4/3) × {alpha:.10f} × {T_17}")
print(f"      = {m_e:.10e} × {(4.0/3.0) * alpha * T_17:.10f}")
print()

m_d = m_e * (4.0 / 3.0) * alpha * T_17
m_d_MeV = m_d * c**2 / e / 1e6

print(f"  m_d = {m_d:.15e} kg")
print(f"      = {m_d_MeV:.10f} MeV/c²")
print()

# Mass ratio to up quark
m_u_MeV = m_u * c**2 / e / 1e6
ratio_d_u = m_d / m_u
print(f"Up quark mass (for comparison):")
print(f"  m_u = {m_u_MeV:.10f} MeV/c²")
print()
print(f"Isospin mass ratio:")
print(f"  m_d / m_u = {ratio_d_u:.10f} = (4/3) / (2/3) = 2.0")
print()

# ========== CALIBRATION CHECKPOINT ==========
m_d_PDG_MS_bar = 4.67  # MeV/c² (MS-bar at 2 GeV, PDG 2022)
m_d_PDG_range = (4.1, 5.3)  # MeV/c² (uncertainty range)
ratio_d_u_PDG = 2.16  # PDG 2022

print("CALIBRATION (PDG 2022)")
print("-" * 70)
print(f"TriPhase m_d        = {m_d_MeV:.10f} MeV/c² (current mass estimate)")
print(f"PDG MS-bar m_d(2GeV)= {m_d_PDG_MS_bar:.2f} MeV/c² (range: {m_d_PDG_range[0]}-{m_d_PDG_range[1]} MeV)")
print()

deviation_percent = abs(m_d_MeV - m_d_PDG_MS_bar) / m_d_PDG_MS_bar * 100
print(f"Deviation           = {deviation_percent:.1f}%")
print()

print(f"TriPhase m_d/m_u    = {ratio_d_u:.2f}")
print(f"PDG m_d/m_u         = {ratio_d_u_PDG:.2f}")
print(f"Ratio deviation     = {abs(ratio_d_u - ratio_d_u_PDG) / ratio_d_u_PDG * 100:.1f}%")
print()

print("NOTE: Quark masses are highly RG-dependent due to QCD running.")
print("  - Current quark mass (MS-bar, 2 GeV): m_d ~ 4.7 MeV (perturbative)")
print("  - Constituent quark mass (IR, 1 GeV): m_d^const ~ 300 MeV (non-perturbative)")
print()
print(f"TriPhase predicts current mass m_d ~ {m_d_MeV:.2f} MeV, within factor ~5-6 of PDG.")
print(f"The isospin ratio m_d/m_u = 2.0 matches observed ratio 2.16 to ~8%.")
print()

print("RENORMALIZATION GROUP INSIGHT")
print("-" * 70)
print("Down quark mass runs with μ like up quark, but with isospin factor 4/3 vs 2/3.")
print("The ratio m_d/m_u ≈ 2 is nearly RG-invariant, an IR fixed point of SU(2)_L isospin.")
print("Both u and d quarks sit at the same RG scale (αT₁₇ suppression), different isospin channels.")
print()
print("=" * 70)

input("Press Enter to exit...")
