"""
TriPhase V16 - Up Quark Mass (QFT Framework)
=============================================

QFT INTERPRETATION:
The up quark (u) is the lightest quark and a fundamental constituent of matter:
- First generation, charge +2e/3, color triplet (3 color states)
- Yukawa coupling to Higgs: m_u = y_u v/√2 with y_u ≈ 1×10⁻⁵
- Constituent of proton: p = uud, neutron: n = udd
- QCD confinement: quarks never observed in isolation
- Current quark mass vs constituent mass: m_u(current) ~ 2 MeV, m_u(constituent) ~ 300 MeV

TriPhase's formula m_u = m_e × 4 × (1 + α/π) gives the current quark mass
(Lagrangian parameter) rather than the constituent mass. The factor 4 may
relate to SU(2) isospin doublet structure, and α/π correction represents
QCD-QED mixing at one-loop level.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
TAG: (D*) - Derived with discrete geometric selection
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

# ========== QFT DERIVATION: UP QUARK MASS ==========
print("=" * 70)
print("TriPhase V16 - Up Quark Mass")
print("QFT Framework: Current Quark Mass & Chiral Symmetry Breaking")
print("=" * 70)
print()

print("QFT CONTEXT:")
print("In QCD, quarks have two mass concepts:")
print()
print("1. CURRENT MASS (m_u): Bare mass parameter in QCD Lagrangian")
print("   - Renormalization scale dependent (typically evaluated at μ = 2 GeV)")
print("   - Light quarks: m_u ≈ 2 MeV, m_d ≈ 5 MeV")
print("   - Breaks chiral symmetry explicitly")
print()
print("2. CONSTITUENT MASS (M_u): Effective mass inside hadrons")
print("   - Includes gluon field energy and quark condensates")
print("   - M_u ≈ 300 MeV from non-perturbative QCD dynamics")
print()
print("TriPhase predicts the current quark mass, the fundamental parameter.")
print()

print("TRIPHASE DERIVATION:")
print("m_u = m_e × 4 × (1 + α/π)")
print()
print(f"Electron mass:        m_e = {m_e:.10e} kg")
print(f"Isospin factor:       4")
print(f"QED correction:       1 + α/π = {1.0 + alpha/math.pi:.10f}")
print()

m_u = m_e * 4.0 * (1.0 + alpha/math.pi)

print(f"m_u (TriPhase):       {m_u:.10e} kg")
print()

# Convert to MeV
m_u_MeV = m_u * c**2 / 1.602176634e-13
print(f"m_u c² (MeV):         {m_u_MeV:.6f} MeV")
print()

# Mass ratios
print(f"Ratio m_u/m_e:        {m_u / m_e:.6f}")
print()

# ========== CALIBRATION CHECKPOINT ==========
pdg_m_u = 2.16  # MeV (PDG 2020, MS-bar scheme at 2 GeV)
deviation_ppm = (m_u_MeV - pdg_m_u) / pdg_m_u * 1e6

print("CALIBRATION CHECKPOINT:")
print(f"PDG 2020 (MS-bar):    {pdg_m_u:.2f} MeV (at μ = 2 GeV)")
print(f"TriPhase:             {m_u_MeV:.2f} MeV")
print(f"Deviation:            {deviation_ppm:+.0f} ppm")
print()
print("Note: Current quark masses are scheme and scale dependent.")
print("Typical range for m_u: 1.8 - 2.4 MeV")
print()

# ========== QFT INSIGHT ==========
print("QFT INSIGHT:")
print("The factor 4 in m_u = 4 m_e × (1 + α/π) suggests a connection to SU(2)")
print("isospin symmetry. In the Standard Model:")
print()
print("  Q_L = (u, d)_L forms an SU(2)_L doublet")
print()
print("The number 4 = 2² may represent:")
print("- Two spin states × two isospin states")
print("- Four-component Dirac spinor structure")
print("- SU(2) Casimir invariant")
print()
print("The QED correction α/π (versus α/2π for leptons) suggests enhanced")
print("electromagnetic-QCD mixing at one loop. This could arise from:")
print()
print("  γ → qq̄ (photon splitting into quark pairs)")
print("  Virtual quark loops modifying photon propagator")
print()
print("The fact that m_u ~ 4 m_e (up to radiative corrections) hints at a deep")
print("unification: quarks and leptons may be different electromagnetic modes of")
print("the same underlying field, distinguished by color and isospin quantum numbers.")
print()
print("The small current quark mass (m_u ≪ M_proton) explains why chiral symmetry")
print("is approximate in QCD. Most of the proton mass comes from gluon field energy,")
print("not Higgs-generated quark masses—a pure QCD phenomenon.")
print()
print("=" * 70)

input("Press Enter to exit...")
