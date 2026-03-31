"""
TriPhase V16 - Down Quark Mass (QFT Framework)
===============================================

QFT INTERPRETATION:
The down quark (d) is the second-lightest quark, partner to the up quark:
- First generation, charge -e/3, color triplet (3 color states)
- Yukawa coupling to Higgs: m_d = y_d v/√2 with y_d ≈ 2×10⁻⁵
- Constituent of proton: p = uud, neutron: n = udd
- Slightly heavier than up quark: m_d/m_u ≈ 2.3
- Current mass m_d ~ 4.7 MeV, constituent mass M_d ~ 300 MeV

TriPhase's formula m_d = m_e × 9 × (1 + α/π) gives current down quark mass.
The factor 9 = 3² suggests a connection to color SU(3) or three generations.
The ratio m_d/m_u = 9/4 = 2.25 emerges naturally from geometric factors,
explaining the small isospin breaking in nuclear physics.

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

# ========== QFT DERIVATION: DOWN QUARK MASS ==========
print("=" * 70)
print("TriPhase V16 - Down Quark Mass")
print("QFT Framework: Isospin Breaking & QCD Chiral Dynamics")
print("=" * 70)
print()

print("QFT CONTEXT:")
print("The down quark is the SU(2) isospin partner of the up quark:")
print()
print("  Q_L = (u, d)_L  forms isospin doublet with I₃ = (+1/2, -1/2)")
print()
print("The mass difference m_d - m_u ≈ 2.5 MeV breaks isospin symmetry,")
print("explaining why neutron is heavier than proton:")
print("  m_n - m_p ≈ 1.3 MeV (dominated by m_d > m_u)")
print()
print("In QCD, the ratio m_d/m_u is crucial for:")
print("- Neutron-proton mass splitting")
print("- Pion mass differences: m_π± vs m_π⁰")
print("- Strong CP problem constraints")
print()

print("TRIPHASE DERIVATION:")
print("m_d = m_e × 9 × (1 + α/π)")
print()
print(f"Electron mass:        m_e = {m_e:.10e} kg")
print(f"Color factor:         9 = 3²")
print(f"QED correction:       1 + α/π = {1.0 + alpha/math.pi:.10f}")
print()

m_d = m_e * 9.0 * (1.0 + alpha/math.pi)

print(f"m_d (TriPhase):       {m_d:.10e} kg")
print()

# Convert to MeV
m_d_MeV = m_d * c**2 / 1.602176634e-13
print(f"m_d c² (MeV):         {m_d_MeV:.6f} MeV")
print()

# Mass ratios
m_u_MeV = m_e * 4.0 * (1.0 + alpha/math.pi) * c**2 / 1.602176634e-13
print(f"Ratio m_d/m_e:        {m_d / m_e:.6f}")
print(f"Ratio m_d/m_u:        {9.0 / 4.0:.6f} = 9/4 (exact)")
print(f"                      {m_d_MeV / m_u_MeV:.6f} (computed)")
print()

# ========== CALIBRATION CHECKPOINT ==========
pdg_m_d = 4.67  # MeV (PDG 2020, MS-bar scheme at 2 GeV)
deviation_ppm = (m_d_MeV - pdg_m_d) / pdg_m_d * 1e6

print("CALIBRATION CHECKPOINT:")
print(f"PDG 2020 (MS-bar):    {pdg_m_d:.2f} MeV (at μ = 2 GeV)")
print(f"TriPhase:             {m_d_MeV:.2f} MeV")
print(f"Deviation:            {deviation_ppm:+.0f} ppm")
print()
print("Note: Current quark masses are scheme and scale dependent.")
print("Typical range for m_d: 4.4 - 5.2 MeV")
print()

# PDG ratio
pdg_ratio = 4.67 / 2.16
print(f"PDG m_d/m_u:          {pdg_ratio:.3f}")
print(f"TriPhase m_d/m_u:     {9.0/4.0:.3f}")
print(f"Agreement:            {abs(pdg_ratio - 9.0/4.0) / pdg_ratio * 100:.1f}% deviation")
print()

# ========== QFT INSIGHT ==========
print("QFT INSIGHT:")
print("The ratio m_d/m_u = 9/4 = 2.25 is remarkably close to the PDG value")
print("m_d/m_u ≈ 2.16, suggesting a geometric origin for isospin breaking.")
print()
print("The factor 9 = 3² has multiple interpretations:")
print()
print("1. COLOR SU(3): Three color states, with 3² = 9 from color-anticolor")
print("   combinations in gluon exchange diagrams")
print()
print("2. GENERATIONAL STRUCTURE: Three quark generations, 3² from nested")
print("   flavor symmetry breaking")
print()
print("3. ISOSPIN + COLOR: I₃ = -1/2 state (down) with 3 colors × 3 combinatorial")
print("   factor from QCD vertex structure")
print()
print("The near-equality of QED correction factors for u and d quarks (both α/π)")
print("suggests electromagnetic contributions to quark masses are flavor-universal,")
print("while the geometric factors {4, 9} distinguish isospin states.")
print()
print("PHYSICAL CONSEQUENCE:")
print("The ratio m_d/m_u = 9/4 determines the neutron-proton mass difference,")
print("which in turn sets the stability of hydrogen (if m_n < m_p, hydrogen would")
print("decay!). This geometric ratio is thus crucial for chemistry and life.")
print()
print("The fact that m_d > m_u by exactly the right amount (neither too large nor")
print("too small) appears as a 'fine-tuning' in the Standard Model but emerges")
print("naturally in TriPhase as a simple ratio of geometric mode numbers: 9/4.")
print()
print("=" * 70)

input("Press Enter to exit...")
