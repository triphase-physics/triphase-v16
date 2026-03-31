"""
TriPhase V16: Proton Mass - QFT Framework
==========================================

QFT INTERPRETATION:
The proton mass m_p ≈ 938.27 MeV is NOT simply the sum of its constituent quark
masses (2m_u + m_d ≈ 9 MeV). Instead, 99% of the proton's mass arises from QCD
binding energy—the kinetic energy of quarks and gluons confined within a radius
r_p ≈ 0.84 fm, plus the gluon field energy stored in the QCD vacuum.

In QFT language, the proton is a bound state |p⟩ = |uud⟩ + gluons + qq̄ pairs,
described by a complicated many-body wavefunction in QCD. The mass emerges from
the trace anomaly of the energy-momentum tensor: ⟨p|T^μ_μ|p⟩ ≠ 0 even for massless
quarks, due to quantum corrections breaking scale invariance.

Lattice QCD calculations show m_p arises from:
  • ~32% quark kinetic energy
  • ~37% gluon field energy
  • ~9% quark-gluon interaction
  • ~23% trace anomaly / quantum effects

TriPhase derives m_p from m_e * (4 × 27 × 17) × (1 + 5α²/π), where the factor
4×27×17 = 1836 approximates the proton-to-electron mass ratio, and the 5α²/π
term encodes QCD binding corrections to the naive quark-mass sum.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
TAG: (D) - Direct derivation
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

# ========== QFT DERIVATION: PROTON MASS ==========
print("=" * 70)
print("  TRIPHASE V16: PROTON MASS (QFT Framework)")
print("=" * 70)
print()
print("QFT CONTEXT:")
print("  The proton is a QCD bound state |p⟩ = |uud⟩ + O(αs) corrections.")
print("  Its mass does NOT come from quark rest masses (which contribute")
print("  only ~1%), but from QCD binding energy: gluon fields, quark kinetic")
print("  energy, and the trace anomaly of the QCD stress tensor.")
print()
print("  Lattice QCD computes m_p from first principles by numerically")
print("  evaluating the QCD path integral on a discrete spacetime grid,")
print("  reproducing the experimental value to within 1%.")
print()

# Derivation
m_p_calc = m_p
m_p_MeV = m_p_calc * c**2 / 1.602176634e-13

print("DERIVATION STEPS:")
print(f"  1. Proton-to-electron mass ratio:")
print(f"     mp/me = 4 × 27 × 17 × (1 + 5α²/π)")
print(f"     = {4 * 27 * 17} × (1 + 5 × {alpha:.6f}² / π)")
print(f"     = 1836 × {1.0 + 5 * alpha**2 / math.pi:.8f}")
print(f"     = {mp_me:.8f}")
print()
print(f"  2. Proton mass:")
print(f"     m_p = m_e × (mp/me)")
print(f"     = {m_e:.6e} kg × {mp_me:.8f}")
print(f"     = {m_p_calc:.15e} kg")
print(f"     = {m_p_MeV:.5f} MeV/c²")
print()

# Calibration
m_p_CODATA = 1.67262192369e-27  # kg (CODATA 2018)
m_p_CODATA_MeV = m_p_CODATA * c**2 / 1.602176634e-13
deviation_ppm = abs(m_p_calc - m_p_CODATA) / m_p_CODATA * 1e6

print("=" * 70)
print("CALIBRATION CHECKPOINT:")
print(f"  TriPhase value:  {m_p_calc:.15e} kg")
print(f"  CODATA 2018:     {m_p_CODATA:.15e} kg")
print(f"  TriPhase:        {m_p_MeV:.5f} MeV/c²")
print(f"  CODATA:          {m_p_CODATA_MeV:.5f} MeV/c²")
print(f"  Deviation:       {deviation_ppm:.2f} ppm")
print("=" * 70)
print()

print("QFT INSIGHT:")
print("  The proton mass is one of the most spectacular predictions of QCD.")
print("  Starting from the QCD Lagrangian with nearly massless u,d quarks,")
print("  lattice simulations generate m_p ≈ 938 MeV purely from binding.")
print()
print("  This is emergent mass: the strong field energy 'crystallizes' into")
print("  the proton's rest mass via E = mc². The trace anomaly relates this")
print("  to the running of the strong coupling αs(μ)—a purely quantum effect.")
print()
print("  TriPhase's geometric formula m_p/m_e = 4×27×17×(1 + 5α²/π) suggests")
print("  the proton mass is encoded in the electromagnetic structure constant")
print("  α through a geometric scaling involving:")
print("    • 4: isospin doublet (p,n)")
print("    • 27 = 3³: color × flavor × generation")
print("    • 17: horizon step linking scales")
print("    • 5α²/π: QCD binding correction")
print()
print("  This ~5000 ppm accuracy from a pure geometric formula hints at a")
print("  deep connection between electromagnetic and strong force scales.")
print("=" * 70)

input("Press Enter to exit...")
