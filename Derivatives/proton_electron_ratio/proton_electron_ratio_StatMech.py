"""
TriPhase V16 — Proton-Electron Mass Ratio (Statistical Mechanics Framework)
========================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
The proton-electron mass ratio emerges from the canonical ensemble of composite
hadronic states. The proton is a bound state of three quarks in QCD, and its mass
arises primarily from the gluon field energy trapped in confinement. This binding
energy can be understood through the partition function Z = Tr[exp(-βH_QCD)], where
the Hamiltonian includes quark kinetic terms, gluon field energy, and the α_s coupling.

The TriPhase formula m_p/m_e = 4·27·17·(1 + 5α²/π) reveals the statistical structure.
The factor 4·27·17 = 1836 counts the degrees of freedom in the hadronic phase space:
4 (spin-isospin), 27 (color states = 3³), and 17 (triangular number T_17 from spatial
wave modes). The correction term 5α²/π arises from virtual photon loops that couple
the electromagnetic and strong sectors—a grand canonical ensemble effect where particle
number fluctuations (virtual e⁺e⁻ pairs) contribute to the effective mass.

TAG: (D) — Direct TriPhase derivation from pure wave mechanics
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
print("TriPhase V16: Proton-Electron Mass Ratio (Statistical Mechanics)")
print("=" * 70)
print()

print("CANONICAL ENSEMBLE OF HADRONIC STATES:")
print("-" * 70)
print("The proton mass emerges from the QCD partition function.")
print("Phase space factors count accessible microstates:")
print()

factor_spin = 4.0
factor_color = 27.0  # 3^3 for three quarks in color SU(3)
factor_spatial = 17.0  # T_17 triangular number, spatial wave modes

print(f"  Spin-isospin degeneracy:      g_spin = {factor_spin:.0f}")
print(f"  Color state multiplicity:     g_color = {factor_color:.0f} (3³)")
print(f"  Spatial mode count:           g_spatial = {factor_spatial:.0f} (T₁₇)")
print()

base_ratio = factor_spin * factor_color * factor_spatial
print(f"Base ratio (microcanonical): m_p/m_e ≈ {base_ratio:.0f}")
print()

print("ELECTROMAGNETIC CORRECTION (GRAND CANONICAL):")
print("-" * 70)
print("Virtual photon loops couple EM and QCD sectors.")
print("Particle number fluctuations contribute to effective mass:")
print()

em_correction = 5.0 * alpha**2 / math.pi
print(f"  Correction factor:  δ = 5α²/π = {em_correction:.8f}")
print(f"  Renormalized ratio: (m_p/m_e) = {base_ratio}·(1 + {em_correction:.6f})")
print()

mp_me_calc = mp_me
print(f"  Final result:  m_p/m_e = {mp_me_calc:.10f}")
print()

print("PARTITION FUNCTION STRUCTURE:")
print("-" * 70)
print("  Z_hadron = g_total · exp(-β·m_p c²)")
print(f"  where g_total = 4·27·17·(1 + 5α²/π) ≈ {mp_me_calc:.0f}")
print()
print("This explains why the proton is ~1836 times heavier than the electron:")
print("it's the ratio of hadronic to leptonic phase space volumes.")
print()

# ========== CALIBRATION CHECKPOINT ==========
mp_me_CODATA = 1836.15267343  # CODATA 2018
deviation_ppm = (mp_me_calc - mp_me_CODATA) / mp_me_CODATA * 1e6

print("=" * 70)
print("CALIBRATION CHECKPOINT")
print("=" * 70)
print(f"CODATA 2018:            m_p/m_e = {mp_me_CODATA:.10f}")
print(f"TriPhase V16 (StatMech):        = {mp_me_calc:.10f}")
print(f"Deviation:                        {deviation_ppm:+.2f} ppm")
print("=" * 70)
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT:")
print("-" * 70)
print("The proton mass is not fundamental—it emerges from the statistical")
print("ensemble of QCD vacuum states. The factor 1836 is the ratio of")
print("accessible microstates in hadronic vs. leptonic phase space.")
print()
print("The decomposition 1836 = 4·27·17 reveals the underlying structure:")
print("  • 4 from spin-isospin (2×2 quantum numbers)")
print("  • 27 from color confinement (3³ for three quarks)")
print("  • 17 from spatial wave modes (triangular number T₁₇)")
print()
print("The correction 5α²/π is a grand canonical effect: virtual e⁺e⁻ pairs")
print("fluctuate in the EM sector and modify the effective hadronic mass.")
print("This is the statistical origin of the proton's identity.")
print("=" * 70)

input("Press Enter to exit...")
