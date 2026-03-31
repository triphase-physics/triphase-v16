"""
========================================================================
TriPhase V16 Derivative: Up Quark Mass (Gauge Theory)
========================================================================

GAUGE THEORY INTERPRETATION:
The up quark mass m_u ≈ 2.2 MeV (current mass) arises from its Yukawa
coupling to the Higgs field in the electroweak sector. As a member of
the first generation, the up quark is the lightest quark and carries
electric charge +2e/3, coupling to the photon via U(1)_EM gauge symmetry.

The up quark participates in all three gauge interactions: U(1)_EM
(photon), SU(2)_L (W, Z bosons), and SU(3)_C (8 gluons). The current
mass m_u ≈ 2 MeV is distinct from the constituent mass ~300 MeV, which
includes QCD binding energy from the SU(3) color gauge field.

The TriPhase formula m_u = 4·m_e·(1 + α/π) gives a mass in the MeV range,
consistent with lattice QCD determinations. The factor 4 may relate to
the four color-spin states of a quark (3 colors × 2 spins / 1.5) or to
the difference in electromagnetic charge between quarks and leptons.
The radiative correction (1 + α/π) accounts for QED vertex corrections.

REFERENCE: PDG 2020 m_u = 2.16⁺⁰·⁴⁹₋₀.₂₆ MeV (MS scheme at 2 GeV)

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
TAG: (D*)
========================================================================
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

print("=" * 70)
print("GAUGE THEORY DERIVATION: Up Quark Current Mass")
print("=" * 70)

# Derive m_u from electron mass and gauge correction
m_u = m_e * 4.0 * (1.0 + alpha/math.pi)

print("\nFirst-Generation Quark Mass (Higgs Yukawa + QCD):")
print(f"  m_e (electron mass):         {m_e:.15e} kg")
print(f"  Factor:                      4")
print(f"  α/π (radiative correction):  {alpha/math.pi:.15e}")
print(f"  Correction factor:           1 + α/π = {1.0 + alpha/math.pi:.15f}")
print(f"  m_u = 4·m_e·(1+α/π):         {m_u:.15e} kg")

# Convert to MeV
m_u_MeV = m_u * c**2 / 1.602176634e-13

print(f"\nUp quark mass in MeV:")
print(f"  m_u c²:                      {m_u_MeV:.6f} MeV")

# Constituent vs current mass
m_u_constituent = 300  # MeV (approximate)
print(f"\nCurrent vs. constituent mass:")
print(f"  Current mass (Higgs):        {m_u_MeV:.2f} MeV")
print(f"  Constituent mass (QCD):      ~{m_u_constituent} MeV")
print(f"  Ratio:                       {m_u_constituent/m_u_MeV:.1f}")

print("\n" + "=" * 70)
print("CALIBRATION CHECKPOINT")
print("=" * 70)

m_u_PDG = 2.16  # MeV (central value)
m_u_PDG_upper = 2.65
m_u_PDG_lower = 1.90

deviation = abs(m_u_MeV - m_u_PDG)
deviation_percent = (deviation / m_u_PDG) * 100

print(f"\nTriPhase m_u:     {m_u_MeV:.6f} MeV")
print(f"PDG 2020:         {m_u_PDG:.2f} +{m_u_PDG_upper-m_u_PDG:.2f} -{m_u_PDG-m_u_PDG_lower:.2f} MeV")
print(f"Deviation:        {deviation:.6f} MeV")
print(f"Deviation (%):    {deviation_percent:.3f} %")

if m_u_MeV >= m_u_PDG_lower and m_u_MeV <= m_u_PDG_upper:
    print("✓ Within PDG uncertainty range")
elif deviation_percent < 50:
    print("✓ Correct order of magnitude")
else:
    print("⚠ Outside PDG range")

print("\n" + "=" * 70)
print("GAUGE THEORY INSIGHTS")
print("=" * 70)

print("""
The up quark mass in gauge theory:

1. QUARK QUANTUM NUMBERS:
   - Charge: Q = +2e/3 (fractional!)
   - Isospin: I₃ = +1/2 (SU(2) weak)
   - Color: r, g, b (SU(3) strong triplet)
   - Baryon number: B = 1/3
   - First generation: (u, d) doublet

2. GAUGE INTERACTIONS:
   - U(1)_EM: Couples to photon with Q = +2/3
   - SU(2)_L: Left-handed doublet (u, d)_L, right-handed singlet u_R
   - SU(3)_C: Color triplet, couples to 8 gluons
   - All three SM gauge groups apply to quarks!

3. YUKAWA COUPLING TO HIGGS:
   - Higgs doublet: φ = (φ⁺, φ⁰)ᵀ
   - Up-type coupling: L_Y = -y_u Q̄_L φ̃ u_R + h.c.
   - φ̃ = iσ₂φ* (conjugate doublet for up quarks)
   - Mass: m_u = y_u v/√2 ≈ 2.2 MeV

4. CURRENT vs. CONSTITUENT MASS:
   - Current mass: m_u ~ 2 MeV (Higgs/Yukawa origin)
   - Constituent mass: M_u ~ 300 MeV (includes QCD binding)
   - Proton mass: m_p ≈ 938 MeV ≈ 2M_u + M_d + binding
   - ~98% of proton mass from QCD, not Higgs!

5. CHIRAL SYMMETRY BREAKING:
   - If m_u = 0: QCD has chiral SU(3)_L × SU(3)_R symmetry
   - QCD vacuum breaks this to SU(3)_V (vector)
   - Generates constituent quark mass dynamically
   - Goldstone bosons: π, K, η (pseudo-Nambu-Goldstone)

6. LATTICE QCD:
   - Non-perturbative calculation of m_u
   - Discretize spacetime, solve QCD numerically
   - MS scheme at μ = 2 GeV: m_u = 2.16 MeV
   - Large systematic uncertainties from QCD

7. FACTOR 4 INTERPRETATION:
   - 4 = 2 × 2 (color factor × isospin?)
   - 4/3 × 3 = 4 (charge ratio × color)
   - (Q_u / Q_e)² = (2/3)² / 1² ≈ 0.44 (not quite 4)
   - May be geometric or group-theoretic

8. RADIATIVE CORRECTION (1 + α/π):
   - QED vertex correction (same as tau)
   - Running Yukawa coupling: y_u(μ) varies with scale
   - Includes EM and weak loop corrections
   - Strong coupling αₛ dominates at low energy

9. QUARK CONFINEMENT:
   - Free quarks never observed (color confinement)
   - Energy to separate quarks → creates new qq̄ pair
   - Asymptotic freedom: αₛ(high E) → 0
   - Confinement scale: Λ_QCD ≈ 200 MeV

10. CKM MATRIX:
    - Quark mixing: Flavor eigenstates ≠ mass eigenstates
    - (d', s', b') = V_CKM (d, s, b)
    - Up quark couples to down, strange, bottom via W
    - V_ud ≈ 0.974 (dominant), V_us ≈ 0.224, V_ub ≈ 0.004

11. PION DECAY:
    - π⁺ = ūd (up antiquark + down quark)
    - π⁰ = (uū - dd̄)/√2 (I=1 superposition)
    - π mass: m_π ≈ 140 MeV ∝ √(m_u + m_d)
    - Goldstone boson (pseudo, since m_q ≠ 0)

12. QUARK MASS HIERARCHY:
    - Up: m_u ~ 2 MeV
    - Down: m_d ~ 5 MeV
    - Strange: m_s ~ 95 MeV
    - Charm: m_c ~ 1.3 GeV
    - Bottom: m_b ~ 4.2 GeV
    - Top: m_t ~ 173 GeV
    - 5 orders of magnitude! (why?)

13. GRAND UNIFICATION:
    - SO(10): Quarks and leptons in 16-plet
    - (u, d, νₑ, e)_L,R + ν_R (sterile)
    - Yukawa unification at GUT scale?
    - y_u(M_GUT) ≈ y_e(M_GUT)? (not observed)

The up quark is the lightest quark and a cornerstone of visible matter
(protons are uud, neutrons are udd). Its current mass m_u ~ 2 MeV from
the Higgs is tiny compared to its constituent mass ~300 MeV from QCD.
The TriPhase formula m_u = 4·m_e·(1 + α/π) suggests a connection between
quark and lepton masses through the factor 4, possibly reflecting charge
or color structure in the gauge theory.
""")

print("=" * 70)
input("Press Enter to exit...")
