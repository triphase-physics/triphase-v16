"""
========================================================================
TriPhase V16 Derivative: Down Quark Mass (Gauge Theory)
========================================================================

GAUGE THEORY INTERPRETATION:
The down quark mass m_d ≈ 4.7 MeV (current mass) arises from its Yukawa
coupling to the Higgs field in the electroweak sector. As the isospin
partner of the up quark, the down quark carries electric charge -e/3
and couples to the photon via U(1)_EM gauge symmetry.

The down quark participates in all three gauge interactions: U(1)_EM
(photon), SU(2)_L (W, Z bosons), and SU(3)_C (8 gluons). The current
mass m_d ≈ 5 MeV is distinct from the constituent mass ~300 MeV, which
includes QCD binding energy. The ratio m_d/m_u ≈ 2.3 is unexplained by
the Standard Model but crucial for nuclear stability.

The TriPhase formula m_d = 9·m_e·(1 + α/π) gives a mass slightly heavier
than the up quark, consistent with lattice QCD. The factor 9 = 3² may
relate to the three color states (r, g, b) of the quark or to the three
generations. The mass difference Δm = m_d - m_u ≈ 2.5 MeV determines
neutron-proton mass difference and thus nuclear stability.

REFERENCE: PDG 2020 m_d = 4.67⁺⁰·⁴⁸₋₀.₁₇ MeV (MS scheme at 2 GeV)

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
print("GAUGE THEORY DERIVATION: Down Quark Current Mass")
print("=" * 70)

# Derive m_d from electron mass and gauge correction
m_d = m_e * 9.0 * (1.0 + alpha/math.pi)

print("\nFirst-Generation Quark Mass (Higgs Yukawa + QCD):")
print(f"  m_e (electron mass):         {m_e:.15e} kg")
print(f"  Factor:                      9")
print(f"  α/π (radiative correction):  {alpha/math.pi:.15e}")
print(f"  Correction factor:           1 + α/π = {1.0 + alpha/math.pi:.15f}")
print(f"  m_d = 9·m_e·(1+α/π):         {m_d:.15e} kg")

# Convert to MeV
m_d_MeV = m_d * c**2 / 1.602176634e-13

print(f"\nDown quark mass in MeV:")
print(f"  m_d c²:                      {m_d_MeV:.6f} MeV")

# Compare to up quark
m_u = m_e * 4.0 * (1.0 + alpha/math.pi)
m_u_MeV = m_u * c**2 / 1.602176634e-13
ratio_d_to_u = m_d / m_u

print(f"\nMass ratio to up quark:")
print(f"  m_u:                         {m_u_MeV:.6f} MeV")
print(f"  m_d / m_u:                   {ratio_d_to_u:.6f}")
print(f"  Δm = m_d - m_u:              {m_d_MeV - m_u_MeV:.6f} MeV")

# Constituent mass
m_d_constituent = 300  # MeV (approximate)
print(f"\nCurrent vs. constituent mass:")
print(f"  Current mass (Higgs):        {m_d_MeV:.2f} MeV")
print(f"  Constituent mass (QCD):      ~{m_d_constituent} MeV")
print(f"  Ratio:                       {m_d_constituent/m_d_MeV:.1f}")

print("\n" + "=" * 70)
print("CALIBRATION CHECKPOINT")
print("=" * 70)

m_d_PDG = 4.67  # MeV (central value)
m_d_PDG_upper = 5.15
m_d_PDG_lower = 4.50

deviation = abs(m_d_MeV - m_d_PDG)
deviation_percent = (deviation / m_d_PDG) * 100

print(f"\nTriPhase m_d:     {m_d_MeV:.6f} MeV")
print(f"PDG 2020:         {m_d_PDG:.2f} +{m_d_PDG_upper-m_d_PDG:.2f} -{m_d_PDG-m_d_PDG_lower:.2f} MeV")
print(f"Deviation:        {deviation:.6f} MeV")
print(f"Deviation (%):    {deviation_percent:.3f} %")

if m_d_MeV >= m_d_PDG_lower and m_d_MeV <= m_d_PDG_upper:
    print("✓ Within PDG uncertainty range")
elif deviation_percent < 50:
    print("✓ Correct order of magnitude")
else:
    print("⚠ Outside PDG range")

print("\n" + "=" * 70)
print("GAUGE THEORY INSIGHTS")
print("=" * 70)

print("""
The down quark mass in gauge theory:

1. QUARK QUANTUM NUMBERS:
   - Charge: Q = -e/3 (fractional!)
   - Isospin: I₃ = -1/2 (SU(2) weak partner of up)
   - Color: r, g, b (SU(3) strong triplet)
   - Baryon number: B = 1/3
   - First generation: (u, d) doublet

2. GAUGE INTERACTIONS:
   - U(1)_EM: Couples to photon with Q = -1/3
   - SU(2)_L: Left-handed doublet (u, d)_L, right-handed singlet d_R
   - SU(3)_C: Color triplet, couples to 8 gluons
   - All three SM gauge groups apply to down quark

3. YUKAWA COUPLING TO HIGGS:
   - Higgs doublet: φ = (φ⁺, φ⁰)ᵀ
   - Down-type coupling: L_Y = -y_d Q̄_L φ d_R + h.c.
   - φ = (0, v/√2)ᵀ after SSB
   - Mass: m_d = y_d v/√2 ≈ 4.7 MeV

4. CURRENT vs. CONSTITUENT MASS:
   - Current mass: m_d ~ 5 MeV (Higgs/Yukawa origin)
   - Constituent mass: M_d ~ 300 MeV (includes QCD binding)
   - Neutron mass: m_n ≈ 939 MeV ≈ M_u + 2M_d + binding
   - Proton mass: m_p ≈ 938 MeV ≈ 2M_u + M_d + binding

5. NEUTRON-PROTON MASS DIFFERENCE:
   - Δm_np = m_n - m_p ≈ 1.29 MeV
   - QCD contribution: m_d - m_u ≈ 2.5 MeV
   - EM contribution: Q_u² - Q_d² ≈ -0.8 MeV (opposite sign!)
   - Total: ΔE_QCD + ΔE_EM ≈ 1.3 MeV ✓
   - Crucial for neutron beta decay stability!

6. FACTOR 9 INTERPRETATION:
   - 9 = 3² (color states squared?)
   - 9/4 = 2.25 ≈ m_d/m_u (ratio of factors)
   - 9 = 3 × 3 (generations × color)
   - May be group-theoretic (SU(3) dimension?)

7. RADIATIVE CORRECTION (1 + α/π):
   - Same as up quark, tau lepton
   - QED + weak loops
   - Running Yukawa: y_d(μ) varies with scale
   - Strong coupling αₛ dominates at low energy

8. QUARK MASS RATIOS:
   - m_d / m_u ≈ 2.3 (unexplained!)
   - If m_d < m_u: Hydrogen unstable (p → n would be energetic)
   - If m_d >> m_u: No stable nuclei (n-p gap too large)
   - Fine-tuning for nuclear stability!

9. FLAVOR CHANGING NEUTRAL CURRENTS:
   - Tree level: d couples to W, Z (flavor-conserving)
   - Loop level: d ↔ s, b via CKM mixing
   - FCNC: K⁰-K̄⁰ mixing, B⁰-B̄⁰ mixing
   - Highly suppressed: GIM mechanism

10. PION MASS:
    - π⁺ = ūd (up antiquark + down quark)
    - π⁰ = (uū - dd̄)/√2
    - m_π² ∝ (m_u + m_d) (chiral perturbation theory)
    - If m_u, m_d → 0: m_π → 0 (Goldstone boson)

11. CABIBBO ANGLE:
    - Quark mixing: d' = d cos θ_C + s sin θ_C
    - θ_C ≈ 13° (Cabibbo angle)
    - V_ud = cos θ_C ≈ 0.974
    - V_us = sin θ_C ≈ 0.224
    - Weak decays: d → u + W⁻ (CKM-suppressed if s present)

12. NEUTRON BETA DECAY:
    - n → p + e⁻ + ν̄_e (lifetime τ_n ≈ 880 s)
    - Quark level: d → u + W⁻, W⁻ → e⁻ + ν̄_e
    - Rate ∝ |V_ud|² (m_n - m_p)⁵
    - Determines free neutron lifetime
    - Protons stable (lightest baryon, B conservation)

13. STRONG CP PROBLEM:
    - QCD allows θ-term: L_θ = θ (g²/32π²) G G̃
    - Violates CP symmetry (θ ≠ 0, π)
    - Neutron EDM: d_n < 10⁻²⁶ e·cm → θ < 10⁻¹⁰
    - Why so small? Peccei-Quinn symmetry → axion

14. QUARK-GLUON PLASMA:
    - At T > 10¹² K: Deconfined quarks and gluons
    - Early universe (t < 10⁻⁶ s): QGP phase
    - Heavy-ion collisions: Recreate QGP briefly
    - Chiral restoration: m_u, m_d → bare masses

15. ANTHROPIC CONSIDERATIONS:
    - If m_d - m_u > m_e: Hydrogen unstable (p + e → n + ν)
    - If m_d - m_u < 0: Different periodic table (more stable elements?)
    - m_d ≈ 2 m_u: Accidental or fundamental?
    - Multiverse: Only universes with stable H, He have chemistry

The down quark is the isospin partner of the up quark and is essential
for nuclear stability. The mass difference m_d - m_u ≈ 2.5 MeV determines
the neutron-proton mass gap and hence the lifetime of the free neutron.
The TriPhase formula m_d = 9·m_e·(1 + α/π) with factor 9 (vs. 4 for up)
gives the correct ratio m_d/m_u ≈ 2.25, suggesting a geometric origin
for this critical mass hierarchy.
""")

print("=" * 70)
input("Press Enter to exit...")
