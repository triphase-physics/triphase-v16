"""
========================================================================
TriPhase V16 Derivative: Electron g-Factor (Gauge Theory)
========================================================================

GAUGE THEORY INTERPRETATION:
The electron g-factor characterizes the magnetic moment μ = (g/2)(eℏ/2m_e)S,
where S is the spin angular momentum. The Dirac equation predicts g = 2
exactly, but quantum electrodynamics (QED) corrections from virtual photon
loops modify this to g ≈ 2.002319.

In gauge theory, the anomalous magnetic moment a_e = (g-2)/2 arises from
radiative corrections to the electron-photon vertex. Each order in the
U(1) gauge coupling α contributes additional Feynman diagrams: one-loop
(Schwinger term ~ α/2π), two-loop (~ α²), three-loop (~ α³), etc.

The TriPhase formula g = 2(1 + α/2π - 0.328(α/π)²) captures the leading
QED corrections. The coefficient 0.328 ≈ 1/3 represents the sum of all
two-loop diagrams involving virtual electron-positron pairs and photon
propagators. This is one of the most precisely calculated and measured
quantities in physics, testing QED to extraordinary precision.

REFERENCE: Experimental g/2 - 1 = 0.001159652181643(764)

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
TAG: (D)
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
print("GAUGE THEORY DERIVATION: Electron g-Factor (QED Corrections)")
print("=" * 70)

# Derive g-2 from QED perturbative expansion
g2 = 2.0 * (1.0 + alpha/(2.0*math.pi) - 0.328*(alpha/math.pi)**2)

print("\nQED Radiative Corrections to Magnetic Moment:")
print(f"  Dirac value:                 g₀ = 2.0 (tree level)")
print(f"  α (EM gauge coupling):       {alpha:.15f}")
print(f"  One-loop (Schwinger):        α/(2π) = {alpha/(2.0*math.pi):.15e}")
print(f"  Two-loop coefficient:        0.328")
print(f"  Two-loop term:               -0.328(α/π)² = {-0.328*(alpha/math.pi)**2:.15e}")
print(f"  g = 2(1 + α/2π - 0.328(α/π)²): {g2:.15f}")

# Anomalous magnetic moment
a_e = (g2 - 2.0) / 2.0

print(f"\nAnomalous magnetic moment:")
print(f"  a_e = (g-2)/2:               {a_e:.15e}")

# Individual contributions
a_e_1loop = alpha / (2.0 * math.pi)
a_e_2loop = -0.328 * (alpha / math.pi)**2

print(f"\nContributions by order:")
print(f"  O(α¹):   {a_e_1loop:.15e}")
print(f"  O(α²):   {a_e_2loop:.15e}")
print(f"  Total:   {a_e:.15e}")

print("\n" + "=" * 70)
print("CALIBRATION CHECKPOINT")
print("=" * 70)

# Experimental value (Harvard 2008, Hanneke et al.)
a_e_expt = 0.001159652181643
a_e_uncertainty = 0.000000000000764

deviation = abs(a_e - a_e_expt)
deviation_sigma = deviation / a_e_uncertainty

print(f"\nTriPhase a_e:     {a_e:.15e}")
print(f"Experimental:     {a_e_expt:.15e}")
print(f"Uncertainty:      ±{a_e_uncertainty:.15e}")
print(f"Deviation:        {deviation:.15e}")
print(f"Deviation (σ):    {deviation_sigma:.3f} σ")

if deviation_sigma < 3:
    print("✓ Within 3σ of experimental value")
elif deviation_sigma < 5:
    print("✓ Within 5σ (possible tension)")
else:
    print("⚠ Exceeds 5σ (significant discrepancy)")

print("\n" + "=" * 70)
print("GAUGE THEORY INSIGHTS")
print("=" * 70)

print("""
The electron g-factor in gauge theory:

1. DIRAC EQUATION PREDICTION:
   - Relativistic quantum mechanics: g = 2 (exact)
   - Magnetic moment: μ = (g/2) μ_B S (μ_B = Bohr magneton)
   - Interaction: H = -μ·B = -(g/2) μ_B S·B
   - No quantum field theory corrections at tree level

2. QED RADIATIVE CORRECTIONS:
   - Virtual photon loops modify electron-photon vertex
   - Anomalous moment: a_e = (g-2)/2
   - Perturbative expansion: a_e = C₂(α/π) + C₄(α/π)² + C₆(α/π)³ + ...
   - Calculable to all orders in principle

3. SCHWINGER TERM (ONE-LOOP):
   - Julian Schwinger (1948): First QED correction
   - C₂ = 1/2 → a_e = α/(2π) ≈ 0.00116
   - Single Feynman diagram: Electron emits and reabsorbs photon
   - Vertex correction modifies magnetic coupling

4. TWO-LOOP CORRECTIONS:
   - Seven Feynman diagrams at O(α²)
   - Electron self-energy insertions
   - Vacuum polarization (virtual e⁺e⁻ pairs)
   - Light-by-light scattering
   - Sum: C₄ ≈ -0.328 (TriPhase) vs. exact -0.32847...

5. HIGHER-ORDER TERMS:
   - Three-loop: 891 diagrams, C₆ ≈ 1.181
   - Four-loop: ~12,672 diagrams, C₈ ≈ -1.909
   - Five-loop: ~10⁷ diagrams, C₁₀ ≈ 9.16
   - Computed numerically to incredible precision

6. HADRONIC CONTRIBUTIONS:
   - Virtual quark loops (not just leptons)
   - Hadronic vacuum polarization: Δa_e^had ≈ 1.7×10⁻¹²
   - Hadronic light-by-light: Δa_e^hll ≈ 0.03×10⁻¹²
   - Small but measurable at current precision

7. WEAK INTERACTION CORRECTIONS:
   - Virtual W, Z boson loops
   - Δa_e^weak ≈ 0.03×10⁻¹² (tiny, m_e << M_W)
   - Grows for heavier leptons (muon, tau)

8. EXPERIMENTAL MEASUREMENT:
   - Penning trap: Single electron in magnetic field
   - Cyclotron frequency ω_c = eB/m_e
   - Spin precession frequency ω_s = (g/2) ω_c
   - Measure ω_s/ω_c to 0.28 parts per trillion!

9. FINE-STRUCTURE CONSTANT DETERMINATION:
   - QED theory: a_e = f(α) (known function)
   - Experiment: a_e = measured value
   - Invert: α⁻¹ = 137.035999084(51) (most precise α determination!)
   - Tests QED self-consistency

10. COMPARISON TO MUON g-2:
    - Muon anomalous moment: a_μ = (g_μ - 2)/2
    - Same QED diagrams, but m_μ = 207 m_e matters
    - Hadronic contributions ~630×10⁻¹⁰ (much larger!)
    - Current tension: 4-5σ between theory and experiment

11. GAUGE INVARIANCE:
    - g-factor is physical observable (gauge-independent)
    - Vertex function Γ_μ(p,p') encodes all corrections
    - Ward-Takahashi identity: Ensures gauge invariance
    - Anomalous moment = gauge-invariant residue

12. TESTING BSM PHYSICS:
    - New particles (SUSY, dark matter) contribute to loops
    - Precision measurement constrains new physics
    - Current agreement: Strong evidence for SM at low energy
    - Any deviation would signal new gauge bosons or fermions

The electron g-2 is the most precisely tested prediction in all of
physics. The agreement between QED theory and experiment to 12 decimal
places is a triumph of quantum field theory and gauge invariance. The
TriPhase formula captures the first two orders (α¹, α²) and demonstrates
how electromagnetic gauge coupling generates observable corrections to
the classical Dirac value.
""")

print("=" * 70)
input("Press Enter to exit...")
