"""
TriPhase V16 — Electron Anomalous Magnetic Moment g-2 (Renormalization Group Framework)
=========================================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

RENORMALIZATION GROUP INTERPRETATION
--------------------------------------
The electron anomalous magnetic moment a_e = (g-2)/2 is THE prototypical example of
renormalization group corrections in QED. The Dirac equation predicts g = 2 exactly,
but radiative corrections from virtual photons and electron-positron loops shift g
away from 2. The formula a_e = α/(2π) - (α/π)²×0.328478965... is a loop expansion
in powers of α: the first term α/(2π) is the one-loop Schwinger correction, and
the second term is the two-loop correction.

In RG language, g-2 IS the anomalous dimension of the electron's magnetic moment
operator. As you integrate out high-energy virtual photons shell by shell (Wilson's
RG procedure), the effective magnetic moment receives corrections at each loop order.
The coefficient α/(2π) comes from the one-loop vertex correction diagram (Schwinger 1948),
while the two-loop coefficient 0.328478965... includes vacuum polarization and vertex
corrections at O(α²).

This is pure RG flow: each loop order corresponds to integrating out one more shell
of virtual momenta. The fact that higher-loop corrections are suppressed by powers
of α ~ 1/137 shows that QED is a weakly coupled theory in the IR (perturbative RG).
The agreement between theory (computed to 5 loops) and experiment (measured to 12
digits) is the most precise test of QED and RG methods in all of physics.

TAG: (D) — Pure derivation (RG loop expansion in α)
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
print("TriPhase V16: Electron g-2 Anomaly (Renormalization Group)")
print("=" * 70)
print()

print("RG LOOP EXPANSION: ANOMALOUS DIMENSION OF MAGNETIC MOMENT")
print("-" * 70)
print("Dirac equation (tree level): g = 2 exactly")
print()
print("One-loop correction (Schwinger 1948):")
print(f"  a_e^(1) = α / (2π)")
print(f"          = {alpha:.15f} / (2π)")
print(f"          = {alpha / (2.0 * math.pi):.15e}")
print()

a_e_1loop = alpha / (2.0 * math.pi)

print("Two-loop correction (Petermann 1957, Sommerfield 1957):")
C_2 = 0.328478965  # Two-loop coefficient
print(f"  a_e^(2) = -(α/π)² × {C_2}")
print(f"          = -({alpha:.10f}/π)² × {C_2}")
print(f"          = {-(alpha / math.pi)**2 * C_2:.15e}")
print()

a_e_2loop = -(alpha / math.pi)**2 * C_2

print("Total anomalous magnetic moment (up to two loops):")
print(f"  a_e = α/(2π) - (α/π)² × {C_2}")
print(f"      = {a_e_1loop:.15e} + {a_e_2loop:.15e}")
print(f"      = {a_e_1loop + a_e_2loop:.15e}")
print()

a_e = a_e_1loop + a_e_2loop

# ========== CALIBRATION CHECKPOINT ==========
a_e_CODATA = 1.15965218128e-3  # CODATA 2018 (includes up to 5-loop QED)
deviation_ppb = abs(a_e - a_e_CODATA) / a_e_CODATA * 1e9

print("CALIBRATION")
print("-" * 70)
print(f"TriPhase a_e (2-loop)     = {a_e:.15e}")
print(f"CODATA 2018 a_e (5-loop)  = {a_e_CODATA:.15e}")
print(f"Deviation                 = {deviation_ppb:.3f} ppb")
print()
print("Note: CODATA includes 3rd, 4th, 5th loop corrections (~10⁻¹² contribution),")
print("      plus hadronic and weak corrections. TriPhase 2-loop is within 0.1%.")
print()

# Higher loop structure (conceptual)
print("HIGHER-LOOP RG STRUCTURE")
print("-" * 70)
print("Full QED expansion (5-loop, Aoyama et al. 2012):")
print(f"  a_e = C₁(α/π) + C₂(α/π)² + C₃(α/π)³ + C₄(α/π)⁴ + C₅(α/π)⁵ + ...")
print()
print(f"  C₁ = 0.5 (Schwinger)")
print(f"  C₂ = {C_2} (Petermann, Sommerfield)")
print(f"  C₃ ≈ 1.181... (3-loop, 891 Feynman diagrams)")
print(f"  C₄ ≈ -1.909... (4-loop, 12,672 diagrams)")
print(f"  C₅ ≈ 9.16... (5-loop, >10⁷ diagrams)")
print()
print("Each loop = one RG step, integrating out virtual photons at higher momenta.")
print()

print("RENORMALIZATION GROUP INSIGHT")
print("-" * 70)
print("g-2 IS the anomalous dimension: RG corrections from virtual photon loops.")
print("Each power of α corresponds to one loop (one RG shell integration).")
print("This is the most precise test of QED and RG methods: theory vs experiment to 12 digits.")
print()
print("=" * 70)

input("Press Enter to exit...")
