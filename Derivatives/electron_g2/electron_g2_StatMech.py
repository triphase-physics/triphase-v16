"""
TriPhase V16 — Electron Anomalous Magnetic Moment g-2 (Statistical Mechanics Framework)
========================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
The electron's anomalous magnetic moment (g-2) quantifies the deviation from the
Dirac prediction g = 2. In QED, this anomaly arises from virtual photon loops that
contribute to the electron's magnetic moment through vacuum fluctuations. The
statistical mechanics interpretation is that (g-2)/2 represents the correction to
the magnetic moment from the grand canonical ensemble of virtual photon states.

The lowest-order QED correction is the Schwinger term: a_e = α/(2π) ≈ 0.001162.
Higher orders include contributions from multiple photon loops, weighted by powers
of α. The partition function for the electron-photon system includes all loop
diagrams: Z = exp(-S_eff), where S_eff contains the renormalized moment.

From the statistical perspective, (g-2)/2 measures the fluctuation-induced shift
in the magnetic moment. Virtual e⁺e⁻ pairs screen the electron's magnetic field,
similar to charge screening in a plasma. The factor α/(2π) emerges from integrating
over the phase space of virtual photons in the canonical ensemble.

TAG: (D*) — TriPhase prediction including QED radiative corrections
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
print("TriPhase V16: Electron g-2 Anomaly (Statistical Mechanics)")
print("=" * 70)
print()

print("DIRAC PREDICTION:")
print("-" * 70)
print("The Dirac equation predicts g = 2 exactly for a point-like spin-1/2 fermion.")
print("The magnetic moment is:")
print("  μ = g·(e/2m_e)·S = 2·μ_B  (where μ_B = eℏ/2m_e)")
print()

mu_B = e * hbar / (2.0 * m_e)  # Bohr magneton
print(f"Bohr magneton:  μ_B = eℏ/(2m_e) = {mu_B:.6e} J/T")
print()

print("QED CORRECTION — SCHWINGER TERM:")
print("-" * 70)
print("Vacuum fluctuations (virtual photon loops) modify g:")
print("  g = 2(1 + a_e)")
print("  where a_e = (g-2)/2 is the anomalous moment")
print()
print("Lowest order (one-loop):")
print()

a_e_Schwinger = alpha / (2.0 * math.pi)

print(f"  a_e^(1) = α/(2π) = {a_e_Schwinger:.10f}")
print()

print("This is the famous Schwinger result (1948)—the first QED prediction.")
print()

print("HIGHER-ORDER QED:")
print("-" * 70)
print("The full QED expansion is:")
print("  a_e = C₁(α/π) + C₂(α/π)² + C₃(α/π)³ + ...")
print()
print("where C₁ = 1/2, C₂ = -0.32848..., C₃ = 1.181..., etc.")
print()

# Simplified higher-order estimate (not exact, just illustrative)
a_e_full = a_e_Schwinger  # We only compute the leading term here

print(f"Leading term only:  a_e ≈ {a_e_full:.10f}")
print()
print("(Full QED calculation includes O(α²), O(α³), ..., O(α⁵) terms,")
print("plus hadronic and weak contributions. Not computed here.)")
print()

print("STATISTICAL MECHANICS INTERPRETATION:")
print("-" * 70)
print("The anomalous moment arises from vacuum fluctuations in the canonical")
print("ensemble of virtual photons.")
print()
print("The partition function is:")
print("  Z = ∫ D[A_μ] exp(-S_QED/ℏ)")
print()
print("The magnetic moment is:")
print("  μ = ⟨ψ̄ σ_μν F^μν ψ⟩_Z / ⟨ψ̄ψ⟩_Z")
print()
print("Virtual photon loops contribute to the numerator:")
print("  Δμ ~ ∫ d⁴k/(k²) · (e²/k²) ~ α/(2π)")
print()
print("The integral over photon momenta k gives the factor 1/(2π) from")
print("phase space, and α from the coupling.")
print()

print("FLUCTUATION-DISSIPATION INTERPRETATION:")
print("-" * 70)
print("The anomalous moment measures the magnetic susceptibility of the vacuum.")
print()
print("  χ_mag = ∂μ/∂B|_B=0 ~ a_e · μ_B")
print()
print("This susceptibility arises from vacuum polarization—virtual pairs align")
print("with the applied field, screening the electron's intrinsic moment.")
print()
print("The factor α/(2π) is the statistical weight for creating a virtual pair")
print("that contributes to screening.")
print()

# ========== CALIBRATION CHECKPOINT ==========
a_e_measured = 0.00115965218073  # CODATA 2018 (experimental)
a_e_QED_full = 0.00115965218178  # Full QED theory (5-loop + hadronic)
a_e_calc = a_e_Schwinger  # Leading order only

deviation_measured = (a_e_calc - a_e_measured) / a_e_measured * 1e6
deviation_QED = (a_e_calc - a_e_QED_full) / a_e_QED_full * 1e6

print("=" * 70)
print("CALIBRATION CHECKPOINT")
print("=" * 70)
print(f"Measured (2018):        a_e = {a_e_measured:.14f}")
print(f"QED theory (5-loop):        = {a_e_QED_full:.14f}")
print(f"Schwinger (1-loop only):    = {a_e_calc:.10f}")
print()
print(f"Deviation from measurement:   {deviation_measured:+.0f} ppm")
print(f"Deviation from full QED:      {deviation_QED:+.0f} ppm")
print()
print("(The Schwinger term is only ~0.1% of a_e. Higher orders are needed")
print("for precision comparison.)")
print("=" * 70)
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT:")
print("-" * 70)
print("The electron's anomalous magnetic moment is a window into vacuum")
print("fluctuations. The value a_e = α/(2π) emerges from the statistical")
print("ensemble of virtual photons that 'dress' the electron.")
print()
print("From the partition function perspective:")
print("  Z_bare = exp(-S_Dirac)  →  g = 2  (no fluctuations)")
print("  Z_QED = exp(-S_QED)     →  g = 2(1 + α/(2π))  (1-loop)")
print()
print("Each additional loop contributes a factor (α/π)^n, weighted by a")
print("numerical coefficient C_n that comes from integrating over the phase")
print("space of virtual particles.")
print()
print("The factor 1/(2π) is universal—it appears in all 1-loop corrections")
print("in QFT. It represents the density of states in momentum space:")
print("  ∫ d⁴k/(2π)⁴ (relativistic phase space volume)")
print()
print("The electron g-2 is the most precisely measured quantity in physics")
print("(parts per trillion). The agreement between theory and experiment")
print("validates QED as the correct statistical mechanics of the EM vacuum.")
print()
print("This is fluctuation-dissipation in action: the vacuum's response to")
print("an external magnetic field (dissipation) is determined by its intrinsic")
print("fluctuations (virtual pairs). The anomalous moment a_e = α/(2π) is")
print("the quantitative link between these two aspects of vacuum statistics.")
print("=" * 70)

input("Press Enter to exit...")
