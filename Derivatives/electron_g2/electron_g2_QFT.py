"""
TriPhase V16 - Electron Anomalous Magnetic Moment (QFT Framework)
==================================================================

QFT INTERPRETATION:
The electron g-factor g ≈ 2.002319... quantifies the magnetic moment μ = gμ_B s/ℏ:
- Dirac equation prediction: g = 2 exactly for point-like spin-1/2 particle
- QED radiative corrections: Δg = g - 2 from virtual photon loops
- Schwinger's 1-loop result: (g-2)/2 = α/(2π) ≈ 0.00116
- Higher loops: 2-loop, 3-loop, ... up to 5-loop calculations
- Most precisely measured quantity in physics: tests QED to 10 parts per trillion

TriPhase's formula g = 2(1 + α/(2π) - 0.328(α/π)²) includes the 1-loop and
approximate 2-loop corrections, demonstrating how perturbative QED builds up
the anomalous magnetic moment from vacuum polarization and vertex corrections.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
TAG: (D) - Direct derivation with QED loop corrections
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

# ========== QFT DERIVATION: ELECTRON g-2 ==========
print("=" * 70)
print("TriPhase V16 - Electron Anomalous Magnetic Moment (g-2)")
print("QFT Framework: Perturbative QED Loop Corrections")
print("=" * 70)
print()

print("QFT CONTEXT:")
print("The electron g-factor measures the magnetic moment in units of the Bohr")
print("magneton μ_B = eℏ/(2m_e). Dirac theory predicts g = 2 exactly, but QED")
print("corrections from virtual particle loops modify this:")
print()
print("  Vertex diagrams:    e⁻ → e⁻ + γ_virtual (1-loop Schwinger term)")
print("  Vacuum polarization: photon → e⁺e⁻ pair (modifies photon propagator)")
print("  Light-by-light:     γγ → γγ scattering (4-loop contribution)")
print()
print("The anomalous moment a_e = (g-2)/2 is the most stringent test of QED.")
print()

print("TRIPHASE DERIVATION:")
print("g = 2 × (1 + α/(2π) - 0.328(α/π)²)")
print()
print(f"Fine structure:       α = {alpha:.12f}")
print(f"1-loop (Schwinger):   α/(2π) = {alpha/(2.0*math.pi):.10e}")
print(f"2-loop term:          0.328(α/π)² = {0.328*(alpha/math.pi)**2:.10e}")
print(f"Correction sum:       {alpha/(2.0*math.pi) - 0.328*(alpha/math.pi)**2:.10e}")
print()

g_e = 2.0 * (1.0 + alpha/(2.0*math.pi) - 0.328*(alpha/math.pi)**2)
a_e = (g_e - 2.0) / 2.0

print(f"g (TriPhase):         {g_e:.14f}")
print(f"a_e = (g-2)/2:        {a_e:.10e}")
print()

# ========== CALIBRATION CHECKPOINT ==========
codata_g = 2.00231930436256  # CODATA 2018
codata_a_e = (codata_g - 2.0) / 2.0
deviation_ppm = (g_e - codata_g) / codata_g * 1e6

print("CALIBRATION CHECKPOINT:")
print(f"CODATA 2018 g:        {codata_g:.14f}")
print(f"TriPhase g:           {g_e:.14f}")
print(f"Deviation:            {deviation_ppm:+.2f} ppm")
print()
print(f"CODATA a_e:           {codata_a_e:.10e}")
print(f"TriPhase a_e:         {a_e:.10e}")
print(f"Δa_e:                 {a_e - codata_a_e:.10e}")
print()

# Show perturbative expansion
print("PERTURBATIVE EXPANSION:")
schwinger = alpha / (2.0 * math.pi)
two_loop = 0.328 * (alpha / math.pi)**2
print(f"O(α):   Schwinger     {schwinger:.10e}")
print(f"O(α²):  2-loop        {-two_loop:.10e} (negative)")
print(f"Sum:                  {schwinger - two_loop:.10e}")
print()
print("Note: Full CODATA includes O(α³), O(α⁴), O(α⁵) terms and hadronic/weak")
print("contributions. TriPhase captures the dominant QED structure.")
print()

# ========== QFT INSIGHT ==========
print("QFT INSIGHT:")
print("The Schwinger term α/(2π) arises from the simplest 1-loop vertex correction:")
print()
print("      e⁻ ────┬──── e⁻")
print("              │ γ_virtual")
print("              └─────┘")
print()
print("This virtual photon loop modifies the electron-photon vertex, changing the")
print("effective magnetic moment. Higher loops add corrections: 2-loop (α/π)² includes")
print("vacuum polarization insertions and multi-photon vertices. The factor 0.328")
print("comes from detailed Feynman diagram integration.")
print()
print("The extraordinary agreement (sub-ppm level) between QED calculations and")
print("experimental measurements validates the entire renormalization program and")
print("confirms the reality of virtual particle clouds around 'bare' particles.")
print()
print("=" * 70)

input("Press Enter to exit...")
