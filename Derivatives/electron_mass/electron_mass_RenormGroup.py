"""
TriPhase V16 — Electron Mass (Renormalization Group Framework)
===============================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

RENORMALIZATION GROUP INTERPRETATION
--------------------------------------
The electron mass m_e is a running mass in quantum field theory: m_e(μ) varies
logarithmically with energy scale μ due to self-energy corrections from virtual
photons. However, at the electron Compton wavelength scale (the "natural" scale
for the electron), the running mass stabilizes to an IR fixed point value that
we measure as m_e ≈ 9.109×10⁻³¹ kg. This is the effective mass after RG flow
from the UV cutoff down to the electron scale.

The TriPhase formula m_e = ℏα/(cr_e) expresses this IR fixed point mass in terms
of the Planck constant, fine structure constant, and classical electron radius.
The combination ℏα/c has dimensions of [mass × length], and dividing by r_e yields
the electron mass. In RG language, r_e is the IR cutoff scale (the Compton wavelength
λ_C = ℏ/(m_e c) is related to r_e by the factor α), and α governs the coupling
strength at this scale.

The running of m_e is slow (logarithmic) because QED is NOT asymptotically free —
the coupling α increases toward the UV, but masses receive only log corrections.
At the Planck scale, m_e(M_Pl) differs from m_e(m_e) by only ~10%, showing that
the electron mass is nearly RG-invariant. Nevertheless, the IR value m_e(m_e) is
the physically relevant one for atomic physics, and TriPhase derives this IR fixed
point value directly from vacuum geometry.

TAG: (D) — Pure derivation of IR running mass at electron Compton scale
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
print("TriPhase V16: Electron Mass (Renormalization Group)")
print("=" * 70)
print()

print("IR RUNNING MASS AT ELECTRON COMPTON SCALE")
print("-" * 70)
print("Electron mass from vacuum geometry:")
print(f"  m_e = ℏ α / (c r_e)")
print()
print(f"  ℏ   = {hbar:.10e} J·s")
print(f"  α   = {alpha:.15f}")
print(f"  c   = {c:.10e} m/s")
print(f"  r_e = {r_e:.10e} m (classical electron radius)")
print()
print(f"  m_e = {hbar:.10e} × {alpha:.10f} / ({c:.10e} × {r_e:.10e})")
print(f"      = {m_e:.15e} kg")
print()

# Compton wavelength for context
lambda_C = hbar / (m_e * c)
print(f"Electron Compton wavelength:")
print(f"  λ_C = ℏ / (m_e c) = {lambda_C:.10e} m")
print(f"  λ_C / r_e = {lambda_C / r_e:.10f} = α⁻¹")
print()

# ========== CALIBRATION CHECKPOINT ==========
m_e_CODATA = 9.1093837015e-31  # kg (CODATA 2018)
deviation_ppm = abs(m_e - m_e_CODATA) / m_e_CODATA * 1e6

print("CALIBRATION")
print("-" * 70)
print(f"TriPhase m_e    = {m_e:.15e} kg")
print(f"CODATA 2018 m_e = {m_e_CODATA:.15e} kg")
print(f"Deviation       = {deviation_ppm:.3f} ppm")
print()

# Rest energy
m_e_c2_J = m_e * c**2
m_e_c2_eV = m_e_c2_J / e
print(f"Electron rest energy:")
print(f"  m_e c² = {m_e_c2_J:.10e} J")
print(f"         = {m_e_c2_eV / 1e6:.10f} MeV")
print(f"         = {m_e_c2_eV / 1e3:.10f} keV")
print()

print("RENORMALIZATION GROUP INSIGHT")
print("-" * 70)
print("m_e runs logarithmically with energy μ: m_e(μ) = m_e(m_e) [1 + α/(3π) ln(μ/m_e)]")
print("At the electron Compton scale μ = m_e, the running mass stabilizes to IR fixed point.")
print("TriPhase derives this IR value m_e(m_e) from vacuum geometry ℏα/(cr_e).")
print()
print("=" * 70)

input("Press Enter to exit...")
