"""
TriPhase V16 — Cosmological Horizon (18-Step) (Statistical Mechanics Framework)
================================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
The cosmological horizon is the maximum comoving distance from which light could
have reached us since the Big Bang, defining the boundary of the observable universe.
In statistical mechanics, this represents the size of the "system" — the region over
which we can compute partition functions and ensemble averages. The Hubble radius
r_H = c/H_0 ≈ 4.4 Gpc marks the scale where cosmic expansion equals the speed of
light. Objects beyond this horizon are causally disconnected from us; their microstates
cannot affect measurements in our cosmic volume. The partition function of the
universe is thus computed over the Hubble volume V_H ~ (c/H_0)³.

In the grand canonical ensemble for cosmology, the horizon sets the 'box size' for
statistical mechanics. The entropy within the horizon scales as S ~ (c/H_0)² / l_P²,
where l_P is the Planck length — this is the holographic bound, suggesting that
entropy scales with area rather than volume. The TriPhase derivation connects the
horizon to fundamental constants via H_0 = π√3 f_e α^18, encoding the universe's
expansion rate through the electron frequency and fine structure constant. This
links cosmological scales to atomic physics through 18 powers of alpha, suggesting
a deep connection between quantum mechanics and cosmic structure.

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
print("TriPhase V16: Cosmological Horizon (Statistical Mechanics)")
print("=" * 70)
print()

print("STATISTICAL MECHANICS FRAMEWORK")
print("--------------------------------")
print("Ensemble: Grand Canonical (entire observable universe)")
print("System size: Hubble volume V_H ~ (c/H_0)³")
print("Boundary: Causal horizon r_H = c/H_0")
print("Observable: Holographic entropy S ~ A_H / (4 l_P²)")
print()

print("HUBBLE RADIUS AS SYSTEM BOUNDARY")
print("---------------------------------")
print(f"Speed of light c = {c:.6e} m/s")
print(f"Electron frequency f_e = {f_e:.6e} Hz")
print(f"Fine structure α = {alpha:.10f}")
print()

# TriPhase Hubble constant
print("TriPhase Hubble constant:")
print(f"  H_0 = π√3 × f_e × α^18")
print()

H_0_calc = math.pi * math.sqrt(3.0) * f_e * alpha**18
print(f"  H_0 (TriPhase) = {H_0_calc:.6e} Hz")
print()

# Convert to km/s/Mpc (standard cosmology units)
# 1 Mpc = 3.0857e22 m
Mpc = 3.0857e22  # meters
H_0_cosmo = H_0_calc * Mpc / 1000.0  # Convert to km/s/Mpc
print(f"  H_0 (TriPhase) = {H_0_cosmo:.4f} km/s/Mpc")
print()

# Hubble radius (cosmological horizon)
r_H = c / H_0_calc
r_H_Gpc = r_H / (1e9 * Mpc)  # Convert to Gpc
print(f"Hubble radius r_H = c / H_0")
print(f"  r_H = {r_H:.6e} m")
print(f"  r_H = {r_H_Gpc:.4f} Gpc")
print()

# Hubble volume and holographic entropy
V_H = (4.0 / 3.0) * math.pi * r_H**3
l_P = math.sqrt(hbar * G / c**3)  # Planck length
A_H = 4.0 * math.pi * r_H**2
S_H = A_H / (4.0 * l_P**2)

print(f"Hubble volume V_H = {V_H:.6e} m³")
print(f"Hubble area A_H = {A_H:.6e} m²")
print(f"Planck length l_P = {l_P:.6e} m")
print(f"Holographic entropy S_H = A_H/(4l_P²) = {S_H:.6e}")
print(f"  S_H / k_B = {S_H / 1.380649e-23:.6e} (in bits ~ {math.log2(math.exp(S_H / 1.380649e-23)):.1e})")
print()

# ========== CALIBRATION CHECKPOINT ==========
print("CALIBRATION CHECKPOINT")
print("----------------------")
# Planck 2018: H_0 = 67.36 ± 0.54 km/s/Mpc (CMB)
# SH0ES 2021: H_0 = 73.04 ± 1.04 km/s/Mpc (local distance ladder)
H_0_planck = 67.36
H_0_shoes = 73.04

deviation_planck = (H_0_cosmo - H_0_planck) / H_0_planck * 1e6
deviation_shoes = (H_0_cosmo - H_0_shoes) / H_0_shoes * 1e6

print(f"Planck 2018 (CMB):       {H_0_planck:.4f} km/s/Mpc")
print(f"SH0ES 2021 (local):      {H_0_shoes:.4f} km/s/Mpc")
print(f"TriPhase:                {H_0_cosmo:.4f} km/s/Mpc")
print()
print(f"Deviation from Planck:   {deviation_planck:.0f} ppm")
print(f"Deviation from SH0ES:    {deviation_shoes:.0f} ppm")
print()
print("Note: Hubble tension (Planck vs SH0ES) ≈ 8% discrepancy remains unresolved.")
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT")
print("-----------------------------")
print("The cosmological horizon defines the maximum system size for statistical")
print("mechanics in an expanding universe. Objects beyond r_H = c/H_0 are receding")
print("faster than light (locally allowed by general relativity) and cannot exchange")
print("information with us. The partition function Z must therefore be computed only")
print("over the Hubble volume — a finite region containing ~10^80 baryons.")
print()
print("The holographic principle suggests entropy is bounded by horizon area:")
print("  S_max = A_H / (4 l_P²) ~ 10^123 k_B")
print()
print("This implies the universe's information content scales as area, not volume —")
print("a profound hint that spacetime may be emergent from lower-dimensional degrees")
print("of freedom. In AdS/CFT correspondence, the bulk partition function equals a")
print("boundary CFT partition function, making this correspondence precise.")
print()
print("The TriPhase formula H_0 = π√3 f_e α^18 connects cosmic expansion to atomic")
print("physics. The 18th power of alpha (≈ 10^-40) naturally generates the enormous")
print("ratio between atomic scales (1/f_e ~ 10^-21 s) and cosmological timescales")
print("(1/H_0 ~ 10^10 yr). This suggests the universe's expansion may be linked to")
print("electromagnetic vacuum structure — a tantalizing hint of deeper unity between")
print("quantum field theory and general relativity.")
print()
print("The Hubble tension (8% discrepancy between CMB and local measurements) may")
print("indicate new physics: early dark energy, varying constants, or systematic")
print("errors. Statistical mechanics provides a framework for testing these scenarios")
print("through ensemble averages over cosmic structure formation simulations.")
print()
print("=" * 70)

input("Press Enter to exit...")
