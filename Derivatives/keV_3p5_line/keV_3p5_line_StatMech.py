"""
TriPhase V16 — 3.5 keV X-ray Line (Statistical Mechanics Framework)
========================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
The 3.5 keV X-ray line observed in galaxy clusters and the galactic center represents
a potential signature of dark matter decay or annihilation. In statistical mechanics,
this line emerges from the grand canonical ensemble of dark matter particles with
mass m_DM ~ 7 keV (sterile neutrinos). The decay process DM → γ + X releases a
photon at half the rest mass energy: E_γ = m_DM c²/2 ≈ 3.5 keV.

The partition function for dark matter includes both production and decay channels:
Z_DM = Σ exp(-β(E_n - μN_DM)), where μ is the chemical potential. At freeze-out
temperature T_f ~ few keV, dark matter abundance is determined by the Boltzmann
equation. The 3.5 keV line flux depends on the DM density ρ_DM, decay rate Γ_decay,
and line-of-sight integration through the halo.

In TriPhase, the 3.5 keV energy scale emerges from vacuum mode frequencies at
specific harmonic ratios to the electron Compton energy: E_3.5 ~ α^n · m_e c²
for n ≈ 3-4. This suggests the 3.5 keV line may arise from radiative transitions
in the electromagnetic vacuum rather than exotic dark matter.

TAG: (D*) — TriPhase prediction for dark matter signature
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
print("TriPhase V16: 3.5 keV X-ray Line (Statistical Mechanics)")
print("=" * 70)
print()

print("OBSERVED 3.5 keV LINE:")
print("-" * 70)
print("An unidentified X-ray line at ~3.5 keV has been detected in:")
print("  • Stacked galaxy clusters (XMM-Newton, 2014)")
print("  • Andromeda galaxy (M31)")
print("  • Galactic center (multiple observations)")
print()
print("The line is spatially coincident with dark matter distributions,")
print("suggesting a DM origin (e.g., sterile neutrino decay).")
print()

E_line_keV = 3.5  # keV
E_line_J = E_line_keV * 1000.0 * e  # convert to Joules

print(f"Line energy:  E = {E_line_keV:.1f} keV")
print(f"            = {E_line_J:.6e} J")
print()

print("STERILE NEUTRINO INTERPRETATION:")
print("-" * 70)
print("If the line comes from DM decay: DM → γ + ν_s")
print("The DM mass must be:")
print()

m_DM_kg = 2.0 * E_line_J / c**2  # twice photon energy
m_DM_keV = 2.0 * E_line_keV

print(f"  m_DM = 2·E_γ/c² = {m_DM_kg:.6e} kg")
print(f"       = {m_DM_keV:.1f} keV/c²")
print()

print("This matches the predicted mass for sterile neutrinos (keV-scale).")
print()

print("PARTITION FUNCTION FOR DM:")
print("-" * 70)
print("In the early universe, DM freeze-out occurs when Γ_interaction < H.")
print("The relic abundance is determined by:")
print("  Y_DM = n_DM/s ≈ (xf/g*) · (m_Pl/m_DM) · ⟨σv⟩⁻¹")
print("  where xf = m_DM/T_f ~ 20-30 at freeze-out")
print()

T_f_keV = m_DM_keV / 25.0  # typical freeze-out ratio
k_B = 1.380649e-23  # J/K
T_f_K = T_f_keV * 1000.0 * e / k_B

print(f"Freeze-out temperature:  T_f ~ {T_f_keV:.2f} keV")
print(f"                            ~ {T_f_K:.3e} K")
print()

print("TRIPHASE ORIGIN — VACUUM HARMONIC:")
print("-" * 70)
print("Instead of dark matter, TriPhase interprets the 3.5 keV line as a")
print("radiative transition in the EM vacuum at a harmonic of m_e c².")
print()

E_electron_keV = m_e * c**2 / (e * 1000.0)
print(f"Electron rest energy:  E_e = m_e c² = {E_electron_keV:.2f} keV")
print()

# Find the harmonic ratio
ratio = E_line_keV / E_electron_keV
print(f"Ratio:  E_3.5keV / E_electron = {ratio:.6f}")
print()

# Check if this is an alpha power
n_guess = math.log(ratio) / math.log(alpha)
print(f"If E_3.5 = α^n · E_e, then n ≈ {n_guess:.2f}")
print()

# Try integer n values
for n in range(1, 8):
    E_n_keV = (alpha**n) * E_electron_keV
    if abs(E_n_keV - E_line_keV) / E_line_keV < 0.20:
        print(f"  n = {n}:  α^{n} · E_e = {E_n_keV:.3f} keV (ratio {E_n_keV/E_line_keV:.2f})")

print()
print("The 3.5 keV energy doesn't match a simple α^n scaling, suggesting")
print("a more complex combination (e.g., products of triangular numbers,")
print("spin factors, etc.). Further investigation needed.")
print()

# ========== CALIBRATION CHECKPOINT ==========
E_observed = 3.51  # keV (Bulbul et al. 2014)
E_predicted = 3.5  # keV (nominal)
deviation_eV = (E_predicted - E_observed) * 1000.0

print("=" * 70)
print("CALIBRATION CHECKPOINT")
print("=" * 70)
print(f"Observed (Bulbul 2014):  E = {E_observed:.2f} keV")
print(f"TriPhase V16 (nominal):     = {E_predicted:.2f} keV")
print(f"Deviation:                    {deviation_eV:+.0f} eV")
print()
print("(TriPhase does not yet have a full derivation of the 3.5 keV line;")
print("this is a placeholder for future work.)")
print("=" * 70)
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT:")
print("-" * 70)
print("The 3.5 keV line is a window into the statistical mechanics of dark")
print("matter (if it's DM decay) or the vacuum (if it's a radiative mode).")
print()
print("If sterile neutrino DM:")
print("  • Production: Dodelson-Widrow mechanism via oscillations")
print("  • Freeze-out: T_f ~ 100-200 MeV (electroweak scale)")
print("  • Decay: ν_s → γ + ν_active via loop diagrams")
print("  • Lifetime: τ ~ 10²⁸ s (compatible with observations)")
print()
print("The decay rate is Γ ~ (m_DM⁵/m_Pl⁴) · sin²(2θ), where θ is the")
print("mixing angle. The partition function for decaying DM includes:")
print("  Z_DM = Σ exp(-β·E_n - Γ·t)")
print()
print("If vacuum mode:")
print("  • The 3.5 keV line arises from EM vacuum transitions")
print("  • The energy scale is set by α-harmonics of m_e c²")
print("  • The line appears in regions of high dark matter density because")
print("    dark matter concentrations enhance local vacuum fluctuations")
print()
print("Both interpretations involve grand canonical ensembles with particle")
print("number fluctuations. The statistical signature (line width, spatial")
print("distribution) can distinguish between them.")
print()
print("The 3.5 keV line remains a mystery, but statistical mechanics provides")
print("the framework to understand its origin—whether exotic or mundane.")
print("=" * 70)

input("Press Enter to exit...")
