"""
TriPhase V16 — Lyman Alpha Transition (Statistical Mechanics Framework)
========================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
The Lyman alpha transition (n=2 → n=1 in hydrogen) emerges from the canonical
ensemble of atomic energy levels. In the partition function Z = Σ_n exp(-βE_n),
the Boltzmann factor assigns statistical weight to each quantum state. The
transition energy ΔE = E_2 - E_1 = (3/4)·α²·m_e c² determines the photon frequency
ω = ΔE/ℏ emitted during de-excitation.

From the statistical mechanics perspective, the Lyman alpha line represents the
most probable transition in the hydrogen spectrum at moderate temperatures. The
Rydberg energy E_∞ = α²·m_e c²/2 sets the ionization threshold, and the n=2 level
sits at E_2 = -E_∞/4. The transition probability is proportional to the matrix
element ⟨1|r|2⟩ squared, which can be computed from the wave function overlap.

The Lyman alpha forest in cosmological spectra arises from absorption by neutral
hydrogen in the intergalactic medium. The forest's statistical properties (line
spacing, equivalent width distribution) encode the density of states of the cosmic
gas in the grand canonical ensemble at redshift z ~ 2-6.

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
print("TriPhase V16: Lyman Alpha Transition (Statistical Mechanics)")
print("=" * 70)
print()

print("HYDROGEN ENERGY LEVELS:")
print("-" * 70)
print("The Bohr model gives energy levels:")
print("  E_n = -E_∞/n², where E_∞ = α²·m_e c²/2 (Rydberg energy)")
print()

E_Rydberg = 0.5 * alpha**2 * m_e * c**2  # ionization energy
E_1 = -E_Rydberg  # ground state (n=1)
E_2 = -E_Rydberg / 4.0  # first excited state (n=2)
Delta_E = E_2 - E_1  # transition energy

print(f"Rydberg energy:      E_∞ = {E_Rydberg:.6e} J")
print(f"                         = {E_Rydberg / e:.6f} eV")
print()
print(f"Ground state (n=1):  E_1 = {E_1:.6e} J = {E_1/e:.6f} eV")
print(f"Excited state (n=2): E_2 = {E_2:.6e} J = {E_2/e:.6f} eV")
print()
print(f"Transition energy:   ΔE = E_2 - E_1 = {Delta_E:.6e} J")
print(f"                        = {Delta_E / e:.6f} eV")
print()

# Photon frequency and wavelength
nu_Lya = Delta_E / h  # frequency in Hz
lambda_Lya = c / nu_Lya  # wavelength in meters
lambda_Lya_nm = lambda_Lya * 1e9  # convert to nanometers

print(f"LYMAN ALPHA PHOTON:")
print("-" * 70)
print(f"Frequency:    ν = ΔE/h = {nu_Lya:.6e} Hz")
print(f"Wavelength:   λ = c/ν = {lambda_Lya:.6e} m")
print(f"                      = {lambda_Lya_nm:.4f} nm")
print()

print("STATISTICAL MECHANICS INTERPRETATION:")
print("-" * 70)
print("The partition function for hydrogen at temperature T is:")
print("  Z = Σ_n g_n exp(-βE_n), where β = 1/(k_B T)")
print()
print("The degeneracy is g_n = 2n² (including spin and orbital AM).")
print(f"  g_1 = 2·1² = 2")
print(f"  g_2 = 2·2² = 8")
print()

k_B = 1.380649e-23  # J/K
T_char = Delta_E / k_B  # characteristic temperature

print(f"For excitation to n=2, the characteristic temperature is:")
print(f"  T_char = ΔE/k_B = {T_char:.6e} K")
print()
print(f"At this temperature, exp(-βΔE) ~ exp(-1) ~ 0.37, so ~37% of atoms")
print(f"are excited to n=2 (ignoring degeneracy).")
print()

print("TRANSITION PROBABILITY:")
print("-" * 70)
print("The Einstein A coefficient for spontaneous emission is:")
print("  A_21 ∝ ω³ · |⟨1|r|2⟩|²")
print()
print("For Lyman alpha, A_21 ~ 6.3×10⁸ s⁻¹ (measured).")
print("The lifetime of the n=2 state is τ ~ 1/A_21 ~ 1.6 ns.")
print()
print("This rapid decay makes Lyman alpha the brightest line in the")
print("UV spectrum of neutral hydrogen.")
print()

# ========== CALIBRATION CHECKPOINT ==========
lambda_Lya_measured = 121.567  # nm (measured)
deviation_nm = lambda_Lya_nm - lambda_Lya_measured
deviation_ppm = (lambda_Lya_nm - lambda_Lya_measured) / lambda_Lya_measured * 1e6

print("=" * 70)
print("CALIBRATION CHECKPOINT")
print("=" * 70)
print(f"Measured Lyman α:       λ = {lambda_Lya_measured:.3f} nm")
print(f"TriPhase V16 (StatMech):  = {lambda_Lya_nm:.4f} nm")
print(f"Deviation:                  {deviation_nm:+.4f} nm ({deviation_ppm:+.0f} ppm)")
print("=" * 70)
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT:")
print("-" * 70)
print("The Lyman alpha line is a window into the statistical mechanics of")
print("the early universe. In the intergalactic medium at redshift z ~ 3,")
print("the gas temperature is T ~ 10⁴ K, which is comparable to T_char.")
print("This means neutral hydrogen is in thermal equilibrium between n=1")
print("and n=2 states.")
print()
print("The Lyman alpha forest—thousands of absorption lines in quasar spectra—")
print("traces the density fluctuations in the cosmic web. Each line represents")
print("a cloud of neutral hydrogen in a different region of the universe.")
print("The line spacing follows Poisson statistics, and the equivalent width")
print("distribution follows a power law—signatures of the grand canonical")
print("ensemble at work.")
print()
print("From the partition function perspective, Lyman alpha is special because")
print("it's the lowest-energy transition that can be accessed thermally. Higher")
print("transitions (Lyman beta, gamma, etc.) are exponentially suppressed by")
print("exp(-βΔE), so Lyman alpha dominates the spectrum.")
print()
print("This is the statistical origin of the Lyman alpha forest: it's the")
print("most probable signature of neutral hydrogen in the cosmic ensemble.")
print("=" * 70)

input("Press Enter to exit...")
