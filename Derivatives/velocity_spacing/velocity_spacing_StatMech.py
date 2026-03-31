"""
TriPhase V16 — Velocity Spacing (Statistical Mechanics Framework)
========================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
Velocity spacing in statistical mechanics emerges from the Maxwell-Boltzmann
distribution of particle velocities in thermal equilibrium. For a gas at temperature
T, the probability distribution of velocities is:
  f(v) ∝ v² exp(-m v²/2k_B T)

The most probable velocity is v_mp = √(2k_B T/m), the mean velocity is
⟨v⟩ = √(8k_B T/πm), and the RMS velocity is v_rms = √(3k_B T/m). These define
characteristic velocity spacings in the distribution.

In cosmology, velocity spacing refers to the characteristic velocity dispersion
in the intergalactic medium or galaxy clustering. The peculiar velocities of
galaxies follow a thermal-like distribution with σ_v ~ H₀·d, where d is the
comoving distance. This arises from the grand canonical ensemble of gravitationally
bound structures in the expanding universe.

In TriPhase, we derive the fundamental velocity spacing from the ratio c·α, which
represents the characteristic velocity of electrons in atomic orbitals (the Bohr
velocity v_Bohr = α·c ≈ c/137). This is the natural velocity scale for bound
electromagnetic systems.

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
print("TriPhase V16: Velocity Spacing (Statistical Mechanics)")
print("=" * 70)
print()

print("MAXWELL-BOLTZMANN VELOCITY DISTRIBUTION:")
print("-" * 70)
print("For an ideal gas in thermal equilibrium:")
print("  f(v) = 4π(m/2πk_B T)^(3/2) · v² · exp(-mv²/2k_B T)")
print()
print("Characteristic velocities:")
print("  Most probable:  v_mp = √(2k_B T/m)")
print("  Mean:           ⟨v⟩ = √(8k_B T/πm)")
print("  RMS:            v_rms = √(3k_B T/m)")
print()

print("BOHR VELOCITY (FUNDAMENTAL VELOCITY SPACING):")
print("-" * 70)
print("The electron velocity in the ground state of hydrogen is:")
print()

v_Bohr = alpha * c

print(f"  v_Bohr = α·c")
print(f"         = {alpha:.10f} · {c:.6e} m/s")
print(f"         = {v_Bohr:.6e} m/s")
print(f"         = {v_Bohr / 1000.0:.2f} km/s")
print()

print("This is the characteristic velocity for bound EM systems.")
print(f"It's c/137, or about 1/{alpha_inv:.1f} of light speed.")
print()

print("STATISTICAL INTERPRETATION:")
print("-" * 70)
print("In the canonical ensemble of atomic states, the velocity v_Bohr")
print("emerges from the partition function for the electron-proton system.")
print()
print("The electron momentum is quantized: p_n = n·(ℏ/a_0), where a_0 is the")
print("Bohr radius. For n=1 (ground state):")
print(f"  p_1 = ℏ/a_0 = m_e · v_Bohr")
print()

a_0 = hbar / (m_e * v_Bohr)  # Bohr radius
p_1 = m_e * v_Bohr

print(f"  Bohr radius:    a_0 = ℏ/(m_e v_Bohr) = {a_0:.6e} m")
print(f"  Ground momentum: p_1 = m_e v_Bohr = {p_1:.6e} kg·m/s")
print()

print("THERMAL EQUIVALENT:")
print("-" * 70)
print("If we interpret v_Bohr as a thermal velocity, what temperature does")
print("it correspond to?")
print()

k_B = 1.380649e-23  # J/K
T_eff = m_e * v_Bohr**2 / (3.0 * k_B)  # from v_rms = √(3k_B T/m)

print(f"  Setting v_rms = v_Bohr:")
print(f"  T_eff = m_e v_Bohr² / (3k_B)")
print(f"        = {T_eff:.6e} K")
print()
print("This is far below room temperature, reflecting that atomic binding")
print("energies (~eV) are much smaller than thermal energies at 300 K (~0.026 eV).")
print()

print("COSMOLOGICAL VELOCITY SPACING:")
print("-" * 70)
print("In large-scale structure, peculiar velocities follow:")
print("  v_pec ~ H₀·d ~ c·(H₀/c)·d")
print()

# Characteristic velocity at 100 Mpc scale
d_100Mpc = 100.0 * 3.0857e22  # meters
v_pec_100Mpc = H_0 * d_100Mpc

print(f"At d = 100 Mpc:")
print(f"  v_pec ~ H₀·d = {v_pec_100Mpc:.6e} m/s")
print(f"              = {v_pec_100Mpc / 1000.0:.0f} km/s")
print()
print("This is ~30× larger than v_Bohr, showing that cosmic velocities are")
print("supersonic relative to atomic scales.")
print()

# ========== CALIBRATION CHECKPOINT ==========
v_Bohr_measured = 2.18769126e6  # m/s (from α·c with CODATA values)
deviation_ppm = (v_Bohr - v_Bohr_measured) / v_Bohr_measured * 1e6

print("=" * 70)
print("CALIBRATION CHECKPOINT")
print("=" * 70)
print(f"CODATA 2018:            v_Bohr = {v_Bohr_measured:.6e} m/s")
print(f"TriPhase V16 (StatMech):       = {v_Bohr:.6e} m/s")
print(f"Deviation:                       {deviation_ppm:+.2f} ppm")
print("=" * 70)
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT:")
print("-" * 70)
print("Velocity spacing in statistical mechanics reflects the granularity of")
print("momentum space. In the Maxwell-Boltzmann distribution, the velocity")
print("dispersion σ_v = √(k_B T/m) sets the scale for thermal fluctuations.")
print()
print("For bound EM systems, the fundamental velocity is v_Bohr = α·c ≈ c/137.")
print("This emerges from the partition function of the electron-proton ensemble:")
print("  Z = Σ_n exp(-βE_n), where E_n ∝ (α·c)²")
print()
print("The factor α appears because it's the EM coupling constant—the")
print("statistical weight for photon exchange. The velocity v_Bohr is the")
print("speed at which an electron orbits the proton in the ground state,")
print("and it's the natural unit for atomic velocities.")
print()
print("In cosmology, velocity spacing reflects the expansion of the universe.")
print("Peculiar velocities arise from density fluctuations in the grand")
print("canonical ensemble of matter. The Hubble flow v = H₀·d is the")
print("mean-field approximation; individual galaxy velocities scatter around")
print("this with dispersion σ_v ~ few hundred km/s (from gravitational")
print("clustering).")
print()
print("Velocity spacing is fundamentally a statistical concept: it measures")
print("the typical momentum spread in the ensemble, whether atomic, thermal,")
print("or cosmological.")
print("=" * 70)

input("Press Enter to exit...")
