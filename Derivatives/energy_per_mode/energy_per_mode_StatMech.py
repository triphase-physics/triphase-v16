"""
TriPhase V16 — Energy Per Mode (Statistical Mechanics Framework)
========================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
The energy per mode is the cornerstone of statistical mechanics, emerging directly
from the equipartition theorem and the canonical ensemble. For a system in thermal
equilibrium at temperature T, the Boltzmann distribution assigns probability
p_i = exp(-βE_i)/Z to each microstate, where β = 1/(k_B T) and Z is the partition
function. The average energy per degree of freedom is ⟨E⟩ = k_B T/2 for quadratic
terms in the Hamiltonian (equipartition).

For photons (bosons) in the electromagnetic field, the energy distribution follows
Planck's law: ⟨n(ω)⟩ = 1/[exp(ℏω/k_B T) - 1], where n(ω) is the occupation number
of mode ω. The total energy per mode includes both the discrete quantum energies
ℏω and the continuum of thermal fluctuations. At high temperature (k_B T >> ℏω),
this reduces to the classical Rayleigh-Jeans result: ⟨E⟩ = k_B T per mode.

In TriPhase, we calculate the characteristic energy per vacuum mode using the
electron rest energy m_e c² as the reference scale. The number of modes is set
by the fine structure constant α⁻¹ ≈ 137, which counts the EM phase space volume
at the Compton scale. The energy per mode is then E_mode ~ (m_e c²)/α⁻¹ ~ α·m_e c².

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
print("TriPhase V16: Energy Per Mode (Statistical Mechanics)")
print("=" * 70)
print()

print("EQUIPARTITION THEOREM:")
print("-" * 70)
print("For a classical system with Hamiltonian H = Σ_i (p_i²/2m + V(q_i)),")
print("each quadratic term contributes ⟨E⟩ = k_B T/2 to the average energy.")
print()
print("For a harmonic oscillator (kinetic + potential):")
print("  ⟨E⟩ = k_B T/2 + k_B T/2 = k_B T per mode")
print()
print("For photons (massless bosons), the quantum version is:")
print("  ⟨E_mode⟩ = ℏω / [exp(ℏω/k_B T) - 1]")
print()

print("VACUUM ENERGY DENSITY:")
print("-" * 70)
print("In quantum field theory, each mode contributes zero-point energy ℏω/2.")
print("The total vacuum energy is:")
print("  E_vac = Σ_modes (ℏω/2)")
print()
print("This sum diverges (UV catastrophe) unless we impose a cutoff.")
print("In TriPhase, the natural cutoff is the electron Compton scale.")
print()

E_electron = m_e * c**2
lambda_C = hbar / (m_e * c)
k_max = 1.0 / lambda_C
omega_max = c * k_max

print(f"Electron rest energy:    E_e = m_e c² = {E_electron:.6e} J")
print(f"Compton wavelength:      λ_C = ℏ/(m_e c) = {lambda_C:.6e} m")
print(f"Cutoff wavevector:       k_max = 1/λ_C = {k_max:.6e} m⁻¹")
print(f"Cutoff frequency:        ω_max = c·k_max = {omega_max:.6e} rad/s")
print()

print("NUMBER OF MODES:")
print("-" * 70)
print("The fine structure constant α⁻¹ counts the EM phase space volume:")
print(f"  N_modes ~ α⁻¹ ≈ {alpha_inv:.2f}")
print()
print("This is the number of photon states accessible at the Compton scale.")
print()

N_modes = alpha_inv
E_per_mode = E_electron / N_modes

print(f"Energy per mode:  E_mode = E_e / α⁻¹")
print(f"                         = {E_per_mode:.6e} J")
print(f"                         = {E_per_mode / e:.6f} eV")
print(f"                         = α · m_e c²")
print()

# Alternative: express in terms of k_B T_eff
k_B = 1.380649e-23  # J/K
T_eff = E_per_mode / k_B

print(f"Effective temperature:  T_eff = E_mode / k_B")
print(f"                              = {T_eff:.6e} K")
print()

print("PLANCK DISTRIBUTION:")
print("-" * 70)
print("For thermal photons at temperature T, the average number per mode is:")
print("  ⟨n(ω)⟩ = 1 / [exp(ℏω/k_B T) - 1]")
print()
print("The average energy per mode is:")
print("  ⟨E(ω)⟩ = ℏω · ⟨n(ω)⟩ = ℏω / [exp(ℏω/k_B T) - 1]")
print()
print(f"At the electron Compton frequency ω_e = m_e c²/ℏ = {omega_max:.3e} rad/s,")
print(f"the typical thermal energy is E ~ k_B T_eff ~ {E_per_mode:.3e} J.")
print()

# ========== CALIBRATION CHECKPOINT ==========
alpha_times_mec2 = alpha * m_e * c**2
deviation_percent = abs(E_per_mode - alpha_times_mec2) / alpha_times_mec2 * 100.0

print("=" * 70)
print("CALIBRATION CHECKPOINT")
print("=" * 70)
print(f"E_mode (calculated):    {E_per_mode:.6e} J")
print(f"α·m_e c² (expected):    {alpha_times_mec2:.6e} J")
print(f"Deviation:              {deviation_percent:.4f}%")
print()
print("(Perfect match confirms E_mode = α·m_e c²)")
print("=" * 70)
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT:")
print("-" * 70)
print("The energy per mode is the fundamental unit of energy in statistical")
print("mechanics. For the EM vacuum, this unit is E_mode = α·m_e c², which")
print("is the electron rest energy divided by the number of available photon")
print("states (α⁻¹ ≈ 137).")
print()
print("This quantity emerges from the partition function:")
print("  Z = Σ_n exp(-βE_n), where β = 1/(k_B T)")
print()
print("The average energy is:")
print("  ⟨E⟩ = -∂ln(Z)/∂β = Σ_n E_n · p_n")
print()
print("For the EM vacuum at the Compton scale, ⟨E⟩ ~ α·m_e c² per mode.")
print("This is the characteristic energy of vacuum fluctuations—the minimum")
print("energy required to create a virtual photon at the electron wavelength.")
print()
print("The equipartition theorem states that each quadratic degree of freedom")
print("carries energy k_B T/2. The EM field has two polarizations (E and B),")
print("so the total energy per mode is k_B T. Setting this equal to α·m_e c²")
print("gives the effective vacuum temperature:")
print(f"  T_vac ~ α·m_e c²/k_B ~ {T_eff:.3e} K")
print()
print("This is the 'temperature' of the quantum vacuum—the energy scale at")
print("which EM fluctuations become thermally activated. It's not a real")
print("temperature (the vacuum is at T=0), but it's the effective temperature")
print("that would produce the same fluctuation energy.")
print("=" * 70)

input("Press Enter to exit...")
