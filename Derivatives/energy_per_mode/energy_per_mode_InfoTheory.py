"""
================================================================================
TriPhase V16: energy_per_mode — Information Theory Framework
================================================================================

INFORMATION INTERPRETATION:
Energy per mode E = ℏω represents the fundamental information cost of
exciting a quantum harmonic oscillator mode. This is directly related to
Landauer's principle and quantum information storage.

1. Landauer Energy:
   - Erasing 1 bit requires E_min = k_B T ln(2)
   - At Compton temperature: E_bit ~ ℏω
   - Energy per mode sets minimum energy cost of information

2. Quantum Channel Capacity:
   - Number of photons in mode: n̄ = 1/(e^(ℏω/k_BT) - 1)
   - Entropy per mode: S = -Σ p_n ln(p_n)
   - ℏω determines thermal information content

3. Holevo Bound:
   - Quantum information capacity: C ≤ S(ρ)
   - For single mode: S = f(ℏω/k_BT)
   - Higher ℏω → fewer thermal excitations → less noise

4. Fisher Information:
   - Precision of frequency measurements
   - F(ω) ∝ 1/(Δω)² ~ t²
   - Energy resolution: ΔE ~ ℏ/Δt

MIS TAG: (D) — Direct from quantum mechanics

AUTHOR:  Christian R. Fuccillo
COMPANY: MIS Magnetic Innovative Solutions LLC
LICENSE: Proprietary
DOI:     10.5281/zenodo.17855383
DATE:    2025-2026

Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
All Rights Reserved.
================================================================================
"""

import math

# ============================================================================
# Anchor constants (TriPhase V16 Standard)
# ============================================================================
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

print("=" * 80)
print("TriPhase V16: Energy Per Mode (E = ℏω)")
print("Information Theory Framework")
print("=" * 80)
print()

k_B = 1.380649e-23  # J/K

# ============================================================================
# STEP 1: Landauer Energy at Different Scales
# ============================================================================
print("-" * 80)
print("STEP 1: Landauer Principle and Energy Per Mode")
print("-" * 80)
print()

print("Landauer principle: E_erase ≥ k_B T ln(2)")
print()

# Different temperature scales
T_room = 300  # K
T_LN2 = 77    # K (liquid nitrogen)
T_LHe = 4     # K (liquid helium)
T_dilution = 0.01  # K (dilution refrigerator)

temps = [T_room, T_LN2, T_LHe, T_dilution]
labels = ["Room temp", "LN₂ (77K)", "LHe (4K)", "Dilution (10mK)"]

for T, label in zip(temps, labels):
    E_Land = k_B * T * math.log(2)
    E_Land_eV = E_Land / e
    omega = E_Land / hbar
    freq = omega / (2.0 * math.pi)
    print(f"{label:20s}: E = {E_Land:.6e} J = {E_Land_eV:.6e} eV")
    print(f"{'':20s}  ω = {omega:.6e} rad/s, f = {freq:.6e} Hz")

print()
print("At room temperature, Landauer energy ~ 10⁻²¹ J ~ 0.02 eV")
print("Corresponds to microwave frequency ~ 6 THz")
print()

# ============================================================================
# STEP 2: Planck Distribution and Entropy
# ============================================================================
print("-" * 80)
print("STEP 2: Thermal Entropy Per Mode")
print("-" * 80)
print()

print("Bose-Einstein distribution: n̄ = 1 / (exp(ℏω/k_BT) - 1)")
print("Entropy per mode: S = -(1-p₀) ln(1-p₀) + ...")
print()

# Example: optical photon at room temperature
lambda_optical = 500e-9  # m (green light)
omega_optical = 2.0 * math.pi * c / lambda_optical
E_photon = hbar * omega_optical
E_photon_eV = E_photon / e

print(f"Optical photon (λ = {lambda_optical*1e9:.0f} nm):")
print(f"  Energy: E = ℏω = {E_photon:.6e} J = {E_photon_eV:.3f} eV")
print()

# Mean occupation at room temp
n_bar_optical = 1.0 / (math.exp(E_photon / (k_B * T_room)) - 1)

print(f"  Mean occupation at {T_room}K: n̄ = {n_bar_optical:.6e}")
print()
print("Optical photons have negligible thermal occupation at room temp.")
print("E >> k_B T → classical limit (low entropy)")
print()

# Microwave photon (thermal)
freq_MW = 10e9  # Hz (10 GHz)
omega_MW = 2.0 * math.pi * freq_MW
E_MW = hbar * omega_MW
n_bar_MW = 1.0 / (math.exp(E_MW / (k_B * T_room)) - 1)

print(f"Microwave photon (f = {freq_MW/1e9:.0f} GHz):")
print(f"  Energy: E = ℏω = {E_MW:.6e} J = {E_MW/e:.6e} eV")
print(f"  Mean occupation at {T_room}K: n̄ = {n_bar_MW:.3f}")
print()
print("Microwave modes have significant thermal occupation.")
print("E ~ k_B T → quantum effects important (high entropy)")
print()

# ============================================================================
# STEP 3: Shannon Information Capacity
# ============================================================================
print("-" * 80)
print("STEP 3: Photon Number States and Information")
print("-" * 80)
print()

print("Fock states |n⟩ with n photons:")
print("  log₂(N_max) bits needed to specify n ∈ {0, 1, ..., N_max}")
print()

# For thermal distribution, effective cutoff
N_max_optical = int(10 * n_bar_optical + 1)
N_max_MW = int(10 * n_bar_MW + 1)

info_optical = math.log2(max(N_max_optical, 2))
info_MW = math.log2(max(N_max_MW, 2))

print(f"Optical mode: Effective N_max ~ {N_max_optical}")
print(f"              Information: ~{info_optical:.3f} bits")
print()
print(f"Microwave mode: Effective N_max ~ {N_max_MW}")
print(f"                Information: ~{info_MW:.3f} bits")
print()
print("Higher temperature or lower frequency → more bits needed")
print()

# ============================================================================
# STEP 4: Holevo Bound for Single Mode
# ============================================================================
print("-" * 80)
print("STEP 4: Quantum Channel Capacity (Holevo Bound)")
print("-" * 80)
print()

print("Holevo bound: C ≤ S(ρ_avg) - Σ p_i S(ρ_i)")
print()
print("For coherent states |α⟩ in thermal noise:")
print("  C ≈ log₂(1 + n_signal / n_thermal)")
print()

n_signal = 100  # photons
n_thermal_optical = n_bar_optical
n_thermal_MW = n_bar_MW

SNR_optical = n_signal / max(n_thermal_optical, 1e-10)
SNR_MW = n_signal / max(n_thermal_MW, 1)

C_optical = math.log2(1.0 + SNR_optical)
C_MW = math.log2(1.0 + SNR_MW)

print(f"Signal: {n_signal} photons")
print()
print(f"Optical (low thermal noise):")
print(f"  SNR = {SNR_optical:.6e}")
print(f"  Capacity: C ≈ {C_optical:.3f} bits")
print()
print(f"Microwave (thermal noise):")
print(f"  SNR = {SNR_MW:.3f}")
print(f"  Capacity: C ≈ {C_MW:.3f} bits")
print()
print("Lower ℏω → higher thermal noise → lower channel capacity")
print()

# ============================================================================
# STEP 5: Energy-Time Uncertainty
# ============================================================================
print("-" * 80)
print("STEP 5: Energy Resolution and Fisher Information")
print("-" * 80)
print()

print("Energy-time uncertainty: ΔE Δt ≥ ℏ/2")
print()

# Measurement time determines energy resolution
delta_t = 1e-6  # s (1 microsecond)
Delta_E = hbar / (2.0 * delta_t)
Delta_E_eV = Delta_E / e

print(f"Measurement time: Δt = {delta_t*1e6:.1f} μs")
print(f"Energy resolution: ΔE ≥ ℏ/(2Δt) = {Delta_E:.6e} J")
print(f"                                 = {Delta_E_eV:.6e} eV")
print()

# Corresponding frequency uncertainty
Delta_omega = Delta_E / hbar
Delta_freq = Delta_omega / (2.0 * math.pi)

print(f"Frequency resolution: Δf = {Delta_freq:.6e} Hz")
print()

# Fisher information for frequency estimation
F_omega = delta_t**2
print(f"Fisher information: F(ω) ∝ t² ∝ {F_omega:.6e}")
print()
print("Longer measurement time → better frequency resolution → higher Fisher info")
print()

# ============================================================================
# STEP 6: Zero-Point Energy
# ============================================================================
print("-" * 80)
print("STEP 6: Zero-Point Energy and Vacuum Fluctuations")
print("-" * 80)
print()

print("Ground state energy: E₀ = ℏω/2")
print()
print("Even at T = 0, each mode has zero-point energy.")
print()

# Example: electron Compton frequency
omega_e = f_e * 2.0 * math.pi
E_ZP_e = 0.5 * hbar * omega_e
E_ZP_e_eV = E_ZP_e / e

print(f"Electron Compton mode:")
print(f"  Frequency: f_e = {f_e:.6e} Hz")
print(f"  Zero-point energy: E₀ = ℏω/2 = {E_ZP_e:.6e} J")
print(f"                                = {E_ZP_e_eV:.6e} eV")
print(f"                                = {E_ZP_e_eV/1e6:.3f} MeV")
print()
print("This is half the electron rest mass energy (0.511 MeV)")
print("Zero-point fluctuations carry information even in vacuum")
print()

# ============================================================================
# STEP 7: Photon Statistics and Entropy
# ============================================================================
print("-" * 80)
print("STEP 7: Sub-Poissonian vs Super-Poissonian Statistics")
print("-" * 80)
print()

print("Fock state |n⟩: ΔN = 0 (zero uncertainty)")
print("  Entropy: S = 0 (pure state)")
print()
print("Coherent state |α⟩: ΔN = √n̄ (Poissonian)")
print("  Entropy: S > 0 (mixed state in photon number basis)")
print()
print("Thermal state: ΔN > √n̄ (super-Poissonian)")
print("  Entropy: S(thermal) = (n̄+1)ln(n̄+1) - n̄ ln(n̄)")
print()

# Thermal entropy
if n_bar_MW > 1e-10:
    S_thermal_MW = (n_bar_MW + 1) * math.log(n_bar_MW + 1) - n_bar_MW * math.log(max(n_bar_MW, 1e-10))
    S_thermal_MW_bits = S_thermal_MW / math.log(2)
else:
    S_thermal_MW_bits = 0

print(f"Microwave thermal mode (n̄ = {n_bar_MW:.3f}):")
print(f"  Entropy: S = {S_thermal_MW_bits:.3f} bits")
print()

# ============================================================================
# STEP 8: Kolmogorov Complexity
# ============================================================================
print("-" * 80)
print("STEP 8: Algorithmic Complexity of E = ℏω")
print("-" * 80)
print()

print("Formula: E = ℏω")
print()
print("Complexity analysis:")
print("  - One fundamental constant: ℏ")
print("  - One parameter: ω (frequency)")
print("  - One operation: multiplication")
print()

K_estimate = math.log2(10) + 5
print(f"Estimated Kolmogorov complexity: K(E) ≈ {K_estimate:.1f} bits")
print()
print("Extremely simple formula — one of the most fundamental in physics")
print()

# ============================================================================
# STEP 9: Mutual Information in Coupled Oscillators
# ============================================================================
print("-" * 80)
print("STEP 9: Mutual Information Between Modes")
print("-" * 80)
print()

print("For two coupled oscillators with coupling g:")
print("  I(Mode₁ ; Mode₂) = f(g²/(ω₁ω₂))")
print()
print("Stronger coupling → higher mutual information")
print("  - Entangled modes share quantum information")
print("  - I quantifies correlation strength")
print()

# Example: weak coupling
g_coupling = 1e6  # Hz
omega_1 = omega_MW
omega_2 = omega_MW * 1.1

coupling_ratio = g_coupling**2 / (omega_1 * omega_2)

print(f"Coupling strength: g = {g_coupling/1e6:.1f} MHz")
print(f"Mode frequencies: ω₁ = {omega_1/(2*math.pi*1e9):.1f} GHz")
print(f"                  ω₂ = {omega_2/(2*math.pi*1e9):.1f} GHz")
print(f"Coupling ratio: g²/(ω₁ω₂) = {coupling_ratio:.6e}")
print()

# Estimate mutual info (weak coupling limit)
I_modes = math.log2(1.0 + coupling_ratio) if coupling_ratio > 0 else 0

print(f"Mutual information: I ≈ {I_modes:.6e} bits")
print()

# ============================================================================
# STEP 10: Calibration Checkpoint
# ============================================================================
print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
print()

# Test with known values
freq_test = 1e15  # Hz (optical)
omega_test = 2.0 * math.pi * freq_test
E_test = hbar * omega_test
E_test_eV = E_test / e

hbar_CODATA = 1.054571817e-34
E_test_CODATA = hbar_CODATA * omega_test
deviation_ppm = abs(E_test - E_test_CODATA) / E_test_CODATA * 1e6

print(f"Test frequency: f = {freq_test:.6e} Hz")
print(f"TriPhase E = ℏω:    {E_test:.6e} J = {E_test_eV:.6f} eV")
print(f"CODATA E = ℏω:      {E_test_CODATA:.6e} J")
print(f"Deviation:          {deviation_ppm:.3f} ppm")
print()

if deviation_ppm < 500:
    print("STATUS: EXCELLENT — TriPhase ℏ matches CODATA within 500 ppm")
else:
    print("STATUS: REVIEW — Check ℏ calibration")

print()
print("=" * 80)
print("Information Theory Summary:")
print("=" * 80)
print(f"Landauer energy (room temp):            {k_B*T_room*math.log(2):.6e} J")
print(f"Optical photon (500nm) energy:          {E_photon_eV:.3f} eV")
print(f"Thermal occupation (optical @ 300K):    {n_bar_optical:.6e}")
print(f"Thermal occupation (10 GHz @ 300K):     {n_bar_MW:.3f}")
print(f"Channel capacity (optical):             {C_optical:.3f} bits")
print(f"Channel capacity (microwave):           {C_MW:.3f} bits")
print(f"Energy resolution (1 μs):               {Delta_E_eV:.6e} eV")
print(f"Zero-point energy (Compton):            {E_ZP_e_eV/1e6:.3f} MeV")
print(f"Thermal entropy (MW mode):              {S_thermal_MW_bits:.3f} bits")
print(f"Kolmogorov complexity K(E=ℏω):          ~{K_estimate:.1f} bits")
print("=" * 80)
print()

input("Press Enter to exit...")
