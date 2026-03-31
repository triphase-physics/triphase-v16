"""
================================================================================
TriPhase V16: lyman_alpha — Information Theory Framework
================================================================================

INFORMATION INTERPRETATION:
The Lyman-α spectral line (121.6 nm, 10.2 eV) encodes information about
the hydrogen energy level structure and serves as a cosmological information
carrier through the Lyman-α forest.

1. Spectral Information Content:
   - Each photon wavelength λ carries log₂(λ_max/Δλ) bits of frequency info
   - Spectroscopic resolution Δλ sets information capacity
   - Lyman-α forest: ~10⁴ absorption features → ~10⁴ bits of cosmological data

2. Channel Capacity:
   - Intergalactic medium (IGM) acts as a filter
   - Lyman-α forest encodes neutral hydrogen column density
   - Shannon capacity limited by line broadening and noise

3. Fisher Information:
   - Precision of redshift measurements from Lyman-α
   - Baryon Acoustic Oscillations (BAO) in Lyman-α forest
   - Constrains dark energy equation of state

4. Mutual Information:
   - I(Quasar spectrum ; IGM properties)
   - Correlation between different absorption systems
   - Cosmic web structure encoded in Lyman-α tomography

TRIPHASE DERIVATION:
Lyman-α transition: 2p → 1s in hydrogen
Energy: E = (3/4) × 13.6 eV = 10.2 eV
Wavelength: λ = hc/E = 121.567 nm

This is derived from Rydberg formula with TriPhase constants.

MIS TAG: (D) — Direct from atomic hydrogen

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
print("TriPhase V16: Lyman-Alpha Transition")
print("Information Theory Framework")
print("=" * 80)
print()

# ============================================================================
# STEP 1: Lyman-Alpha Energy and Wavelength
# ============================================================================
print("-" * 80)
print("STEP 1: Lyman-α Spectral Line Parameters")
print("-" * 80)
print()

# Rydberg constant from TriPhase
m_r = m_e * m_p / (m_e + m_p)  # Reduced mass
Ry_inf = m_e * c * alpha**2 / 2.0  # Rydberg energy (infinity)
Ry_H = m_r * c * alpha**2 / 2.0    # Rydberg energy (hydrogen)

# Lyman-alpha: n=2 to n=1
E_Lya = Ry_H * (1.0/1**2 - 1.0/2**2)
E_Lya *= 2.0  # Factor of 2 from Rydberg formula
E_Lya_eV = E_Lya / e
lambda_Lya = h * c / E_Lya
freq_Lya = c / lambda_Lya

print(f"Rydberg energy (hydrogen): Ry = {Ry_H/e:.6f} eV")
print()
print(f"Lyman-α transition (2p → 1s):")
print(f"  Energy: E = {E_Lya_eV:.6f} eV")
print(f"  Wavelength: λ = {lambda_Lya*1e9:.3f} nm")
print(f"  Frequency: f = {freq_Lya:.6e} Hz")
print()

# ============================================================================
# STEP 2: Spectral Resolution and Information Content
# ============================================================================
print("-" * 80)
print("STEP 2: Spectroscopic Information Capacity")
print("-" * 80)
print()

print("Information in wavelength measurement:")
print("  I = log₂(λ_range / Δλ)")
print()

# Typical spectroscopic resolution
R_spec = 10000  # Resolving power λ/Δλ
Delta_lambda = lambda_Lya / R_spec

print(f"Resolving power: R = λ/Δλ = {R_spec}")
print(f"Wavelength resolution: Δλ = {Delta_lambda*1e12:.3f} pm")
print()

# Information content
lambda_range = 300e-9  # UV-visible range
I_wavelength = math.log2(lambda_range / Delta_lambda)

print(f"Wavelength range: {lambda_range*1e9:.0f} nm")
print(f"Information content: I = log₂({lambda_range*1e9:.0f}nm / {Delta_lambda*1e12:.1f}pm)")
print(f"                      = {I_wavelength:.1f} bits")
print()
print("High-resolution spectroscopy can encode ~20-30 bits per line.")
print()

# ============================================================================
# STEP 3: Lyman-Alpha Forest Information
# ============================================================================
print("-" * 80)
print("STEP 3: Lyman-α Forest as Cosmological Barcode")
print("-" * 80)
print()

print("Quasar spectra show Lyman-α forest:")
print("  - Thousands of absorption lines")
print("  - Each line: redshift z, column density N_HI, line width b")
print()

# Typical forest statistics
N_lines = 5000  # absorption features
bits_per_line = 3  # Rough estimate: z (~10 bits), N_HI (~5 bits), b (~3 bits)
                   # But highly correlated, so effective ~3 bits

I_forest_total = N_lines * bits_per_line

print(f"Number of absorption lines: ~{N_lines}")
print(f"Effective bits per line: ~{bits_per_line}")
print(f"Total information: I ≈ {I_forest_total:.0f} bits")
print()
print("The Lyman-α forest encodes:")
print("  - IGM temperature and density")
print("  - Dark matter distribution (cosmic web)")
print("  - Cosmological parameters (H₀, Ω_m, etc.)")
print()

# ============================================================================
# STEP 4: Channel Capacity of IGM
# ============================================================================
print("-" * 80)
print("STEP 4: IGM as Information Channel")
print("-" * 80)
print()

print("Intergalactic medium acts as a filter:")
print("  - Input: Quasar continuum spectrum")
print("  - Filter: Lyman-α absorption (neutral hydrogen)")
print("  - Output: Observed spectrum with forest")
print()
print("Channel capacity (Shannon):")
print("  C = B log₂(1 + SNR)")
print()

# Typical spectral SNR in quasar observations
SNR_qso = 20  # per resolution element
bandwidth_spec = freq_Lya * 0.1  # 10% of central frequency

C_IGM = bandwidth_spec * math.log2(1.0 + SNR_qso)

print(f"Spectral SNR: {SNR_qso}")
print(f"Effective bandwidth: {bandwidth_spec:.6e} Hz")
print(f"Channel capacity: C ≈ {C_IGM:.6e} bits/s")
print()
print("(In practice, integration time limits total information)")
print()

# ============================================================================
# STEP 5: Fisher Information for BAO
# ============================================================================
print("-" * 80)
print("STEP 5: Fisher Information from Lyman-α BAO")
print("-" * 80)
print()

print("Baryon Acoustic Oscillations in Lyman-α forest:")
print("  - BAO scale: r_d ≈ 150 Mpc (comoving)")
print("  - Appears as correlation peak in Lyman-α clustering")
print()
print("Fisher information for cosmological parameters:")
print("  F(H₀, Ω_m, ...) from 3D correlation function")
print()

# Simplified estimate: Fisher info scales with number of modes
N_modes_BAO = 1000  # Effective number of independent modes
F_BAO_normalized = math.sqrt(N_modes_BAO)

print(f"Independent Fourier modes: N ≈ {N_modes_BAO}")
print(f"Fisher information: F ∝ √N ∝ {F_BAO_normalized:.1f}")
print()
print("Lyman-α BAO provides percent-level constraints on H₀ and Ω_m")
print("Comparable to CMB and galaxy surveys")
print()

# ============================================================================
# STEP 6: Mutual Information with Cosmic Web
# ============================================================================
print("-" * 80)
print("STEP 6: Mutual Information I(Lyman-α ; Cosmic Web)")
print("-" * 80)
print()

print("Lyman-α absorption traces large-scale structure:")
print("  - Overdense regions: strong absorption")
print("  - Voids: weak absorption")
print()
print("Mutual information quantifies correlation:")
print("  I(Lyα ; δ_matter) = H(Lyα) - H(Lyα | δ_matter)")
print()

# Estimate: strong correlation but not perfect
H_Lya = I_wavelength  # bits (from step 2)
H_Lya_given_matter = H_Lya * 0.3  # 70% reduction in uncertainty

I_mutual_cosmic = H_Lya - H_Lya_given_matter

print(f"Entropy H(Lyα): {H_Lya:.1f} bits")
print(f"Conditional entropy H(Lyα | matter): {H_Lya_given_matter:.1f} bits")
print(f"Mutual information: I = {I_mutual_cosmic:.1f} bits")
print()
print("Knowing dark matter distribution reduces Lyα uncertainty by ~70%")
print()

# ============================================================================
# STEP 7: Redshift as Information Carrier
# ============================================================================
print("-" * 80)
print("STEP 7: Redshift Encoding in Lyman-α")
print("-" * 80)
print()

print("Observed wavelength: λ_obs = λ_emit × (1 + z)")
print()
print("Each absorption system has redshift z:")
print()

# Example redshifts
z_values = [2.0, 2.5, 3.0, 3.5, 4.0]

print("z      λ_obs (nm)   Information content")
print("-" * 50)

for z in z_values:
    lambda_obs = lambda_Lya * (1.0 + z)
    # Information to specify z with precision Δz ~ 0.001
    Delta_z = 0.001
    I_redshift = math.log2((z + 1) / Delta_z) if z > 0 else 0
    print(f"{z:.1f}    {lambda_obs*1e9:6.1f}       {I_redshift:.1f} bits")

print()
print("Higher redshift → more information needed to specify z precisely")
print()

# ============================================================================
# STEP 8: Damping Wing and Information Loss
# ============================================================================
print("-" * 80)
print("STEP 8: Line Broadening and Information Degradation")
print("-" * 80)
print()

print("Absorption lines are broadened by:")
print("  - Thermal motion (Doppler)")
print("  - Hubble flow")
print("  - Natural line width (Lorentzian)")
print()

# Thermal broadening
T_IGM = 10000  # K (typical IGM temperature)
k_B = 1.380649e-23
m_H = m_p  # Hydrogen mass
v_thermal = math.sqrt(2.0 * k_B * T_IGM / m_H)
Delta_lambda_thermal = lambda_Lya * v_thermal / c

print(f"IGM temperature: T = {T_IGM:.0f} K")
print(f"Thermal velocity: v_th = {v_thermal/1000:.1f} km/s")
print(f"Thermal broadening: Δλ = {Delta_lambda_thermal*1e12:.3f} pm")
print()

# Information loss due to broadening
I_loss = math.log2(Delta_lambda_thermal / (lambda_Lya / R_spec))

print(f"Information loss: ΔI ≈ {I_loss:.1f} bits")
print()
print("Broadening degrades information — limits precision of z and N_HI")
print()

# ============================================================================
# STEP 9: Kolmogorov Complexity
# ============================================================================
print("-" * 80)
print("STEP 9: Algorithmic Complexity of Lyman-α")
print("-" * 80)
print()

print("Rydberg formula: 1/λ = R_∞ (1/n₁² - 1/n₂²)")
print("For Lyman-α: n₁ = 1, n₂ = 2")
print()
print("Complexity analysis:")
print("  - Rydberg constant R_∞ (from α, m_e, c, h)")
print("  - Two quantum numbers: 1, 2")
print("  - Formula: 1/λ = R(1 - 1/4) = (3/4)R")
print()

K_estimate = math.log2(4) + math.log2(3) + 10
print(f"Estimated Kolmogorov complexity: K(Lyα) ≈ {K_estimate:.1f} bits")
print()
print("Low complexity — Lyman-α is algorithmically simple")
print()

# ============================================================================
# STEP 10: Calibration Checkpoint
# ============================================================================
print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
print()

lambda_Lya_CODATA = 121.567e-9  # m (NIST/CODATA)
E_Lya_CODATA_eV = h * c / (lambda_Lya_CODATA * e)

deviation_nm = abs(lambda_Lya - lambda_Lya_CODATA) * 1e9
deviation_ppm = abs(lambda_Lya - lambda_Lya_CODATA) / lambda_Lya_CODATA * 1e6

print(f"TriPhase Lyman-α:")
print(f"  Wavelength: {lambda_Lya*1e9:.3f} nm")
print(f"  Energy: {E_Lya_eV:.6f} eV")
print()
print(f"CODATA Lyman-α:")
print(f"  Wavelength: {lambda_Lya_CODATA*1e9:.3f} nm")
print(f"  Energy: {E_Lya_CODATA_eV:.6f} eV")
print()
print(f"Deviation: {deviation_nm:.6f} nm ({deviation_ppm:.0f} ppm)")
print()

if deviation_ppm < 500:
    print("STATUS: EXCELLENT — TriPhase matches CODATA within 500 ppm")
elif deviation_ppm < 2000:
    print("STATUS: GOOD — TriPhase within 0.2% of CODATA")
else:
    print("STATUS: REVIEW — Deviation exceeds expected tolerance")

print()
print("=" * 80)
print("Information Theory Summary:")
print("=" * 80)
print(f"Lyman-α wavelength:                     {lambda_Lya*1e9:.3f} nm")
print(f"Lyman-α energy:                         {E_Lya_eV:.6f} eV")
print(f"Spectroscopic information (R=10⁴):      {I_wavelength:.1f} bits")
print(f"Lyman-α forest total information:       ~{I_forest_total:.0f} bits")
print(f"IGM channel capacity:                   {C_IGM:.6e} bits/s")
print(f"BAO Fisher information:                 ∝ {F_BAO_normalized:.1f}")
print(f"Mutual info (Lyα ; cosmic web):         {I_mutual_cosmic:.1f} bits")
print(f"Thermal broadening:                     {Delta_lambda_thermal*1e12:.3f} pm")
print(f"Information loss (broadening):          {I_loss:.1f} bits")
print(f"Kolmogorov complexity K(Lyα):           ~{K_estimate:.1f} bits")
print("=" * 80)
print()

input("Press Enter to exit...")
