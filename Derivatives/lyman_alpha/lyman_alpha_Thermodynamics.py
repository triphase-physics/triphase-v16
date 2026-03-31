"""
================================================================================
TriPhase V16 Derivative: Lyman-Alpha Wavelength
Framework: THERMODYNAMICS
Tag: (D) — Pure derivation from thermodynamic principles
================================================================================

THERMODYNAMIC INTERPRETATION:
------------------------------
λ_Lyα = 4/(3R_∞) ≈ 121.5 nm

The Lyman-alpha transition (n=2 → n=1 in hydrogen) is not just a quantum
mechanical transition — it's a thermal emission line at the characteristic
temperature of atomic hydrogen.

In thermodynamics, emission spectra arise from thermal excitation and
relaxation. The Lyman-alpha energy corresponds to a temperature:

    T_Lyα = E_Lyα / k_B = (3/4) × 13.6 eV / k_B ≈ 1.18 × 10⁵ K

This is the temperature at which thermal collisions can efficiently excite
the n=2 level.

THERMODYNAMIC DERIVATION:
The Rydberg constant is:

    R_∞ = α² m_e c / (2h)

The Lyman-alpha transition energy is:

    E_Lyα = (1/1² - 1/2²) × 13.6 eV = (3/4) × 13.6 eV = 10.2 eV

The wavelength is:

    λ_Lyα = h c / E_Lyα = 4 / (3 R_∞)

THERMAL INTERPRETATION:
At temperature T_Lyα, the Boltzmann distribution gives the population ratio:

    n₂/n₁ = (g₂/g₁) exp(-E_Lyα/k_BT)

where g₂ = 8, g₁ = 2 are statistical weights. At T = T_Lyα, significant
population exists in the n=2 state, making Lyman-alpha emission efficient.

This is why Lyman-alpha is prominent in astrophysical plasmas: it's the
thermal signature of hydrogen at ~10⁵ K.

--------------------------------------------------------------------------------
Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
License: CC BY-NC-SA 4.0 International
Website: www.triPhase.org
--------------------------------------------------------------------------------
"""

import math

print("="*80)
print("TriPhase V16 Derivative: Lyman-Alpha Wavelength")
print("Framework: THERMODYNAMICS")
print("Tag: (D) — Pure derivation")
print("="*80)
print()

# ============================================================================
# STANDARD ANCHOR CHAIN
# ============================================================================
print("Building anchor chain from TriPhase fundamentals...")
print()

epsilon_0 = 8.8541878128e-12   # F/m (exact SI)
mu_0      = 1.25663706212e-6   # H/m (exact SI)
e         = 1.602176634e-19    # C (exact SI)

c   = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0 = math.sqrt(mu_0 / epsilon_0)

alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha     = 1.0 / alpha_inv

hbar = Z_0 * e**2 / (4.0 * math.pi * alpha)
h    = 2.0 * math.pi * hbar

r_e = 2.8179403262e-15
m_e = hbar * alpha / (c * r_e)

T_17 = 17 * 18 // 2
k_B = m_e * c**2 * alpha**2 / T_17

print(f"c         = {c:.10e} m/s")
print(f"α         = {alpha:.15e}")
print(f"ℏ         = {hbar:.15e} J·s")
print(f"h         = {h:.15e} J·s")
print(f"m_e       = {m_e:.15e} kg")
print(f"k_B       = {k_B:.15e} J/K")
print()

# ============================================================================
# RYDBERG CONSTANT
# ============================================================================
print("RYDBERG CONSTANT:")
print("-" * 80)
print()
print("The Rydberg constant for infinite mass:")
print()
print("    R_∞ = α² m_e c / (2h)")
print()

R_inf = alpha**2 * m_e * c / (2.0 * h)

print(f"R_∞ = {R_inf:.15e} m⁻¹")
print()

# Rydberg energy
E_Rydberg = h * c * R_inf

print(f"Rydberg energy:            E_Ry = hc R_∞")
print(f"                           E_Ry = {E_Rydberg:.15e} J")
print(f"                                = {E_Rydberg/e:.10f} eV")
print()

# ============================================================================
# LYMAN-ALPHA TRANSITION
# ============================================================================
print("LYMAN-ALPHA TRANSITION:")
print("-" * 80)
print()
print("Hydrogen energy levels:")
print()
print("    E_n = -13.6 eV / n²")
print()
print("Lyman-alpha: n=2 → n=1")
print()
print("    ΔE = E_2 - E_1 = 13.6 eV (1/1² - 1/2²)")
print("       = 13.6 eV × (3/4)")
print("       = 10.2 eV")
print()

# Energy levels
E_1 = -E_Rydberg
E_2 = -E_Rydberg / 4.0

E_Lyalpha = E_2 - E_1

print(f"Ground state:              E_1 = {E_1/e:.6f} eV")
print(f"First excited state:       E_2 = {E_2/e:.6f} eV")
print(f"Lyman-alpha energy:        ΔE = {E_Lyalpha/e:.10f} eV")
print()

# Wavelength
lambda_Lyalpha = h * c / E_Lyalpha

print(f"Lyman-alpha wavelength:    λ_Lyα = hc/ΔE")
print(f"                           λ_Lyα = {lambda_Lyalpha:.15e} m")
print(f"                                 = {lambda_Lyalpha*1e9:.6f} nm")
print()

# Alternative formula
lambda_Lyalpha_alt = 4.0 / (3.0 * R_inf)

print(f"Alternative formula:       λ_Lyα = 4/(3R_∞)")
print(f"                           λ_Lyα = {lambda_Lyalpha_alt:.15e} m")
print(f"                                 = {lambda_Lyalpha_alt*1e9:.6f} nm")
print()

# Frequency
f_Lyalpha = c / lambda_Lyalpha

print(f"Frequency:                 f = c/λ = {f_Lyalpha:.6e} Hz")
print()

# ============================================================================
# THERMODYNAMIC INTERPRETATION
# ============================================================================
print("THERMODYNAMIC INTERPRETATION:")
print("-" * 80)
print()
print("The Lyman-alpha energy sets a characteristic temperature:")
print()

T_Lyalpha = E_Lyalpha / k_B

print(f"Lyman-alpha temperature:   T_Lyα = ΔE/k_B")
print(f"                           T_Lyα = {T_Lyalpha:.6e} K")
print()

print("This is the temperature at which thermal collisions can efficiently")
print("excite hydrogen to the n=2 state.")
print()

# ============================================================================
# BOLTZMANN POPULATION
# ============================================================================
print("BOLTZMANN POPULATION:")
print("-" * 80)
print()
print("The population ratio between levels n=2 and n=1 is:")
print()
print("    n₂/n₁ = (g₂/g₁) exp(-ΔE/k_BT)")
print()

# Statistical weights
g_1 = 2  # 2s+1 = 2 for n=1
g_2 = 8  # 2n² = 8 for n=2

print(f"Statistical weight n=1:    g₁ = {g_1}")
print(f"Statistical weight n=2:    g₂ = {g_2}")
print()

# Population ratios at various temperatures
T_test = [1e3, 1e4, T_Lyalpha, 1e6]

print("Population ratios at various temperatures:")
print()
for T in T_test:
    ratio = (g_2 / g_1) * math.exp(-E_Lyalpha / (k_B * T))
    print(f"  T = {T:.3e} K:  n₂/n₁ = {ratio:.6e}")

print()

# Temperature for 10% excitation
n2_n1_target = 0.1
T_10percent = E_Lyalpha / (k_B * math.log(g_2 / (g_1 * n2_n1_target)))

print(f"For 10% excitation:        T = {T_10percent:.6e} K")
print()

# ============================================================================
# THERMAL EMISSION
# ============================================================================
print("THERMAL EMISSION:")
print("-" * 80)
print()
print("At temperature T, the emission coefficient is:")
print()
print("    j_ν ∝ n₂ A_21 exp(-hν/k_BT)")
print()

# Einstein A coefficient (approximate)
A_21 = 6.27e8  # s⁻¹ for Lyman-alpha

print(f"Einstein A coefficient:    A_21 = {A_21:.3e} s⁻¹")
print(f"Spontaneous lifetime:      τ = 1/A_21 = {1/A_21:.3e} s")
print()

# Emission rate at T_Lyalpha
n_H = 1e20  # Typical stellar atmosphere density (m⁻³)
n_2 = n_H * (g_2 / g_1) * math.exp(-E_Lyalpha / (k_B * T_Lyalpha))
emission_rate = n_2 * A_21

print(f"Density (T={T_Lyalpha:.2e} K):  n_H = {n_H:.6e} m⁻³")
print(f"Excited state density:     n₂ = {n_2:.6e} m⁻³")
print(f"Emission rate:             j = n₂ A_21 = {emission_rate:.6e} s⁻¹ m⁻³")
print()

# ============================================================================
# ASTROPHYSICAL APPLICATIONS
# ============================================================================
print("ASTROPHYSICAL APPLICATIONS:")
print("-" * 80)
print()

# Solar chromosphere
T_solar_chrom = 2e4  # K
n_2_solar = n_H * (g_2 / g_1) * math.exp(-E_Lyalpha / (k_B * T_solar_chrom))
ratio_solar = n_2_solar / n_H

print(f"Solar chromosphere:")
print(f"  Temperature:             T = {T_solar_chrom:.2e} K")
print(f"  Excitation fraction:     n₂/n_H = {ratio_solar:.6e}")
print()

# Intergalactic medium (z~3)
T_IGM = 1e4  # K
n_2_IGM = n_H * (g_2 / g_1) * math.exp(-E_Lyalpha / (k_B * T_IGM))
ratio_IGM = n_2_IGM / n_H

print(f"Intergalactic medium:")
print(f"  Temperature:             T = {T_IGM:.2e} K")
print(f"  Excitation fraction:     n₂/n_H = {ratio_IGM:.6e}")
print()

# Lyman-alpha forest
print("The Lyman-alpha forest in quasar spectra traces neutral hydrogen")
print("in the IGM at temperatures T ~ 10⁴ K — the thermal signature!")
print()

# ============================================================================
# PARTITION FUNCTION
# ============================================================================
print("PARTITION FUNCTION:")
print("-" * 80)
print()
print("The partition function for hydrogen (first few levels):")
print()
print("    Z = Σ g_n exp(-E_n/k_BT)")
print()

# Calculate at T_Lyalpha
Z_levels = []
for n in range(1, 6):
    g_n = 2 * n**2
    E_n = -E_Rydberg / n**2
    term = g_n * math.exp(-E_n / (k_B * T_Lyalpha))
    Z_levels.append((n, g_n, E_n, term))

Z_total = sum([term for _, _, _, term in Z_levels])

print(f"At T = {T_Lyalpha:.3e} K:")
print()
for n, g_n, E_n, term in Z_levels:
    print(f"  n={n}: g={g_n:2d}, E={E_n/e:+.3f} eV, contrib={term:.6e}")

print()
print(f"Total partition function:  Z = {Z_total:.6e}")
print()

# Helmholtz free energy
F_hydrogen = -k_B * T_Lyalpha * math.log(Z_total)

print(f"Free energy:               F = -k_B T ln(Z)")
print(f"                           F = {F_hydrogen:.6e} J")
print(f"                                = {F_hydrogen/e:.3f} eV")
print()

# ============================================================================
# THERMAL BROADENING
# ============================================================================
print("THERMAL BROADENING:")
print("-" * 80)
print()
print("The Lyman-alpha line is broadened by thermal motion (Doppler):")
print()

m_H = 1.67262192369e-27  # Hydrogen atom mass (kg)

# Thermal velocity
v_thermal = math.sqrt(2.0 * k_B * T_Lyalpha / m_H)

print(f"Hydrogen mass:             m_H = {m_H:.6e} kg")
print(f"Thermal velocity:          v_th = √(2k_BT/m_H)")
print(f"                           v_th = {v_thermal:.6e} m/s")
print(f"                           v_th/c = {v_thermal/c:.6e}")
print()

# Doppler width
delta_lambda_thermal = lambda_Lyalpha * v_thermal / c

print(f"Doppler width:             Δλ = λ v_th/c")
print(f"                           Δλ = {delta_lambda_thermal:.6e} m")
print(f"                                = {delta_lambda_thermal*1e12:.3f} pm")
print()

# Frequency width
delta_f_thermal = f_Lyalpha * v_thermal / c

print(f"Frequency width:           Δf = f v_th/c")
print(f"                           Δf = {delta_f_thermal:.6e} Hz")
print()

# ============================================================================
# CALIBRATION COMPARISON
# ============================================================================
print("="*80)
print("CALIBRATION COMPARISON")
print("="*80)
print()

# NIST/CODATA value for Lyman-alpha
lambda_Lyalpha_NIST = 121.567e-9  # m

deviation = lambda_Lyalpha - lambda_Lyalpha_NIST
rel_error = abs(deviation / lambda_Lyalpha_NIST)

print(f"TriPhase λ_Lyα:           {lambda_Lyalpha:.15e} m")
print(f"                          {lambda_Lyalpha*1e9:.10f} nm")
print()
print(f"NIST λ_Lyα:               {lambda_Lyalpha_NIST:.15e} m")
print(f"                          {lambda_Lyalpha_NIST*1e9:.10f} nm")
print()
print(f"Absolute deviation:       {deviation:+.15e} m")
print(f"                          {deviation*1e12:+.6f} pm")
print(f"Relative error:           {rel_error:.6e} ({rel_error*100:.4e}%)")
print()

if rel_error < 1e-6:
    print("✓ EXCELLENT agreement (< 1 ppm)")
elif rel_error < 1e-4:
    print("✓ Good agreement (< 100 ppm)")
else:
    print("⚠ Moderate deviation")

print()
print("NOTE: TriPhase derives λ_Lyα from R_∞ = α²m_ec/(2h).")
print()

# ============================================================================
# PHYSICAL SUMMARY
# ============================================================================
print("="*80)
print("PHYSICAL SUMMARY")
print("="*80)
print()
print(f"The Lyman-alpha wavelength λ_Lyα = {lambda_Lyalpha*1e9:.3f} nm is:")
print()
print("1. THERMAL EMISSION LINE:")
print(f"   Temperature: T_Lyα = {T_Lyalpha:.2e} K")
print("   At this temperature, thermal collisions efficiently populate n=2.")
print()
print("2. ASTROPHYSICAL TRACER:")
print("   Prominent in stellar chromospheres, HII regions, and IGM.")
print("   Traces hydrogen at T ~ 10⁴-10⁵ K.")
print()
print("3. BOLTZMANN DISTRIBUTION:")
print("   Population ratio n₂/n₁ = (g₂/g₁)exp(-10.2 eV/k_BT)")
print("   Significant excitation requires T ≳ 10⁴ K.")
print()
print("4. THERMAL BROADENING:")
print(f"   Doppler width Δλ ~ {delta_lambda_thermal*1e12:.1f} pm")
print(f"   From thermal motion at v_th = {v_thermal/1000:.1f} km/s")
print()
print("5. COSMOLOGICAL PROBE:")
print("   Lyman-alpha forest maps neutral hydrogen in the early universe.")
print("   Each absorption line is a thermal signature at specific redshift.")
print()
print("Lyman-alpha is not just a spectral line — it's a thermometer")
print("revealing the temperature of cosmic hydrogen.")
print()

print("="*80)
print("Derivation complete.")
print("="*80)

input("Press Enter to exit...")
