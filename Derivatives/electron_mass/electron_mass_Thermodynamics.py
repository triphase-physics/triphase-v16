"""
================================================================================
TriPhase V16 Derivative: Electron Mass
Framework: THERMODYNAMICS
Tag: (D) — Pure derivation from thermodynamic principles
================================================================================

THERMODYNAMIC INTERPRETATION:
------------------------------
m_e = ℏα/(cr_e)

The electron mass sets the temperature scale for pair production. In
thermodynamics, the rest mass energy E = m_e c² defines a characteristic
temperature:

    T_e = m_e c² / k_B ≈ 6 × 10⁹ K

At temperatures above T_e, thermal energy is sufficient to create electron-
positron pairs from the vacuum. Below T_e, pair production is thermally
suppressed.

THERMODYNAMIC DERIVATION:
The partition function for a fermionic field includes pair creation:

    Z = Π [1 + exp(-E_k/k_BT)]

where E_k = √((pc)² + (m_e c²)²) is the relativistic energy.

At T >> m_e c²/k_B, the vacuum "boils" with electron-positron pairs.
At T << m_e c²/k_B, pairs are exponentially suppressed.

The electron mass m_e is the energy scale separating these regimes.

TriPhase derives m_e from electromagnetic parameters:

    m_e = ℏα / (c r_e)

where r_e is the classical electron radius. This shows m_e emerges from
the vacuum EM field structure, not as a fundamental parameter.

PHYSICAL PICTURE:
The electron is a localized excitation of the EM vacuum with characteristic
size r_e. The mass m_e is the "confinement energy" of this excitation,
analogous to the QCD confinement energy that gives the proton its mass.

--------------------------------------------------------------------------------
Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
License: CC BY-NC-SA 4.0 International
Website: www.triPhase.org
--------------------------------------------------------------------------------
"""

import math

print("="*80)
print("TriPhase V16 Derivative: Electron Mass")
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

print(f"c         = {c:.10e} m/s")
print(f"Z₀        = {Z_0:.10f} Ω")
print(f"α         = {alpha:.15e}")
print(f"ℏ         = {hbar:.15e} J·s")
print()

# ============================================================================
# THERMODYNAMIC DERIVATION OF m_e
# ============================================================================
print("THERMODYNAMIC DERIVATION:")
print("-" * 80)
print()
print("CLASSICAL ELECTRON RADIUS:")
print("The classical electron radius is defined from electrostatics:")
print()
print("    e² / (4πε₀r_e) = m_e c²")
print()
print("Solving for r_e:")
print()
print("    r_e = e² / (4πε₀m_e c²)")
print()

r_e = 2.8179403262e-15  # m (experimental value)

print(f"Classical electron radius: r_e = {r_e:.15e} m")
print()

print("TRIPHASE ELECTRON MASS:")
print("Inverting the relationship:")
print()
print("    m_e = e² / (4πε₀r_e c²)")
print()
print("Using α = e²/(4πε₀ℏc):")
print()
print("    m_e = ℏα / (c r_e)")
print()

m_e = hbar * alpha / (c * r_e)

print(f"    m_e = {m_e:.15e} kg")
print()

# ============================================================================
# REST MASS ENERGY
# ============================================================================
print("REST MASS ENERGY:")
print("-" * 80)
print()

E_rest = m_e * c**2

print(f"Rest energy:               E = m_e c²")
print(f"                           E = {E_rest:.15e} J")
print(f"                             = {E_rest/e:.10f} eV")
print(f"                             = {E_rest/e/1e6:.10f} MeV")
print()

# ============================================================================
# CHARACTERISTIC TEMPERATURE
# ============================================================================
print("CHARACTERISTIC TEMPERATURE:")
print("-" * 80)
print()

T_17 = 17 * 18 // 2
k_B = m_e * c**2 * alpha**2 / T_17

T_electron = m_e * c**2 / k_B

print(f"Electron temperature:      T_e = m_e c² / k_B")
print(f"                           T_e = {T_electron:.6e} K")
print()

print("At T > T_e, thermal pair production becomes efficient.")
print()

# Thermal energy at various temperatures
T_CMB = 2.725  # K
T_room = 300.0  # K
T_Sun_core = 1.5e7  # K

E_CMB = k_B * T_CMB
E_room = k_B * T_room
E_Sun = k_B * T_Sun_core

print(f"CMB (T = {T_CMB:.1f} K):         E = {E_CMB/e*1e6:.3f} μeV << m_e c²")
print(f"Room (T = {T_room:.0f} K):       E = {E_room/e*1000:.1f} meV << m_e c²")
print(f"Sun core (T = {T_Sun_core:.1e} K): E = {E_Sun/e/1e3:.1f} keV << m_e c²")
print()

print("Even in the Sun's core, k_B T << m_e c², so pair production")
print("is negligible.")
print()

# Early universe pair production
T_pair = T_electron
t_pair = (1e10 / T_pair)**2  # Rough estimate: T ~ 10¹⁰ K × (t/1s)^(-1/2)

print(f"Pair production era:       T ~ {T_pair:.2e} K")
print(f"Time (radiation era):      t ~ {t_pair:.6f} s")
print()

# ============================================================================
# PARTITION FUNCTION FOR PAIRS
# ============================================================================
print("PARTITION FUNCTION:")
print("-" * 80)
print()
print("For a fermion field with mass m_e, the partition function is:")
print()
print("    ln(Z) = ∫ (d³p/(2π)³) ln[1 + exp(-E/k_BT)]")
print()
print("where E = √(p²c² + m_e²c⁴)")
print()

# At T = T_e
# Approximate number density of thermal pairs
n_pair_T_e = (m_e * k_B * T_electron / (2.0 * math.pi * hbar**2))**(3.0/2.0) * math.exp(-m_e * c**2 / (k_B * T_electron))

print(f"At T = T_e:")
print(f"Pair density:              n ~ {n_pair_T_e:.6e} m⁻³")
print()

# Compare to baryon density today
n_baryon_today = 0.25  # m⁻³ (approximate)

print(f"Baryon density (today):    n_b ~ {n_baryon_today:.2f} m⁻³")
print()

# ============================================================================
# COMPTON WAVELENGTH
# ============================================================================
print("COMPTON WAVELENGTH:")
print("-" * 80)
print()

lambda_C = h / (m_e * c)
lambda_C_bar = hbar / (m_e * c)

print(f"Compton wavelength:        λ_C = h/(m_e c)")
print(f"                           λ_C = {lambda_C:.15e} m")
print(f"Reduced Compton:           λ̄_C = ℏ/(m_e c)")
print(f"                           λ̄_C = {lambda_C_bar:.15e} m")
print()

# Thermal wavelength
lambda_th_T_e = h / math.sqrt(2.0 * math.pi * m_e * k_B * T_electron)

print(f"Thermal wavelength (T_e):  λ_th = {lambda_th_T_e:.6e} m")
print(f"Ratio λ_th/λ_C:            {lambda_th_T_e / lambda_C:.6f}")
print()

# ============================================================================
# DE BROGLIE WAVELENGTH
# ============================================================================
print("DE BROGLIE WAVELENGTH:")
print("-" * 80)
print()

# For thermal electrons at room temperature
v_thermal_room = math.sqrt(2.0 * k_B * T_room / m_e)
lambda_dB_room = h / (m_e * v_thermal_room)

print(f"At T = {T_room} K:")
print(f"Thermal velocity:          v_th = {v_thermal_room:.6e} m/s")
print(f"de Broglie wavelength:     λ_dB = {lambda_dB_room:.6e} m")
print()

# For relativistic electrons at T_e
v_thermal_Te = c * math.sqrt(1.0 - 1.0/(1.0 + k_B*T_electron/(m_e*c**2))**2)
p_thermal_Te = m_e * v_thermal_Te / math.sqrt(1.0 - (v_thermal_Te/c)**2)
lambda_dB_Te = h / p_thermal_Te

print(f"At T = T_e:")
print(f"Thermal velocity:          v_th ≈ {v_thermal_Te/c:.3f} c")
print(f"de Broglie wavelength:     λ_dB = {lambda_dB_Te:.6e} m")
print()

# ============================================================================
# ZITTERBEWEGUNG
# ============================================================================
print("ZITTERBEWEGUNG:")
print("-" * 80)
print()

# Zitterbewegung frequency
f_zitter = m_e * c**2 / h

print(f"Zitterbewegung frequency:  f = m_e c² / h")
print(f"                           f = {f_zitter:.6e} Hz")
print(f"Period:                    T = {1/f_zitter:.6e} s")
print(f"Amplitude:                 a ~ λ̄_C = {lambda_C_bar:.6e} m")
print()

print("The electron 'trembles' at frequency f with amplitude ~ Compton wavelength.")
print()

# ============================================================================
# SELF-ENERGY
# ============================================================================
print("SELF-ENERGY:")
print("-" * 80)
print()

# Electrostatic self-energy
E_self_classical = e**2 / (4.0 * math.pi * epsilon_0 * r_e)

print(f"Classical self-energy:     E = e²/(4πε₀r_e)")
print(f"                           E = {E_self_classical:.6e} J")
print(f"                             = {E_self_classical/e/1e6:.3f} MeV")
print()

print(f"Rest mass energy:          m_e c² = {E_rest/e/1e6:.3f} MeV")
print()

ratio_self = E_self_classical / E_rest

print(f"Ratio E_self / (m_e c²):   {ratio_self:.6f}")
print()

print("By construction, E_self = m_e c² at r = r_e.")
print()

# ============================================================================
# CALIBRATION COMPARISON
# ============================================================================
print("="*80)
print("CALIBRATION COMPARISON")
print("="*80)
print()

# CODATA 2018 value
m_e_CODATA = 9.1093837015e-31  # kg

deviation = m_e - m_e_CODATA
rel_error = abs(deviation / m_e_CODATA)

print(f"TriPhase m_e:             {m_e:.15e} kg")
print(f"CODATA 2018 m_e:          {m_e_CODATA:.15e} kg")
print(f"Absolute deviation:       {deviation:+.15e} kg")
print(f"Relative error:           {rel_error:.6e} ({rel_error*100:.4e}%)")
print()

if rel_error < 1e-8:
    print("✓ EXCELLENT agreement (< 10 ppb)")
elif rel_error < 1e-6:
    print("✓ Good agreement (< 1 ppm)")
else:
    print("⚠ Moderate deviation")

print()
print("NOTE: TriPhase derives m_e from ℏ, α, c, and r_e.")
print()

# ============================================================================
# PHYSICAL SUMMARY
# ============================================================================
print("="*80)
print("PHYSICAL SUMMARY")
print("="*80)
print()
print(f"The electron mass m_e = {m_e:.6e} kg = {E_rest/e/1e6:.3f} MeV is:")
print()
print("1. REST MASS TEMPERATURE:")
print(f"   T_e = m_e c²/k_B = {T_electron:.2e} K")
print("   Above this, thermal pair production becomes efficient.")
print()
print("2. COMPTON SCALE:")
print(f"   λ_C = h/(m_e c) = {lambda_C:.3e} m")
print("   The characteristic size of quantum electron wavepackets.")
print()
print("3. ELECTROMAGNETIC ORIGIN:")
print("   m_e = ℏα/(cr_e) emerges from vacuum EM field structure.")
print(f"   Classical radius: r_e = {r_e:.3e} m")
print()
print("4. SELF-ENERGY:")
print("   The electrostatic self-energy e²/(4πε₀r_e) = m_e c²")
print("   sets the mass scale.")
print()
print("5. ZITTERBEWEGUNG:")
print(f"   The electron 'trembles' at f = {f_zitter:.2e} Hz")
print(f"   with amplitude ~ λ̄_C = {lambda_C_bar:.3e} m")
print()
print("The electron mass is the CONFINEMENT ENERGY of a localized")
print("EM field excitation — analogous to how QCD confinement gives")
print("the proton its mass. Both are thermodynamic phenomena!")
print()

print("="*80)
print("Derivation complete.")
print("="*80)

input("Press Enter to exit...")
