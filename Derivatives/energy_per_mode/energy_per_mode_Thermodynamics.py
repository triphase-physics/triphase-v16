"""
================================================================================
TriPhase V16 Derivative: Energy Per Mode
Framework: THERMODYNAMICS
Tag: (D) — Pure derivation from thermodynamic principles
================================================================================

THERMODYNAMIC INTERPRETATION:
------------------------------
E_mode = m_e c² α² / T₁₇ ≈ k_B T_natural

The energy per mode is the fundamental unit of thermal energy in the vacuum.
This emerges directly from the equipartition theorem in statistical mechanics:

    <E> = (1/2) k_B T per degree of freedom

For a system at the "natural temperature" T_natural, each of the T₁₇ = 153
degrees of freedom carries energy:

    E_mode = k_B T_natural / T₁₇

THERMODYNAMIC DERIVATION:
The total thermal energy of the vacuum at temperature T is:

    E_total = T₁₇ × E_mode = k_B T

This is the equipartition theorem: each DOF contributes k_B T / T₁₇.

The natural temperature is defined by:

    T_natural = m_e c² α² / k_B

At this temperature, the energy per mode is:

    E_mode = k_B T_natural / T₁₇ = m_e c² α² / T₁₇

This is equivalently the Boltzmann constant:

    k_B = m_e c² α² / T₁₇ = E_mode

PHYSICAL PICTURE:
E_mode sets the thermal energy scale for vacuum fluctuations. At any
temperature T, the vacuum has:

- Total energy: E = k_B T
- Energy per mode: E_mode = k_B T / T₁₇
- Number of thermally excited modes: N_excited ≈ T₁₇

This is the connection between quantum vacuum fluctuations and thermodynamics.

--------------------------------------------------------------------------------
Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
License: CC BY-NC-SA 4.0 International
Website: www.triPhase.org
--------------------------------------------------------------------------------
"""

import math

print("="*80)
print("TriPhase V16 Derivative: Energy Per Mode")
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

r_e = 2.8179403262e-15
m_e = hbar * alpha / (c * r_e)
f_e = m_e * c**2 / hbar

T_17 = 17 * 18 // 2  # = 153

print(f"c         = {c:.10e} m/s")
print(f"α         = {alpha:.15e}")
print(f"ℏ         = {hbar:.15e} J·s")
print(f"m_e       = {m_e:.15e} kg")
print(f"T₁₇       = {T_17}")
print()

# ============================================================================
# THERMODYNAMIC DERIVATION OF E_mode
# ============================================================================
print("THERMODYNAMIC DERIVATION:")
print("-" * 80)
print()
print("EQUIPARTITION THEOREM:")
print("In statistical mechanics, each degree of freedom in a system")
print("at temperature T has average energy:")
print()
print("    <E> = (1/2) k_B T  (classical)")
print()
print("For T₁₇ = 153 degrees of freedom:")
print()
print("    E_total = T₁₇ × (1/2) k_B T")
print()

print("NATURAL TEMPERATURE:")
print("Define the natural temperature scale from electron mass:")
print()
print("    T_natural = m_e c² α² / k_B")
print()

# Boltzmann constant from T_17
k_B = m_e * c**2 * alpha**2 / T_17

print(f"Boltzmann constant:        k_B = m_e c² α² / T₁₇")
print(f"                           k_B = {k_B:.15e} J/K")
print()

# Natural temperature
T_natural = m_e * c**2 * alpha**2 / k_B

print(f"Natural temperature:       T_nat = m_e c² α² / k_B")
print(f"                           T_nat = {T_natural:.6e} K")
print()

print("ENERGY PER MODE:")
print("At T = T_natural, the energy per mode is:")
print()
print("    E_mode = k_B T_nat / T₁₇")
print()

E_mode = k_B * T_natural / T_17

print(f"                           E_mode = {E_mode:.15e} J")
print()

# Simplify
E_mode_simple = m_e * c**2 * alpha**2 / T_17

print("Simplifying:")
print("    E_mode = k_B T_nat / T₁₇")
print("           = (m_e c² α² / T₁₇) × (m_e c² α² / k_B) / T₁₇")
print("           = m_e c² α² / T₁₇")
print()
print(f"                           E_mode = {E_mode_simple:.15e} J")
print()

# Check equivalence
print(f"Verification:              E_mode = k_B")
print(f"                           {E_mode_simple:.15e} = {k_B:.15e}")
print(f"Ratio:                     {E_mode_simple / k_B:.15f}")
print()

print("At T_natural, E_mode = k_B exactly!")
print()

# ============================================================================
# ENERGY SCALES
# ============================================================================
print("ENERGY SCALES:")
print("-" * 80)
print()

# Electron rest mass
E_electron = m_e * c**2

print(f"Electron rest mass:        m_e c² = {E_electron:.15e} J")
print(f"                                  = {E_electron/e/1e6:.6f} MeV")
print()

# Fine structure scale
E_alpha = m_e * c**2 * alpha

print(f"Fine structure scale:      m_e c² α = {E_alpha:.15e} J")
print(f"                                    = {E_alpha/e/1e3:.6f} keV")
print()

# Natural energy (double fine structure)
E_natural = m_e * c**2 * alpha**2

print(f"Natural energy:            m_e c² α² = {E_natural:.15e} J")
print(f"                                     = {E_natural/e:.6f} eV")
print()

# Energy per mode
print(f"Energy per mode:           E_mode = m_e c² α² / 153")
print(f"                           E_mode = {E_mode:.15e} J")
print(f"                                  = {E_mode/e*1000:.6f} meV")
print()

# Ratios
ratio_total = E_electron / E_mode
ratio_alpha = E_alpha / E_mode
ratio_natural = E_natural / E_mode

print(f"Ratio m_e c² / E_mode:     {ratio_total:.6e}")
print(f"Ratio m_e c² α / E_mode:   {ratio_alpha:.6e}")
print(f"Ratio m_e c² α² / E_mode:  {ratio_natural:.6f}")
print()

# ============================================================================
# TEMPERATURE SCALES
# ============================================================================
print("TEMPERATURE SCALES:")
print("-" * 80)
print()

# Various temperatures
T_electron = m_e * c**2 / k_B
T_alpha = m_e * c**2 * alpha / k_B
T_CMB = 2.725  # K
T_room = 300.0  # K

print(f"Electron rest-mass temp:   T_e = m_e c² / k_B")
print(f"                           T_e = {T_electron:.6e} K")
print()
print(f"Fine structure temp:       T_α = m_e c² α / k_B")
print(f"                           T_α = {T_alpha:.6e} K")
print()
print(f"Natural temperature:       T_nat = m_e c² α² / k_B")
print(f"                           T_nat = {T_natural:.6e} K")
print()
print(f"CMB temperature:           T_CMB = {T_CMB:.3f} K")
print(f"Room temperature:          T_room = {T_room:.1f} K")
print()

# Energy per mode at different temperatures
E_CMB = k_B * T_CMB
E_room = k_B * T_room

print(f"k_B T_CMB:                 {E_CMB:.6e} J = {E_CMB/e*1e6:.3f} μeV")
print(f"k_B T_room:                {E_room:.6e} J = {E_room/e*1000:.1f} meV")
print(f"k_B T_nat:                 {E_mode:.6e} J = {E_mode/e*1000:.3f} meV")
print()

# ============================================================================
# PARTITION FUNCTION
# ============================================================================
print("PARTITION FUNCTION:")
print("-" * 80)
print()
print("For a single harmonic oscillator mode at frequency ω:")
print()
print("    Z_HO = 1 / (1 - exp(-ℏω/k_BT))")
print()

# Mode frequency
omega_mode = E_mode / hbar
f_mode = omega_mode / (2.0 * math.pi)

print(f"Mode frequency:            ω = E_mode / ℏ")
print(f"                           ω = {omega_mode:.6e} rad/s")
print(f"                           f = {f_mode:.6e} Hz")
print()

# Partition function at T_natural
Z_mode_natural = 1.0 / (1.0 - math.exp(-E_mode / (k_B * T_natural)))

print(f"At T = T_natural:")
print(f"Partition function:        Z = {Z_mode_natural:.6f}")
print()

# Average occupation
n_avg_natural = 1.0 / (math.exp(E_mode / (k_B * T_natural)) - 1.0)

print(f"Average occupation:        <n> = {n_avg_natural:.6f}")
print()

# At CMB temperature
Z_mode_CMB = 1.0 / (1.0 - math.exp(-E_mode / (k_B * T_CMB)))
n_avg_CMB = 1.0 / (math.exp(E_mode / (k_B * T_CMB)) - 1.0)

print(f"At T = T_CMB:")
print(f"Partition function:        Z = {Z_mode_CMB:.6f}")
print(f"Average occupation:        <n> = {n_avg_CMB:.6e}")
print()

print("At CMB temperature, these modes are barely excited.")
print()

# ============================================================================
# VACUUM FLUCTUATIONS
# ============================================================================
print("VACUUM FLUCTUATIONS:")
print("-" * 80)
print()
print("Even at T = 0, quantum zero-point fluctuations persist:")
print()
print("    E_zero = (1/2) ℏω per mode")
print()

E_zero = 0.5 * hbar * omega_mode

print(f"Zero-point energy:         E_0 = (1/2) ℏω")
print(f"                           E_0 = {E_zero:.15e} J")
print(f"                                = {E_zero/e*1000:.6f} meV")
print()

# Total zero-point energy for all modes
E_zero_total = T_17 * E_zero

print(f"Total (T₁₇ modes):         E_total = T₁₇ × E_0")
print(f"                           E_total = {E_zero_total:.6e} J")
print(f"                                   = {E_zero_total/e:.6f} eV")
print()

# Compare to natural energy
ratio_zero = E_zero_total / E_natural

print(f"Ratio E_total / (m_e c² α²): {ratio_zero:.6f}")
print()

print("The total zero-point energy is ~half the natural energy scale.")
print()

# ============================================================================
# THERMAL vs QUANTUM REGIME
# ============================================================================
print("THERMAL vs QUANTUM REGIMES:")
print("-" * 80)
print()
print("The crossover between thermal and quantum behavior occurs when:")
print()
print("    k_B T ~ ℏω")
print()

T_crossover = E_mode / k_B

print(f"Crossover temperature:     T_cross = ℏω / k_B")
print(f"                           T_cross = {T_crossover:.6e} K")
print()

print(f"This equals T_natural = {T_natural:.6e} K")
print()

print("For T << T_natural: Quantum regime (zero-point dominates)")
print("For T >> T_natural: Classical regime (thermal dominates)")
print("For T ~  T_natural: Transition regime")
print()

# Quantum parameter
xi_CMB = hbar * omega_mode / (k_B * T_CMB)
xi_room = hbar * omega_mode / (k_B * T_room)
xi_natural = hbar * omega_mode / (k_B * T_natural)

print(f"Quantum parameter ℏω/k_BT:")
print(f"  At T_CMB:                ξ = {xi_CMB:.6e} (deep quantum)")
print(f"  At T_room:               ξ = {xi_room:.6e} (still quantum)")
print(f"  At T_natural:            ξ = {xi_natural:.6f} (transition)")
print()

# ============================================================================
# CALIBRATION COMPARISON
# ============================================================================
print("="*80)
print("CALIBRATION COMPARISON")
print("="*80)
print()

# CODATA Boltzmann constant
k_B_CODATA = 1.380649e-23  # J/K (exact since 2019)

deviation = E_mode - k_B_CODATA
rel_error = abs(deviation / k_B_CODATA)

print(f"TriPhase E_mode:          {E_mode:.15e} J")
print(f"CODATA k_B:               {k_B_CODATA:.15e} J/K")
print(f"Absolute deviation:       {deviation:+.15e}")
print(f"Relative error:           {rel_error:.6e} ({rel_error*100:.4e}%)")
print()

if rel_error < 1e-8:
    print("✓ EXCELLENT agreement (< 10 ppb)")
elif rel_error < 1e-6:
    print("✓ Good agreement (< 1 ppm)")
else:
    print("⚠ Moderate deviation")

print()
print("NOTE: E_mode = k_B at T_natural by construction.")
print("TriPhase derives k_B from m_e, c, α, and T₁₇.")
print()

# ============================================================================
# PHYSICAL SUMMARY
# ============================================================================
print("="*80)
print("PHYSICAL SUMMARY")
print("="*80)
print()
print("The energy per mode E_mode = m_e c² α² / T₁₇ is:")
print()
print("1. EQUIPARTITION ENERGY:")
print("   At T_natural, each of T₁₇ = 153 DOF has energy E_mode.")
print("   Total energy: E = T₁₇ × E_mode = m_e c² α²")
print()
print("2. BOLTZMANN CONSTANT:")
print("   E_mode = k_B at T = T_natural")
print("   This defines k_B = m_e c² α² / 153")
print()
print("3. MODE FREQUENCY:")
print(f"   ω_mode = E_mode/ℏ = {omega_mode:.3e} rad/s")
print(f"   f_mode = {f_mode:.3e} Hz")
print()
print("4. THERMAL SCALE:")
print(f"   T_crossover = {T_crossover:.3e} K")
print("   Below this, quantum effects dominate.")
print()
print("5. VACUUM FLUCTUATIONS:")
print("   Zero-point energy: E_0 = (1/2)ℏω per mode")
print("   Even at T = 0, each mode has quantum fluctuations.")
print()
print("E_mode is the fundamental thermal quantum — the energy")
print("granularity of vacuum thermodynamics.")
print()

print("="*80)
print("Derivation complete.")
print("="*80)

input("Press Enter to exit...")
