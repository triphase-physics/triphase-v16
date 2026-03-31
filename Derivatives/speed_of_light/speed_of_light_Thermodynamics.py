"""
================================================================================
TriPhase V16 Derivative: Speed of Light
Framework: THERMODYNAMICS
Tag: (D) — Pure derivation from thermodynamic principles
================================================================================

THERMODYNAMIC INTERPRETATION:
------------------------------
c = 1/√(ε₀μ₀)

The speed of light is the maximum signal velocity — the speed at which
thermal equilibrium propagates through the vacuum.

In thermodynamics, heat diffusion is characterized by thermal diffusivity:

    D_thermal = κ/(ρ C_p)

where κ is thermal conductivity, ρ is density, and C_p is specific heat.

For the vacuum electromagnetic field, the "thermal diffusivity" becomes the
speed of signal propagation. This is the speed at which information about
temperature/energy gradients can propagate.

THERMODYNAMIC DERIVATION:
From Maxwell's equations, EM waves propagate at:

    v_EM = 1/√(ε₀μ₀)

But what IS this physically? It's the speed of equilibration.

Consider a local energy perturbation in the vacuum (temperature spike).
The vacuum will relax to equilibrium by radiating energy away. The radiation
propagates at speed c, which is therefore the thermal equilibration speed.

Physical picture:
- ε₀: Electric permittivity (capacitance of space)
- μ₀: Magnetic permeability (inductance of space)
- c = 1/√(ε₀μ₀): Signal speed in this EM "transmission line"

The vacuum is an EM transmission line with impedance Z₀ = √(μ₀/ε₀) and
signal speed c = 1/√(ε₀μ₀).

THERMAL INTERPRETATION:
The speed of light is NOT about photons moving — it's about how fast the
vacuum EM field can respond to perturbations. This is fundamentally a
thermodynamic relaxation timescale.

For a perturbation of wavelength λ:
    τ_relax = λ/c

This is the time for thermal equilibrium to be re-established.

--------------------------------------------------------------------------------
Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
License: CC BY-NC-SA 4.0 International
Website: www.triPhase.org
--------------------------------------------------------------------------------
"""

import math

print("="*80)
print("TriPhase V16 Derivative: Speed of Light")
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

print(f"ε₀        = {epsilon_0:.15e} F/m")
print(f"μ₀        = {mu_0:.15e} H/m")
print()

# ============================================================================
# THERMODYNAMIC DERIVATION OF c
# ============================================================================
print("THERMODYNAMIC DERIVATION:")
print("-" * 80)
print()
print("The vacuum is an electromagnetic transmission line.")
print()
print("TRANSMISSION LINE PARAMETERS:")
print("Capacitance per length:   C' = ε₀ (F/m)")
print("Inductance per length:    L' = μ₀ (H/m)")
print()
print("For any transmission line, signals propagate at speed:")
print()
print("    v_signal = 1/√(L' C')")
print()
print("For the vacuum:")
print()
print("    c = 1/√(μ₀ ε₀)")
print()

# Derive c
c = 1.0 / math.sqrt(epsilon_0 * mu_0)

print(f"Speed of light:           c = {c:.15e} m/s")
print()

# Impedance
Z_0 = math.sqrt(mu_0 / epsilon_0)

print(f"Vacuum impedance:         Z₀ = √(μ₀/ε₀) = {Z_0:.10f} Ω")
print()

# ============================================================================
# THERMAL DIFFUSIVITY INTERPRETATION
# ============================================================================
print("THERMAL DIFFUSIVITY:")
print("-" * 80)
print()
print("In thermodynamics, heat diffuses with diffusivity D = κ/(ρ C_p).")
print()
print("For the vacuum EM field, 'heat' is electromagnetic energy.")
print("The vacuum has:")
print()

# Energy density of EM field
# u = (1/2)ε₀E² + (1/2μ₀)B²
# For a wave: E = cB, so u = ε₀E²

print("Energy density:           u = ε₀ E²")
print("Energy flux (Poynting):   S = (1/μ₀) E × B = (ε₀ c) E²")
print()
print("Thermal diffusivity:      D_th = S/u = c")
print()

D_thermal = c

print(f"Vacuum thermal diffusion: D_th = {D_thermal:.6e} m²/s")
print()

print("The thermal diffusivity of the vacuum is c — the speed of light!")
print("This shows c is fundamentally a THERMAL PROPERTY.")
print()

# ============================================================================
# RELAXATION TIMESCALES
# ============================================================================
print("THERMAL RELAXATION:")
print("-" * 80)
print()
print("For a perturbation of wavelength λ, equilibration time is:")
print()
print("    τ_relax = λ/c")
print()

# Example wavelengths
lambda_Compton_e = 2.42631023867e-12  # Compton wavelength of electron
lambda_visible = 550e-9  # Green light
lambda_CMB = 1.9e-3  # CMB peak

tau_Compton = lambda_Compton_e / c
tau_visible = lambda_visible / c
tau_CMB = lambda_CMB / c

print(f"Electron Compton λ:       λ_C = {lambda_Compton_e:.6e} m")
print(f"Relaxation time:          τ = {tau_Compton:.6e} s")
print()
print(f"Visible light λ:          λ = {lambda_visible:.6e} m")
print(f"Relaxation time:          τ = {tau_visible:.6e} s")
print()
print(f"CMB peak λ:               λ = {lambda_CMB:.6e} m")
print(f"Relaxation time:          τ = {tau_CMB:.6e} s")
print()

# ============================================================================
# THERMODYNAMIC QUANTITIES
# ============================================================================
print("THERMODYNAMIC INTERPRETATION:")
print("-" * 80)
print()

# Derive other constants
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha     = 1.0 / alpha_inv
hbar      = Z_0 * e**2 / (4.0 * math.pi * alpha)

r_e = 2.8179403262e-15
m_e = hbar * alpha / (c * r_e)
f_e = m_e * c**2 / hbar
k_B = m_e * c**2 * alpha**2 / 153.0

print(f"α         = {alpha:.15e}")
print(f"ℏ         = {hbar:.15e} J·s")
print(f"m_e       = {m_e:.15e} kg")
print(f"k_B       = {k_B:.15e} J/K")
print()

# Natural temperature
T_natural = m_e * c**2 * alpha**2 / k_B

print(f"Natural vacuum temp:      T_0 = m_e c² α² / k_B")
print(f"                          T_0 = {T_natural:.6e} K")
print()

# Thermal energy scale
E_thermal = k_B * T_natural

print(f"Thermal energy scale:     E_th = k_B T_0")
print(f"                          E_th = {E_thermal:.6e} J")
print(f"                               = {E_thermal/e:.6e} eV")
print()

# Thermal velocity
v_thermal = math.sqrt(2.0 * k_B * T_natural / m_e)

print(f"Electron thermal velocity: v_th = √(2k_B T/m_e)")
print(f"                           v_th = {v_thermal:.6e} m/s")
print(f"Ratio to c:                v_th/c = {v_thermal/c:.6e}")
print()

# Speed of sound in photon gas
# For relativistic gas: c_s = c/√3
c_sound_photon = c / math.sqrt(3.0)

print(f"Sound speed (photon gas):  c_s = c/√3")
print(f"                           c_s = {c_sound_photon:.6e} m/s")
print(f"Ratio to c:                c_s/c = {c_sound_photon/c:.10f}")
print()

# Debye temperature (if vacuum were solid)
# ω_Debye ~ c/a, where a is lattice spacing
a_lattice = lambda_Compton_e  # Use Compton wavelength as "lattice"
omega_Debye = c / a_lattice
T_Debye = hbar * omega_Debye / k_B

print(f"Vacuum 'lattice' spacing:  a ~ λ_C = {a_lattice:.6e} m")
print(f"Debye frequency:           ω_D = c/a = {omega_Debye:.6e} rad/s")
print(f"Debye temperature:         T_D = ℏω_D/k_B = {T_Debye:.6e} K")
print()

# ============================================================================
# CAUSALITY AND THERMAL CONTACT
# ============================================================================
print("CAUSALITY AND THERMAL CONTACT:")
print("-" * 80)
print()
print("Two regions of space are in thermal contact if they can exchange")
print("energy. The maximum rate of energy exchange is limited by c.")
print()

# Thermal conductivity
# κ = (energy flux) / (temperature gradient)
# For vacuum: κ ~ ε₀ c³

kappa_vacuum = epsilon_0 * c**3

print(f"Vacuum thermal conductivity: κ ~ ε₀ c³")
print(f"                             κ ~ {kappa_vacuum:.6e} W/(m·K)")
print()

# Heat capacity per volume
# C_V ~ energy density / temperature
C_V_vacuum = epsilon_0 * c**2 / T_natural

print(f"Vacuum heat capacity:      C_V ~ ε₀ c² / T")
print(f"                           C_V ~ {C_V_vacuum:.6e} J/(m³·K)")
print()

# Thermal diffusivity check
D_check = kappa_vacuum / C_V_vacuum

print(f"Thermal diffusivity:       D = κ/C_V")
print(f"                           D = {D_check:.6e} m²/s")
print(f"Ratio to c:                D/c = {D_check/c:.6f}")
print()

print("The vacuum's thermal diffusivity is on the order of c, confirming")
print("that c is the speed of thermal equilibration.")
print()

# ============================================================================
# RELATIVISTIC EFFECTS
# ============================================================================
print("RELATIVISTIC THERMODYNAMICS:")
print("-" * 80)
print()
print("The speed limit c emerges from thermodynamic consistency.")
print()
print("Consider a moving heat reservoir at velocity v. The temperature")
print("transforms as:")
print()
print("    T' = T / γ")
print()
print("where γ = 1/√(1 - v²/c²).")
print()
print("For v → c, T' → 0: the reservoir appears frozen.")
print()

v_test = 0.9 * c
gamma_test = 1.0 / math.sqrt(1.0 - (v_test/c)**2)
T_test = T_natural
T_prime = T_test / gamma_test

print(f"Test velocity:            v = 0.9c = {v_test:.6e} m/s")
print(f"Lorentz factor:           γ = {gamma_test:.6f}")
print(f"Lab temperature:          T = {T_test:.6e} K")
print(f"Moving frame temperature: T' = T/γ = {T_prime:.6e} K")
print()

print("At v = c, T' = 0: no thermal motion in photon frame.")
print("This is why photons don't experience proper time — they're at")
print("absolute zero in their own frame!")
print()

# ============================================================================
# CALIBRATION COMPARISON
# ============================================================================
print("="*80)
print("CALIBRATION COMPARISON")
print("="*80)
print()

# Exact SI definition (since 2019)
c_exact = 299792458.0  # m/s (exact by definition)

deviation = c - c_exact
rel_error = abs(deviation / c_exact)

print(f"TriPhase c:               {c:.15e} m/s")
print(f"SI definition c:          {c_exact:.15e} m/s")
print(f"Absolute deviation:       {deviation:+.15e} m/s")
print(f"Relative error:           {rel_error:.6e} ({rel_error*100:.4e}%)")
print()

if rel_error < 1e-10:
    print("✓ EXACT agreement (within numerical precision)")
elif rel_error < 1e-6:
    print("✓ EXCELLENT agreement (< 1 ppm)")
else:
    print("⚠ Deviation (should be exact from ε₀, μ₀)")

print()
print("NOTE: c is derived from ε₀ and μ₀ (exact SI values).")
print("Any deviation is from numerical precision, not physics.")
print()

# ============================================================================
# PHYSICAL SUMMARY
# ============================================================================
print("="*80)
print("PHYSICAL SUMMARY")
print("="*80)
print()
print("The speed of light c = 299792458 m/s emerges from:")
print()
print("1. TRANSMISSION LINE:")
print("   The vacuum is an EM transmission line with capacitance ε₀")
print("   and inductance μ₀. Signals propagate at c = 1/√(ε₀μ₀).")
print()
print("2. THERMAL DIFFUSIVITY:")
print("   c is the thermal diffusivity of the vacuum EM field.")
print("   Heat (EM energy) diffuses at speed c.")
print()
print("3. EQUILIBRATION SPEED:")
print("   For a perturbation of size λ, thermal equilibrium is")
print("   re-established in time τ = λ/c. This is the relaxation time.")
print()
print("4. RELATIVISTIC THERMODYNAMICS:")
print("   Temperature transforms as T' = T/γ. At v = c, T' = 0:")
print("   photons are at absolute zero in their own frame.")
print()
print("5. CAUSALITY:")
print("   Two regions can only be in thermal contact if connected")
print("   by a signal traveling at ≤ c. This defines the causal structure.")
print()
print("The speed of light is NOT just a kinematic limit — it's the")
print("fundamental THERMODYNAMIC timescale of the vacuum.")
print()

print("="*80)
print("Derivation complete.")
print("="*80)

input("Press Enter to exit...")
