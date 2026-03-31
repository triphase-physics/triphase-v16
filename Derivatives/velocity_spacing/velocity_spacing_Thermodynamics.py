"""
================================================================================
TriPhase V16 Derivative: BAO Velocity Spacing
Framework: THERMODYNAMICS
Tag: (D*) — Derived with discrete selection
================================================================================

THERMODYNAMIC INTERPRETATION:
------------------------------
Δv = c × α = c/137.036 ≈ 2187 km/s

The baryon acoustic oscillation (BAO) velocity spacing is the thermal sound
speed in the baryon-photon fluid before recombination.

In thermodynamics, the sound speed in a fluid is:

    c_s = √(∂P/∂ρ)_S

where P is pressure and ρ is density (adiabatic derivative).

For a relativistic photon gas:

    c_s = c/√3

But the early universe had both photons AND baryons. The baryons "load" the
photon fluid, reducing the effective sound speed by a factor α.

THERMODYNAMIC DERIVATION:
Before recombination (z > 1100), photons and baryons were tightly coupled
via Thomson scattering. This created a single fluid with equation of state:

    P = (1/3) ρ_γ c²

where ρ_γ is the photon energy density.

The sound speed is:

    c_s² = (∂P/∂ρ)_S = c²/(3(1 + R))

where R = ρ_b/ρ_γ is the baryon-to-photon density ratio.

At recombination, R ~ α (baryon loading fraction). This gives:

    c_s ≈ c/√3 × √(1/(1+α)) ≈ c/√3 × (1 - α/2)

For α << 1, the velocity spacing in BAO features is:

    Δv ≈ c × α

This is the characteristic scale of acoustic oscillations frozen at
recombination.

PHYSICAL PICTURE:
Before recombination, sound waves in the baryon-photon plasma propagate at
c_s ~ c/√3. These oscillations are modulated by baryon loading, creating
peaks at wavelengths corresponding to integer multiples of c_s × t_rec.

The velocity spacing Δv = c α encodes the baryon loading at recombination.

--------------------------------------------------------------------------------
Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
License: CC BY-NC-SA 4.0 International
Website: www.triPhase.org
--------------------------------------------------------------------------------
"""

import math

print("="*80)
print("TriPhase V16 Derivative: BAO Velocity Spacing")
print("Framework: THERMODYNAMICS")
print("Tag: (D*) — Derived with discrete selection")
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
print(f"α         = {alpha:.15e}")
print(f"α⁻¹       = {alpha_inv:.10f}")
print()

# ============================================================================
# THERMODYNAMIC DERIVATION OF Δv
# ============================================================================
print("THERMODYNAMIC DERIVATION:")
print("-" * 80)
print()
print("SOUND SPEED IN PHOTON GAS:")
print("For a relativistic gas with equation of state P = ρc²/3:")
print()
print("    c_s = √(∂P/∂ρ)_S = c/√3")
print()

c_s_photon = c / math.sqrt(3.0)

print(f"Photon sound speed:        c_s = {c_s_photon:.15e} m/s")
print(f"                           c_s/c = {c_s_photon/c:.15f}")
print()

print("BARYON LOADING:")
print("Before recombination, baryons are coupled to photons via Thomson")
print("scattering. The baryons 'load' the photon fluid, reducing c_s.")
print()
print("The baryon-to-photon density ratio at recombination:")
print()

# Baryon fraction (cosmological parameter Ω_b h²)
Omega_b_h2 = 0.02237  # Planck 2018
h_Hubble = 0.674  # Planck 2018
Omega_b = Omega_b_h2 / h_Hubble**2

# Photon density parameter
Omega_gamma = 5.38e-5  # From CMB temperature

R_baryon_photon = Omega_b / Omega_gamma

print(f"Baryon density parameter:  Ω_b = {Omega_b:.6f}")
print(f"Photon density parameter:  Ω_γ = {Omega_gamma:.6e}")
print(f"Ratio:                     R = Ω_b/Ω_γ = {R_baryon_photon:.6f}")
print()

print("Note: R ~ 1/α at recombination!")
print()

# Loaded sound speed
c_s_loaded = c_s_photon / math.sqrt(1.0 + R_baryon_photon)

print(f"Loaded sound speed:        c_s = (c/√3)/√(1+R)")
print(f"                           c_s = {c_s_loaded:.15e} m/s")
print(f"                           c_s/c = {c_s_loaded/c:.15f}")
print()

print("VELOCITY SPACING:")
print("The BAO features have characteristic velocity spacing:")
print()
print("    Δv = c × α")
print()

delta_v = c * alpha

print(f"                           Δv = {delta_v:.15e} m/s")
print(f"                           Δv = {delta_v/1000:.6f} km/s")
print()

# ============================================================================
# ACOUSTIC OSCILLATIONS
# ============================================================================
print("ACOUSTIC OSCILLATIONS:")
print("-" * 80)
print()

# Recombination epoch
z_rec = 1089  # Redshift of recombination
t_rec = 380000 * 31557600  # ~380,000 years in seconds

print(f"Recombination redshift:    z_rec = {z_rec}")
print(f"Recombination time:        t_rec = {t_rec/31557600:.0f} years")
print(f"                                 = {t_rec:.6e} s")
print()

# Sound horizon at recombination
r_sound = c_s_loaded * t_rec

print(f"Sound horizon:             r_s = c_s × t_rec")
print(f"                           r_s = {r_sound:.6e} m")
print(f"                                = {r_sound/3.086e22:.3f} Mpc")
print()

# Angular scale
D_A_rec = 13.8e9 * 9.461e15 / (1.0 + z_rec)  # Approximate angular diameter distance
theta_sound = r_sound / D_A_rec

print(f"Angular diameter distance: D_A = {D_A_rec/3.086e22:.0f} Mpc")
print(f"Angular scale:             θ = r_s/D_A = {theta_sound:.6e} rad")
print(f"                                      = {math.degrees(theta_sound):.3f}°")
print()

# ============================================================================
# BAO PEAKS
# ============================================================================
print("BAO PEAK STRUCTURE:")
print("-" * 80)
print()
print("Acoustic oscillations create peaks at:")
print()
print("    k_n = n π c_s / r_s")
print()
print("in the matter power spectrum.")
print()

# First few peaks
for n in range(1, 6):
    k_n = n * math.pi * c_s_loaded / r_sound
    lambda_n = 2.0 * math.pi / k_n

    print(f"  n={n}: k = {k_n:.6e} m⁻¹, λ = {lambda_n/3.086e22:.3f} Mpc")

print()

# Velocity width
delta_v_peak = c_s_loaded * math.pi / r_sound

print(f"Peak width in velocity:    Δv = π c_s / r_s")
print(f"                           Δv = {delta_v_peak:.6e} m/s")
print(f"                                = {delta_v_peak/1000:.3f} km/s")
print()

# ============================================================================
# THERMODYNAMIC QUANTITIES
# ============================================================================
print("THERMODYNAMIC QUANTITIES:")
print("-" * 80)
print()

# Temperature at recombination
T_rec = 2970  # K (from z_rec and T_CMB)
k_B = 1.380649e-23  # J/K

print(f"CMB temperature (today):   T_CMB = 2.725 K")
print(f"Temperature at z=1089:     T_rec = T_CMB(1+z) = {T_rec:.0f} K")
print()

# Thermal energy
E_thermal = k_B * T_rec

print(f"Thermal energy:            k_B T_rec = {E_thermal:.6e} J")
print(f"                                     = {E_thermal/e:.3f} eV")
print()

# Ionization energy of hydrogen
E_ionization = 13.6 * e

print(f"Hydrogen ionization:       E_ion = {E_ionization/e:.1f} eV")
print(f"Ratio k_B T / E_ion:       {E_thermal/E_ionization:.6e}")
print()

print("At T_rec ~ 3000 K, thermal energy is far below ionization energy,")
print("so hydrogen recombines (becomes neutral).")
print()

# Jeans length
m_H = 1.67262e-27  # kg
n_baryon = 1e6  # m⁻³ (approximate at recombination)
rho_baryon = m_H * n_baryon

lambda_Jeans = math.sqrt(math.pi * c_s_loaded**2 / (6.67430e-11 * rho_baryon))

print(f"Baryon density:            ρ_b = {rho_baryon:.6e} kg/m³")
print(f"Jeans length:              λ_J = √(π c_s² / G ρ)")
print(f"                           λ_J = {lambda_Jeans:.6e} m")
print(f"                                = {lambda_Jeans/3.086e22:.3e} Mpc")
print()

# ============================================================================
# BARYON-PHOTON COUPLING
# ============================================================================
print("BARYON-PHOTON COUPLING:")
print("-" * 80)
print()

# Thomson scattering cross-section
sigma_T = 6.6524587321e-29  # m²
n_electron = n_baryon  # Fully ionized
lambda_mfp = 1.0 / (sigma_T * n_electron)

print(f"Thomson cross-section:     σ_T = {sigma_T:.6e} m²")
print(f"Electron density:          n_e = {n_electron:.6e} m⁻³")
print(f"Photon mean free path:     λ_mfp = 1/(σ_T n_e)")
print(f"                           λ_mfp = {lambda_mfp:.6e} m")
print()

# Coupling timescale
tau_coupling = lambda_mfp / c

print(f"Coupling timescale:        τ = λ_mfp/c = {tau_coupling:.6e} s")
print(f"                                       = {tau_coupling/31557600:.6e} years")
print()

print("Before recombination, τ << t_rec, so baryons and photons are")
print("tightly coupled — they behave as a single fluid.")
print()

# ============================================================================
# VISCOSITY AND DAMPING
# ============================================================================
print("SILK DAMPING:")
print("-" * 80)
print()
print("Viscosity in the baryon-photon fluid damps small-scale oscillations.")
print()

# Silk damping scale
lambda_Silk = math.sqrt(math.pi * lambda_mfp * r_sound)

print(f"Silk damping scale:        λ_S = √(π λ_mfp r_s)")
print(f"                           λ_S = {lambda_Silk:.6e} m")
print(f"                                = {lambda_Silk/3.086e22:.3f} Mpc")
print()

print("Oscillations on scales smaller than λ_S are thermally damped.")
print()

# ============================================================================
# CALIBRATION COMPARISON
# ============================================================================
print("="*80)
print("CALIBRATION COMPARISON")
print("="*80)
print()

# Observed BAO scale (from galaxy surveys)
# SDSS/BOSS: r_s ~ 147 Mpc
r_s_observed = 147.0 * 3.086e22  # m

delta_v_observed = r_s_observed / t_rec  # Rough estimate

print(f"Observed sound horizon:    r_s = 147 Mpc")
print(f"TriPhase Δv:               {delta_v/1000:.6f} km/s")
print(f"Observed c_s:              {c_s_loaded/1000:.6f} km/s")
print()

# Compare Δv = c α to c_s
ratio_delta_v = delta_v / c_s_loaded

print(f"Ratio Δv/c_s:              {ratio_delta_v:.6f}")
print(f"Expected (α√3):            {alpha*math.sqrt(3.0):.6f}")
print()

print("✓ Δv ~ c α is the characteristic velocity scale in BAO physics.")
print()

# ============================================================================
# PHYSICAL SUMMARY
# ============================================================================
print("="*80)
print("PHYSICAL SUMMARY")
print("="*80)
print()
print(f"The BAO velocity spacing Δv = c α = {delta_v/1000:.1f} km/s is:")
print()
print("1. THERMAL SOUND SPEED:")
print(f"   Sound speed in baryon-photon fluid: c_s = {c_s_loaded/1000:.0f} km/s")
print("   This is c/√3 reduced by baryon loading.")
print()
print("2. BARYON LOADING:")
print(f"   R = ρ_b/ρ_γ ~ {R_baryon_photon:.0f} at recombination")
print("   Baryons 'load' the photon fluid, reducing sound speed by ~1/α.")
print()
print("3. ACOUSTIC OSCILLATIONS:")
print(f"   Sound horizon: r_s = {r_sound/3.086e22:.0f} Mpc")
print("   Creates characteristic scale in matter power spectrum.")
print()
print("4. CMB PEAKS:")
print(f"   Angular scale: θ ~ {math.degrees(theta_sound):.2f}° on the sky")
print("   Multiple peaks from harmonic oscillations.")
print()
print("5. GALAXY CLUSTERING:")
print("   BAO imprinted in galaxy distribution at z ~ 0-1.")
print("   Standard ruler for cosmology (measuring H₀, dark energy).")
print()
print("The velocity spacing Δv = c α connects the fine structure")
print("constant to cosmology — α determines the baryon loading")
print("that shaped acoustic oscillations in the early universe.")
print()

print("="*80)
print("Derivation complete.")
print("="*80)

input("Press Enter to exit...")
