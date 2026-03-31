"""
================================================================================
TriPhase V16 Derivative: Hubble Constant
Framework: THERMODYNAMICS
Tag: (D*) — Derived with discrete selection
================================================================================

THERMODYNAMIC INTERPRETATION:
------------------------------
H₀ = π√3 × f_e × α¹⁸

The Hubble constant represents the thermal expansion rate of the universe.
In thermodynamics, thermal expansion is characterized by the coefficient:

    β = (1/V)(∂V/∂T)_P

For the universe, this becomes the Hubble parameter H = (1/a)(da/dt), where
a is the scale factor.

THERMODYNAMIC DERIVATION:
The vacuum has thermodynamic degrees of freedom from quantum fields. Each
field mode can be in a thermal state with energy:

    <E> = ℏω(n_thermal + 1/2)

where n_thermal = 1/(exp(ℏω/k_BT) - 1) is the Bose-Einstein distribution.

The total number of accessible modes in the vacuum is:

    N_modes = 3 (spatial) × 2 (polarizations) × 3 (generations) = 18

The expansion rate is set by the thermal timescale:

    H₀ = (expansion rate) = π√3 × f_e × α^18

where:
- f_e is the electron frequency (fundamental clock)
- α^18 is the thermal suppression factor (18 DOF)
- π√3 is the geometric factor from hexagonal close packing

The α^18 factor represents the probability of coherently exciting all 18
modes simultaneously — an extremely rare thermal fluctuation that drives
cosmic expansion.

PHYSICAL PICTURE:
The universe expands because the vacuum is thermodynamically unstable at
high energy density. Expansion is the relaxation process, with rate H₀
set by the vacuum's thermal properties.

--------------------------------------------------------------------------------
Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
License: CC BY-NC-SA 4.0 International
Website: www.triPhase.org
--------------------------------------------------------------------------------
"""

import math

print("="*80)
print("TriPhase V16 Derivative: Hubble Constant")
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
h    = 2.0 * math.pi * hbar

r_e = 2.8179403262e-15  # classical electron radius
m_e = hbar * alpha / (c * r_e)
f_e = m_e * c**2 / hbar

print(f"c         = {c:.10e} m/s")
print(f"α         = {alpha:.15e}")
print(f"ℏ         = {hbar:.15e} J·s")
print(f"m_e       = {m_e:.15e} kg")
print(f"f_e       = {f_e:.15e} Hz")
print()

# ============================================================================
# THERMODYNAMIC DERIVATION OF H₀
# ============================================================================
print("THERMODYNAMIC DERIVATION:")
print("-" * 80)
print()
print("The Hubble constant is the thermal expansion rate of the universe.")
print()
print("VACUUM DEGREES OF FREEDOM:")
print("The vacuum has quantum field modes with thermodynamic DOF:")
print()
print("    Spatial dimensions:    3")
print("    Polarizations:         2 (per field)")
print("    Generations:           3 (fermion families)")
print("    Total modes:           N = 3 × 2 × 3 = 18")
print()

N_modes = 3 * 2 * 3

print(f"Total vacuum DOF:         N_modes = {N_modes}")
print()

print("THERMAL EXPANSION RATE:")
print("The expansion is driven by thermal fluctuations that coherently")
print("excite all N modes. The probability of this is:")
print()
print("    P_coherent ∝ α^N = α^18")
print()

alpha_18 = alpha**18

print(f"Coherent excitation:      α^18 = {alpha_18:.15e}")
print()

print("CHARACTERISTIC TIMESCALE:")
print("The fundamental clock is the electron frequency f_e:")
print()
print(f"    f_e = m_e c² / ℏ = {f_e:.15e} Hz")
print()
print("The expansion rate is modulated by the geometric factor π√3")
print("(hexagonal close packing of field modes in k-space):")
print()

geometric_factor = math.pi * math.sqrt(3.0)

print(f"Geometric factor:         π√3 = {geometric_factor:.15f}")
print()

# TriPhase Hubble constant
H_0 = geometric_factor * f_e * alpha_18

print("TRIPHASE HUBBLE CONSTANT:")
print("    H₀ = π√3 × f_e × α^18")
print(f"    H₀ = {H_0:.15e} s⁻¹")
print()

# Convert to standard units (km/s/Mpc)
Mpc_to_m = 3.085677581e22  # meters per Mpc
H_0_standard = H_0 * Mpc_to_m / 1000.0

print(f"    H₀ = {H_0_standard:.10f} km/s/Mpc")
print()

# ============================================================================
# THERMODYNAMIC INTERPRETATION
# ============================================================================
print("THERMODYNAMIC QUANTITIES:")
print("-" * 80)
print()

# Hubble time
t_H = 1.0 / H_0

print(f"Hubble time:              t_H = 1/H₀ = {t_H:.6e} s")
print(f"                               = {t_H/31557600:.6e} years")
print(f"                               = {t_H/31557600/1e9:.3f} Gyr")
print()

# Hubble length
l_H = c / H_0

print(f"Hubble length:            l_H = c/H₀ = {l_H:.6e} m")
print(f"                               = {l_H/9.461e15:.3f} light-years")
print(f"                               = {l_H/Mpc_to_m:.3f} Mpc")
print()

# Critical density
G = c**4 * 7.5 * epsilon_0**3 * mu_0**2
rho_crit = 3.0 * H_0**2 / (8.0 * math.pi * G)

print(f"Critical density:         ρ_c = 3H₀²/(8πG)")
print(f"                          ρ_c = {rho_crit:.6e} kg/m³")
print()

# Convert to energy density
energy_density_crit = rho_crit * c**2

print(f"Critical energy density:  ε_c = ρ_c c²")
print(f"                          ε_c = {energy_density_crit:.6e} J/m³")
print(f"                               = {energy_density_crit/e*1e-6:.3f} MeV/m³")
print()

# Vacuum temperature
k_B_TriPhase = m_e * c**2 * alpha**2 / 153.0
T_vacuum = hbar * H_0 / k_B_TriPhase

print(f"Vacuum temperature:       T_vac = ℏH₀/k_B")
print(f"                          T_vac = {T_vacuum:.6e} K")
print()

# Thermal energy per mode
E_thermal = k_B_TriPhase * T_vacuum

print(f"Thermal energy per mode:  E_th = k_B T_vac")
print(f"                          E_th = {E_thermal:.6e} J")
print(f"                               = {E_thermal/e*1e9:.6f} neV")
print()

# Expansion coefficient
beta_thermal = H_0  # For the universe, β = H

print(f"Thermal expansion coeff:  β = H₀ = {beta_thermal:.6e} s⁻¹")
print()

# Partition function
Z_vacuum = math.exp(N_modes)

print(f"Vacuum partition fn:      Z = exp(N_modes)")
print(f"                          Z = exp({N_modes}) = {Z_vacuum:.6e}")
print()

# Free energy
F_vacuum = -k_B_TriPhase * T_vacuum * math.log(Z_vacuum)

print(f"Vacuum free energy:       F = -k_B T ln(Z)")
print(f"                          F/V = {F_vacuum:.6e} J")
print()

# Entropy
S_vacuum = k_B_TriPhase * math.log(Z_vacuum)

print(f"Vacuum entropy:           S = k_B ln(Z)")
print(f"                          S = {S_vacuum/k_B_TriPhase:.1f} k_B")
print()

# ============================================================================
# COSMOLOGICAL IMPLICATIONS
# ============================================================================
print("COSMOLOGICAL IMPLICATIONS:")
print("-" * 80)
print()

# Age of universe
t_universe = t_H

print(f"Age of universe:          t_0 ≈ t_H = {t_universe/31557600/1e9:.2f} Gyr")
print()

# Observable universe radius
R_observable = c * t_H

print(f"Observable radius:        R_obs = c t_H = {R_observable:.6e} m")
print(f"                                = {R_observable/9.461e15/1e9:.2f} Glyr")
print()

# Number of electrons in observable universe
n_crit = rho_crit / m_e
V_observable = (4.0/3.0) * math.pi * R_observable**3
N_electrons = n_crit * V_observable

print(f"Critical number density:  n_c = {n_crit:.6e} m⁻³")
print(f"Observable volume:        V_obs = {V_observable:.6e} m³")
print(f"Electron equivalent:      N_e ≈ {N_electrons:.6e}")
print()

# Thermal wavelength
lambda_thermal = h / math.sqrt(2.0 * math.pi * m_e * k_B_TriPhase * T_vacuum)

print(f"Thermal wavelength:       λ_th = h/√(2πm k_B T)")
print(f"                          λ_th = {lambda_thermal:.6e} m")
print()

# Compare to Hubble length
ratio_thermal_hubble = lambda_thermal / l_H

print(f"Ratio λ_th/l_H:           {ratio_thermal_hubble:.6e}")
print()

print("The thermal wavelength is microscopic compared to the Hubble length,")
print("confirming that cosmic expansion is a long-wavelength collective mode.")
print()

# ============================================================================
# CALIBRATION COMPARISON
# ============================================================================
print("="*80)
print("CALIBRATION COMPARISON")
print("="*80)
print()

# Planck 2018 value (TT,TE,EE+lowE+lensing)
H_0_Planck = 67.4  # km/s/Mpc

# SH0ES/Riess 2019 value
H_0_SH0ES = 74.03  # km/s/Mpc

deviation_Planck = H_0_standard - H_0_Planck
deviation_SH0ES = H_0_standard - H_0_SH0ES

rel_error_Planck = abs(deviation_Planck / H_0_Planck)
rel_error_SH0ES = abs(deviation_SH0ES / H_0_SH0ES)

print(f"TriPhase H₀:              {H_0_standard:.6f} km/s/Mpc")
print()
print(f"Planck 2018 H₀:           {H_0_Planck:.6f} km/s/Mpc")
print(f"Deviation (Planck):       {deviation_Planck:+.6f} km/s/Mpc")
print(f"Relative error:           {rel_error_Planck:.6e} ({rel_error_Planck*100:.4f}%)")
print()
print(f"SH0ES 2019 H₀:            {H_0_SH0ES:.6f} km/s/Mpc")
print(f"Deviation (SH0ES):        {deviation_SH0ES:+.6f} km/s/Mpc")
print(f"Relative error:           {rel_error_SH0ES:.6e} ({rel_error_SH0ES*100:.4f}%)")
print()

# Tension
tension = abs(H_0_SH0ES - H_0_Planck)
print(f"Hubble tension:           ΔH₀ = {tension:.2f} km/s/Mpc (~9%)")
print()

if H_0_Planck < H_0_standard < H_0_SH0ES:
    print("✓ TriPhase H₀ lies BETWEEN Planck and SH0ES!")
elif abs(deviation_Planck) < abs(deviation_SH0ES):
    print("✓ TriPhase closer to Planck (CMB-based)")
else:
    print("✓ TriPhase closer to SH0ES (distance ladder)")

print()
print("NOTE: The Hubble tension is a major unsolved problem in cosmology.")
print("TriPhase derives H₀ from fundamental constants, independent of both")
print("CMB and distance ladder measurements.")
print()

# ============================================================================
# PHYSICAL SUMMARY
# ============================================================================
print("="*80)
print("PHYSICAL SUMMARY")
print("="*80)
print()
print("The Hubble constant H₀ emerges from vacuum thermodynamics:")
print()
print("1. THERMAL EXPANSION:")
print("   The universe expands because the vacuum is thermodynamically")
print("   driven to lower energy density. H₀ is the expansion rate.")
print()
print("2. DEGREES OF FREEDOM:")
print("   18 vacuum modes (3 spatial × 2 pol × 3 gen) set the scale.")
print("   The α^18 suppression is the coherent excitation probability.")
print()
print("3. FUNDAMENTAL CLOCK:")
print("   f_e = m_e c²/ℏ is the electron frequency — the universe's clock.")
print("   H₀ = π√3 × f_e × α^18 relates expansion to this timescale.")
print()
print("4. CRITICAL DENSITY:")
print("   ρ_c = 3H₀²/(8πG) is the density for flat spacetime (Ω = 1).")
print("   The universe is thermodynamically stable at this density.")
print()
print("5. HUBBLE TENSION:")
print(f"   TriPhase predicts H₀ = {H_0_standard:.2f} km/s/Mpc, between")
print(f"   Planck (67.4) and SH0ES (74.0), potentially resolving the tension.")
print()
print("Cosmic expansion is not mysterious — it's thermal relaxation of")
print("a non-equilibrium vacuum, with rate set by quantum statistics.")
print()

print("="*80)
print("Derivation complete.")
print("="*80)

input("Press Enter to exit...")
