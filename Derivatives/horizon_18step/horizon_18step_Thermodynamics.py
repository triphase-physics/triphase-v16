"""
================================================================================
TriPhase V16 - Hubble Horizon (18-Step) Derivative
Framework: THERMODYNAMICS
Tag: (D*) - Derived with discrete selection
================================================================================

THERMODYNAMICS FRAMEWORK:
Interprets each physical quantity through statistical mechanics / thermodynamic
concepts: partition functions Z and free energy F = -k_BT ln(Z), entropy
S = -∂F/∂T, equipartition theorem (½k_BT per degree of freedom), Boltzmann
distributions, Maxwell-Boltzmann statistics, phase transitions, order parameters,
critical phenomena, equations of state, thermodynamic potentials (U, H, F, G),
degrees of freedom counting, mode counting, heat capacity, specific heat,
adiabatic processes, black-body radiation, Planck distribution, thermodynamic
stability conditions, Stefan-Boltzmann law, Wien displacement, chemical potential,
Gibbs free energy.

PHYSICAL DERIVATION:
The Hubble horizon R_H = c/H₀ is the maximum distance over which thermal
equilibrium has been established in the universe. It is the thermodynamic
correlation length — the largest scale over which information (entropy) has
been exchanged since the Big Bang.

Thermodynamically, R_H divides the universe into two regimes:
  - r < R_H: Causally connected, thermally equilibrated (entropy maximized)
  - r > R_H: Causally disconnected, potentially out of equilibrium

The Hubble constant H₀ = π√3 × f_e × α^18 emerges from TriPhase's α^18 scaling,
connecting the cosmological expansion rate to atomic physics.

The horizon is the thermodynamic "sphere of influence" — the volume within which
the universe has had time to thermalize. Beyond R_H, the universe may be in a
different thermodynamic state.

Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# ============================================================================
# STANDARD ANCHOR CHAIN
# ============================================================================
epsilon_0 = 8.8541878128e-12   # F/m
mu_0      = 1.25663706212e-6   # H/m
e         = 1.602176634e-19    # C (exact SI)
c         = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0       = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha     = 1.0 / alpha_inv
hbar      = Z_0 * e**2 / (4.0 * math.pi * alpha)
h         = 2.0 * math.pi * hbar
G         = c**4 * 7.5 * epsilon_0**3 * mu_0**2
r_e       = 2.8179403262e-15
m_e       = hbar * alpha / (c * r_e)
f_e       = m_e * c**2 / hbar
T_17      = 17 * 18 // 2   # = 153
mp_me     = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p       = m_e * mp_me
H_0       = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r      = c**4 / (8.0 * math.pi * G)

# ============================================================================
# THERMODYNAMIC DERIVATION - HUBBLE HORIZON
# ============================================================================

print("=" * 80)
print("HUBBLE HORIZON - THERMODYNAMICS FRAMEWORK")
print("=" * 80)
print()
print("THERMODYNAMIC INTERPRETATION:")
print("The Hubble horizon R_H = c/H₀ is the thermal correlation length.")
print("It is the maximum distance over which the universe has thermalized.")
print("Beyond R_H, causal connection and thermal equilibrium are not guaranteed.")
print()

# Hubble constant from TriPhase
print("HUBBLE CONSTANT DERIVATION:")
print(f"  Electron frequency f_e      = {f_e:.6e} Hz")
print(f"  Fine structure α            = {alpha:.10f}")
print(f"  α^18                        = {alpha**18:.6e}")
print()
print(f"  H₀ = π√3 × f_e × α^18")
print(f"     = {math.pi * math.sqrt(3.0):.6f} × {f_e:.6e} × {alpha**18:.6e}")
print(f"     = {H_0:.6e} Hz")
print()

# Convert to standard units (km/s/Mpc)
H_0_SI = H_0  # Already in Hz = 1/s
km_per_m = 1e-3
Mpc_per_m = 1.0 / 3.085677581e22  # 1 Mpc in meters
H_0_standard = H_0_SI / (km_per_m / Mpc_per_m)  # km/s/Mpc

print(f"  H₀ (standard units)         = {H_0_standard:.3f} km/s/Mpc")
print()

# Hubble horizon
R_H = c / H_0
R_H_Mpc = R_H * Mpc_per_m
R_H_Gly = R_H / 9.461e24  # Convert to billion light-years

print("HUBBLE HORIZON:")
print(f"  R_H = c / H₀")
print(f"      = {c:.6e} m/s / {H_0:.6e} Hz")
print(f"      = {R_H:.6e} m")
print(f"      = {R_H_Mpc:.2f} Mpc")
print(f"      = {R_H_Gly:.2f} Gly (billion light-years)")
print()

# Thermodynamic correlation length
print("=" * 80)
print("THERMODYNAMIC CORRELATION LENGTH")
print("=" * 80)
print()

print("CAUSAL STRUCTURE:")
print(f"  The Hubble horizon divides spacetime into causal regions.")
print(f"  Light travel time across R_H:")
print(f"    t_H = R_H / c = 1 / H₀ = {1.0/H_0:.6e} seconds")
print(f"                          = {1.0/H_0 / (365.25*24*3600):.2e} years")
print(f"                          ≈ {1.0/H_0 / (1e9*365.25*24*3600):.2f} Gyr")
print()

# Thermalization timescale
t_H_Gyr = 1.0 / H_0 / (1e9 * 365.25 * 24 * 3600)
print("THERMALIZATION:")
print(f"  Within R_H: The universe has had {t_H_Gyr:.1f} Gyr to exchange")
print(f"              entropy and reach thermal equilibrium.")
print(f"  Beyond R_H: Regions may have different temperatures, densities,")
print(f"              or even different thermodynamic phases.")
print()

# Entropy and information
S_Hubble_kB = 2.0 * math.pi * (R_H / 1.616255e-35)**2  # Holographic entropy (order of magnitude)
print("HOLOGRAPHIC ENTROPY:")
print(f"  The Hubble horizon has surface area A_H = 4πR_H²")
print(f"  Holographic entropy bound:")
print(f"    S_H ≤ (A_H / 4ℓ_P²) × k_B")
print(f"  where ℓ_P is the Planck length.")
print(f"  This entropy S_H ∼ {S_Hubble_kB:.2e} k_B is the maximum")
print(f"  thermodynamic information content of the observable universe.")
print()

# Temperature scales
T_CMB_K = 2.725  # CMB temperature today
T_CMB_eV = T_CMB_K * 8.617333e-5  # Convert to eV
k_B = 1.380649e-23  # J/K

print("THERMAL STATE:")
print(f"  Current CMB temperature     = {T_CMB_K:.3f} K")
print(f"                              = {T_CMB_eV*1e3:.3f} meV")
print(f"  Thermal wavelength λ_th     = hc/(k_BT)")
print(f"                              = {h*c/(k_B*T_CMB_K):.3e} m")
print(f"                              ≈ {h*c/(k_B*T_CMB_K)*1e3:.1f} mm")
print()
print(f"  The CMB photons are in thermal equilibrium within R_H.")
print(f"  This equilibrium was established at recombination (T ≈ 3000 K).")
print()

# Horizon problem
print("=" * 80)
print("HORIZON PROBLEM (THERMODYNAMIC PUZZLE)")
print("=" * 80)
print()

print("THE PUZZLE:")
print(f"  The CMB has the same temperature T = {T_CMB_K:.3f} K to 1 part in 10⁵")
print(f"  in all directions. This implies thermal equilibrium.")
print()
print(f"  BUT: Opposite sides of the sky (separated by 2R_H) were never")
print(f"       in causal contact in standard Big Bang cosmology.")
print()
print(f"  How did they reach the same temperature (thermalize)?")
print()

print("THERMODYNAMIC EXPLANATION (INFLATION):")
print(f"  Inflationary cosmology solves this by expanding R_H faster than c:")
print(f"    During inflation: H ≈ constant → R_H = c/H fixed")
print(f"                      But scale factor a(t) ∝ exp(Ht) expands")
print(f"    Result: Regions now separated by > R_H were once within")
print(f"            a much smaller horizon → had time to thermalize")
print()
print(f"  This is a thermodynamic causality argument for inflation.")
print()

# Free energy and horizon
print("=" * 80)
print("FREE ENERGY AND THE HORIZON")
print("=" * 80)
print()

print("COSMOLOGICAL FREE ENERGY:")
print(f"  The horizon R_H sets the scale of the cosmological free energy.")
print(f"  Energy within the Hubble volume:")
print(f"    E_H = ρ_c × V_H = ρ_c × (4π/3)R_H³")
print()

# Critical density
rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)
V_H = (4.0/3.0) * math.pi * R_H**3
E_H = rho_c * V_H
M_H = E_H / c**2

print(f"  Critical density ρ_c        = {rho_c:.6e} kg/m³")
print(f"  Hubble volume V_H           = {V_H:.6e} m³")
print(f"  Energy E_H                  = {E_H:.6e} J")
print(f"  Equivalent mass M_H         = {M_H:.6e} kg")
print(f"                              = {M_H/1.989e30:.2e} M_☉")
print()

print("THERMODYNAMIC EQUILIBRIUM:")
print(f"  The horizon R_H is the maximum scale for equilibrium thermodynamics.")
print(f"  Statistical mechanics applies within R_H (canonical ensemble).")
print(f"  Beyond R_H, we enter the realm of non-equilibrium cosmology.")
print()

# ============================================================================
# CALIBRATION COMPARISON
# ============================================================================

print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)

# Planck 2018 Hubble constant (tension exists!)
H_0_Planck = 67.4  # km/s/Mpc (CMB-based)
H_0_SH0ES = 73.0  # km/s/Mpc (local distance ladder)

R_H_Planck_Mpc = c / ((H_0_Planck * 1e3) / 3.085677581e22)
R_H_SH0ES_Mpc = c / ((H_0_SH0ES * 1e3) / 3.085677581e22)

deviation_Planck = ((H_0_standard - H_0_Planck) / H_0_Planck) * 100
deviation_SH0ES = ((H_0_standard - H_0_SH0ES) / H_0_SH0ES) * 100

print()
print(f"TriPhase derived H₀         = {H_0_standard:.3f} km/s/Mpc")
print(f"Planck 2018 (CMB)           = {H_0_Planck:.1f} ± 0.5 km/s/Mpc")
print(f"SH0ES (local)               = {H_0_SH0ES:.1f} ± 1.0 km/s/Mpc")
print()
print(f"Deviation from Planck       = {deviation_Planck:+.2f}%")
print(f"Deviation from SH0ES        = {deviation_SH0ES:+.2f}%")
print()
print(f"TriPhase R_H                = {R_H_Mpc:.2f} Mpc = {R_H_Gly:.2f} Gly")
print(f"Planck-implied R_H          = {R_H_Planck_Mpc*Mpc_per_m*1e-6:.2f} Gly")
print(f"SH0ES-implied R_H           = {R_H_SH0ES_Mpc*Mpc_per_m*1e-6:.2f} Gly")
print()
print("NOTE: Hubble tension is a major unsolved problem in cosmology.")
print("      TriPhase value is anchored to atomic physics (α^18 scaling).")
print()

print("THERMODYNAMIC SIGNIFICANCE:")
print("  The Hubble horizon is the ultimate thermodynamic boundary.")
print("  It defines the largest system for which statistical mechanics applies.")
print("  The horizon entropy S_H is the information content of the cosmos.")
print()

print("=" * 80)
print("END HUBBLE HORIZON THERMODYNAMICS DERIVATIVE")
print("=" * 80)

input("Press Enter to exit...")
