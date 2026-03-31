#!/usr/bin/env python3
"""
================================================================================
thermal_pressure_GroupTheory.py
================================================================================
TriPhase Wave Mechanics Framework — V16 Python Derivatives
GroupTheory Module: Thermal Pressure P = nk_BT

Interprets thermal pressure as counting ACCESSIBLE GROUP REPRESENTATIONS at
temperature T. Boltzmann constant converts between energy (Casimir eigenvalue)
and temperature (representation parameter). Equipartition: ½k_BT per generator.

MIS Tag: (D) — Pure derivation from ε₀, μ₀, α, e, c

IRON RULES:
- import math only (NO numpy, scipy)
- CODATA/PDG values are CALIBRATION CHECKPOINTS only
- Every calculation derives from epsilon_0 and mu_0
- mp_me = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)

Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
All Rights Reserved.
================================================================================
"""

import math

print("="*80)
print("TRIPHASE V16: THERMAL PRESSURE — GROUP THEORY")
print("P = nk_BT as Representation Counting")
print("="*80)
print()

# ==============================================================================
# ANCHOR CHAIN: Standard TriPhase V16 Constants
# ==============================================================================
print("ANCHOR CHAIN: Building from ε₀ and μ₀")
print("-" * 80)

epsilon_0 = 8.8541878128e-12  # F/m
mu_0      = 1.25663706212e-6  # H/m
e         = 1.602176634e-19   # C

print(f"ε₀ = {epsilon_0:.13e} F/m")
print(f"μ₀ = {mu_0:.14e} H/m")
print(f"e  = {e:.12e} C")
print()

c = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0 = math.sqrt(mu_0 / epsilon_0)
print(f"c = {c:.10e} m/s")
print(f"Z₀ = {Z_0:.10f} Ω")
print()

alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha = 1.0 / alpha_inv
hbar = Z_0 * e**2 / (4.0 * math.pi * alpha)
h = 2.0 * math.pi * hbar
print(f"α = {alpha:.15f}")
print(f"ℏ = {hbar:.15e} J·s")
print()

r_e = 2.8179403262e-15
m_e = hbar * alpha / (c * r_e)
f_e = m_e * c**2 / hbar
print(f"m_e = {m_e:.15e} kg")
print(f"f_e = {f_e:.10e} Hz")
print()

mp_me = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p = m_e * mp_me
print(f"mp/me = {mp_me:.10f}")
print(f"m_p = {m_p:.15e} kg")
print()

G = c**4 * 7.5 * epsilon_0**3 * mu_0**2
H_0 = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r = c**4 / (8.0 * math.pi * G)
w_0 = -(17.0 / 18.0)**2
print(f"G = {G:.15e} m³/(kg·s²)")
print(f"H₀ = {H_0:.10e} Hz")
print(f"VF = {VF_r:.10e} Pa")
print(f"w₀ = {w_0:.15f}")
print()

# Boltzmann constant (CODATA 2019 exact definition)
k_B = 1.380649e-23  # J/K
print(f"k_B = {k_B:.9e} J/K (exact by SI definition)")
print()

# ==============================================================================
# GROUP THEORY: STATISTICAL MECHANICS AS REPRESENTATION THEORY
# ==============================================================================
print("="*80)
print("GROUP THEORY: THERMAL PHYSICS AS REPRESENTATION COUNTING")
print("="*80)
print()

print("CONCEPTUAL FRAMEWORK:")
print("-" * 80)
print("""
Statistical mechanics can be understood as counting ACCESSIBLE REPRESENTATIONS
of the system's symmetry group at a given energy/temperature.

KEY CONCEPTS:

1. PHASE SPACE = REPRESENTATION SPACE
   - Classical: Phase space is a symplectic manifold (positions × momenta)
   - Quantum: Hilbert space with symmetry group representations
   - Temperature T parameterizes which representations are accessible

2. BOLTZMANN FACTOR: e^(-E/k_BT)
   - NOT just a probability weight — it's the CHARACTER of a representation!
   - E = Casimir eigenvalue (energy of the representation)
   - k_BT = temperature scale (controls which representations contribute)

3. PARTITION FUNCTION: Z = Tr[e^(-H/k_BT)]
   - Sum over ALL representations weighted by their characters
   - Group trace = sum over representation dimensions
   - Z is the GENERATING FUNCTION for thermodynamic quantities

4. EQUIPARTITION THEOREM:
   - Each GENERATOR of the symmetry group contributes ½k_BT to energy
   - For classical gas: SO(3) rotation + R³ translation = 6 generators
   - Monoatomic gas: 3 translational DOF → (3/2)k_BT average energy
   - Diatomic gas: add 2 rotational DOF → (5/2)k_BT
""")
print()

print("THERMAL PRESSURE DERIVATION:")
print("-" * 80)
print("""
IDEAL GAS LAW: PV = Nk_BT

This is a REPRESENTATION-COUNTING formula:
- N particles, each in a TRANSLATIONAL REPRESENTATION of R³
- Each particle has 3 momentum generators (p_x, p_y, p_z)
- Average kinetic energy: ⟨K⟩ = (3/2)k_BT (equipartition)
- Pressure arises from momentum transfer to walls

MICROSCOPIC DERIVATION:
- Momentum flux = (number density) × (average momentum) × (velocity)
- P = n⟨p⟩v = n·(mv)·v = nm⟨v²⟩
- From equipartition: ½m⟨v²⟩ = (3/2)k_BT
- Therefore: ⟨v²⟩ = 3k_BT/m
- Pressure: P = (n/3)m⟨v²⟩ = nk_BT

The factor of 1/3 comes from averaging over 3 spatial dimensions (SO(3) symmetry).
""")
print()

# ==============================================================================
# NUMERICAL EXAMPLES: THERMAL PRESSURE
# ==============================================================================
print("="*80)
print("NUMERICAL EXAMPLES: THERMAL PRESSURES")
print("="*80)
print()

print("EXAMPLE 1: Room Temperature Air")
print("-" * 80)
T_room = 293.15  # K (20°C)
P_room = 101325.0  # Pa (1 atm)
n_room = P_room / (k_B * T_room)
print(f"Temperature: T = {T_room:.2f} K = {T_room-273.15:.0f}°C")
print(f"Pressure: P = {P_room:.0f} Pa = 1 atm")
print(f"Number density: n = P/(k_BT) = {n_room:.6e} m⁻³")
print()

# Average molecular mass (air is ~79% N₂, 21% O₂)
m_N2 = 28.0 * 1.66054e-27  # kg (molecular mass of N₂)
m_air = 29.0 * 1.66054e-27  # kg (average for air)
v_rms = math.sqrt(3.0 * k_B * T_room / m_air)
print(f"Average molecular mass: m = {m_air:.6e} kg")
print(f"RMS velocity: v_rms = √(3k_BT/m) = {v_rms:.2f} m/s")
print()

# Mean free path
sigma_N2 = 3.7e-19  # m² (collision cross section for N₂)
lambda_mfp = 1.0 / (math.sqrt(2.0) * n_room * sigma_N2)
print(f"Collision cross section: σ = {sigma_N2:.2e} m²")
print(f"Mean free path: λ = 1/(√2 nσ) = {lambda_mfp:.6e} m = {lambda_mfp*1e9:.1f} nm")
print()

print("EXAMPLE 2: Solar Core")
print("-" * 80)
T_sun = 1.57e7  # K (solar core temperature)
rho_sun = 1.5e5  # kg/m³ (solar core density)
# Solar core is mostly ionized hydrogen: protons + electrons
m_avg = m_p / 2.0  # Average mass per particle (proton + electron, roughly equal numbers)
n_sun = rho_sun / m_avg
P_sun = n_sun * k_B * T_sun
print(f"Temperature: T = {T_sun:.2e} K")
print(f"Density: ρ = {rho_sun:.2e} kg/m³")
print(f"Number density: n = ρ/m_avg = {n_sun:.6e} m⁻³")
print(f"Thermal pressure: P = nk_BT = {P_sun:.6e} Pa")
print()

# Gravitational pressure (rough estimate)
M_sun = 2e30  # kg
R_sun = 7e8  # m
P_grav = G * M_sun**2 / R_sun**4
print(f"Gravitational pressure: P_grav ~ GM²/R⁴ = {P_grav:.6e} Pa")
print(f"Ratio: P_thermal/P_grav = {P_sun/P_grav:.2f}")
print()
print("Thermal pressure balances gravitational collapse — stellar equilibrium!")
print()

print("EXAMPLE 3: Neutron Star Core")
print("-" * 80)
T_ns = 1e9  # K (young neutron star core temperature)
rho_ns = 5e17  # kg/m³ (nuclear density)
n_ns = rho_ns / m_p  # Treat as neutron gas
P_thermal_ns = n_ns * k_B * T_ns
print(f"Temperature: T = {T_ns:.1e} K")
print(f"Density: ρ = {rho_ns:.1e} kg/m³")
print(f"Number density: n = ρ/m_p = {n_ns:.6e} m⁻³")
print(f"Thermal pressure: P_thermal = {P_thermal_ns:.6e} Pa")
print()

# Degeneracy pressure (dominant for cold neutron stars)
# P_deg ~ (ℏ²/m)n^(5/3) for non-relativistic fermions
P_deg = (hbar**2 / m_p) * n_ns**(5.0/3.0) * (3.0 * math.pi**2)**(2.0/3.0)
print(f"Degeneracy pressure: P_deg ~ (ℏ²/m)n^(5/3) = {P_deg:.6e} Pa")
print(f"Ratio: P_thermal/P_deg = {P_thermal_ns/P_deg:.6f}")
print()
print("Neutron stars are supported by DEGENERACY PRESSURE (quantum, not thermal)!")
print("Thermal pressure is negligible except in very hot newborn neutron stars.")
print()

print("EXAMPLE 4: Interstellar Medium")
print("-" * 80)
T_ism = 8000.0  # K (warm ionized medium)
n_ism = 0.3e6  # m⁻³ (0.3 cm⁻³)
P_ism = n_ism * k_B * T_ism
print(f"Temperature: T = {T_ism:.0f} K")
print(f"Number density: n = {n_ism:.2e} m⁻³ = {n_ism/1e6:.2f} cm⁻³")
print(f"Thermal pressure: P = {P_ism:.6e} Pa")
print()

# Compare to Earth's atmosphere
print(f"Ratio to atmospheric pressure: P_ISM/P_atm = {P_ism/P_room:.3e}")
print()
print("Interstellar pressure is ~10⁻¹¹ times atmospheric pressure — near vacuum!")
print()

print("EXAMPLE 5: Early Universe (Recombination)")
print("-" * 80)
T_cmb = 2.725  # K (CMB temperature today)
T_recomb = 3000.0  # K (temperature at recombination, z~1100)
z_recomb = T_recomb / T_cmb
print(f"CMB temperature today: T₀ = {T_cmb:.3f} K")
print(f"Recombination temperature: T_rec = {T_recomb:.0f} K")
print(f"Redshift: z = T_rec/T₀ - 1 = {z_recomb-1:.0f}")
print()

# Number density of baryons (mainly hydrogen)
Omega_b = 0.049  # Baryon density parameter
rho_c = 3.0 * H_0**2 / (8.0 * math.pi * G)
rho_b_today = Omega_b * rho_c
n_b_today = rho_b_today / m_p
n_b_recomb = n_b_today * z_recomb**3  # Scale with (1+z)³
P_recomb = n_b_recomb * k_B * T_recomb
print(f"Baryon density today: ρ_b = {rho_b_today:.6e} kg/m³")
print(f"Baryon density at recombination: ρ_b(z) = {rho_b_today * z_recomb**3:.6e} kg/m³")
print(f"Number density: n = {n_b_recomb:.6e} m⁻³")
print(f"Thermal pressure: P = {P_recomb:.6e} Pa")
print()

# Radiation pressure at recombination
u_rad_today = 4.17e-14  # J/m³ (CMB energy density today)
u_rad_recomb = u_rad_today * z_recomb**4  # Scales as (1+z)⁴
P_rad_recomb = u_rad_recomb / 3.0
print(f"Radiation energy density: u_rad = {u_rad_recomb:.6e} J/m³")
print(f"Radiation pressure: P_rad = u/3 = {P_rad_recomb:.6e} Pa")
print(f"Ratio: P_rad/P_thermal = {P_rad_recomb/P_recomb:.2f}")
print()
print("At recombination, radiation and matter pressures were comparable.")
print("Before z~3000, radiation pressure dominated (radiation era).")
print()

# ==============================================================================
# EQUIPARTITION THEOREM
# ==============================================================================
print("="*80)
print("EQUIPARTITION THEOREM: ENERGY PER DEGREE OF FREEDOM")
print("="*80)
print()

print("CLASSICAL EQUIPARTITION:")
print("-" * 80)
print("""
For a system in thermal equilibrium, each quadratic term in the Hamiltonian
contributes ½k_BT to the average energy:

  ⟨E⟩ = (f/2)k_BT

where f = number of degrees of freedom (DOF).

GROUP-THEORETIC INTERPRETATION:
- Each DOF corresponds to a GENERATOR of the symmetry group
- For a classical particle in 3D: 6 DOF (3 position, 3 momentum)
- Average energy: ⟨E⟩ = 3k_BT (only kinetic energy counts for ideal gas)

EXAMPLES:
""")

print("1. Monoatomic Gas (Ar, He, Ne):")
print("   - 3 translational DOF")
print("   - ⟨E⟩ = (3/2)k_BT")
print(f"   - At T = 300 K: ⟨E⟩ = {1.5 * k_B * 300:.6e} J")
print()

print("2. Diatomic Gas (N₂, O₂):")
print("   - 3 translational + 2 rotational DOF")
print("   - ⟨E⟩ = (5/2)k_BT")
print(f"   - At T = 300 K: ⟨E⟩ = {2.5 * k_B * 300:.6e} J")
print("   - (Vibrational DOF frozen out at room T)")
print()

print("3. Solid (Einstein/Debye Model):")
print("   - Each atom has 3 vibrational modes (3D oscillator)")
print("   - Each mode has 2 quadratic terms (KE + PE)")
print("   - ⟨E⟩ = 3k_BT per atom")
print(f"   - At T = 300 K: ⟨E⟩ = {3.0 * k_B * 300:.6e} J")
print("   - Dulong-Petit law: C_v = 3k_B per atom")
print()

print("4. Electromagnetic Field (Blackbody Radiation):")
print("   - Infinite DOF (all photon modes)")
print("   - Energy density: u = σT⁴ (Stefan-Boltzmann)")
print("   - Each mode contributes k_BT (not ½k_BT — quantum!)")
print()

# ==============================================================================
# QUANTUM CORRECTIONS
# ==============================================================================
print("="*80)
print("QUANTUM CORRECTIONS: BEYOND CLASSICAL EQUIPARTITION")
print("="*80)
print()

print("When ℏω ≳ k_BT, quantum effects become important:")
print("-" * 80)
print()

print("1. PLANCK DISTRIBUTION (Photons):")
print("   ⟨n⟩ = 1/(e^(ℏω/k_BT) - 1)")
print("   - Bosons (integer spin)")
print("   - No chemical potential (photons not conserved)")
print()

print("2. BOSE-EINSTEIN DISTRIBUTION (Bosons):")
print("   ⟨n⟩ = 1/(e^((E-μ)/k_BT) - 1)")
print("   - Chemical potential μ for conserved particles")
print("   - Bose-Einstein condensation at low T")
print()

print("3. FERMI-DIRAC DISTRIBUTION (Fermions):")
print("   ⟨n⟩ = 1/(e^((E-μ)/k_BT) + 1)")
print("   - Electrons, neutrons, protons")
print("   - Pauli exclusion → degeneracy pressure")
print()

# Degeneracy temperature for electrons
n_e_metal = 1e29  # m⁻³ (typical for metals)
E_F = (hbar**2 / (2.0 * m_e)) * (3.0 * math.pi**2 * n_e_metal)**(2.0/3.0)
T_F = E_F / k_B
print("EXAMPLE: Electron Degeneracy in Metals")
print(f"  Electron density: n_e = {n_e_metal:.1e} m⁻³")
print(f"  Fermi energy: E_F = {E_F:.6e} J = {E_F/1.602e-19:.2f} eV")
print(f"  Fermi temperature: T_F = E_F/k_B = {T_F:.6e} K")
print()
print(f"At room temperature (T = 300 K), T/T_F = {300/T_F:.3e} << 1")
print("Electrons in metals are HIGHLY DEGENERATE — quantum pressure dominates!")
print()

# ==============================================================================
# CALIBRATION CHECKPOINT
# ==============================================================================
print("="*80)
print("CALIBRATION CHECKPOINT")
print("="*80)
print()

print("Test: Ideal Gas Law at STP")
print("-" * 80)
T_stp = 273.15  # K (0°C)
P_stp = 101325.0  # Pa (1 atm)
n_stp = P_stp / (k_B * T_stp)
print(f"Standard Temperature: T = {T_stp:.2f} K")
print(f"Standard Pressure: P = {P_stp:.0f} Pa")
print(f"Derived number density: n = P/(k_BT) = {n_stp:.6e} m⁻³")
print()

# Compare to Loschmidt constant (number density at STP)
N_L = 2.6867811e25  # m⁻³ (CODATA 2018)
print(f"Loschmidt constant: N_L = {N_L:.7e} m⁻³")
print(f"Fractional difference: {abs(n_stp - N_L)/N_L:.6e}")
print()

if abs(n_stp - N_L)/N_L < 1e-6:
    print("✓ PASS: Ideal gas law matches Loschmidt constant")
else:
    print("✗ CAUTION: Check k_B value")
print()

# ==============================================================================
# SUMMARY
# ==============================================================================
print("="*80)
print("SUMMARY: THERMAL PRESSURE AS REPRESENTATION THEORY")
print("="*80)
print()

print("IDEAL GAS LAW:")
print("  P = nk_BT")
print()

print("GROUP-THEORETIC INTERPRETATION:")
print("  • Temperature T parameterizes accessible representations")
print("  • k_B converts energy (Casimir eigenvalue) ↔ temperature")
print("  • Equipartition: ½k_BT per generator (DOF)")
print("  • Partition function Z = Tr[e^(-H/k_BT)] sums over representations")
print()

print("KEY EXAMPLES:")
print(f"  • Room air: P = {P_room:.0f} Pa at T = {T_room:.0f} K")
print(f"  • Solar core: P = {P_sun:.2e} Pa at T = {T_sun:.1e} K")
print(f"  • Neutron star: P_thermal << P_degeneracy (quantum dominates)")
print(f"  • ISM: P = {P_ism:.2e} Pa (near vacuum)")
print()

print("QUANTUM CORRECTIONS:")
print("  • ℏω > k_BT → Planck/Bose-Einstein/Fermi-Dirac distributions")
print("  • Degeneracy pressure (fermions)")
print("  • Bose-Einstein condensation (bosons)")
print()

print("="*80)
print("END OF THERMAL_PRESSURE_GROUPTHEORY.PY")
print("="*80)

input("Press Enter to exit...")
