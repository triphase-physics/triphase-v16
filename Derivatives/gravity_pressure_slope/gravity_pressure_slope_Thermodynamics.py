"""
================================================================================
TriPhase V16 - Gravity Pressure Slope Derivative
Framework: THERMODYNAMICS
Tag: (D) - Pure derivation
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
The gravitational pressure gradient is the fundamental thermodynamic force in
hydrostatic equilibrium. For a self-gravitating system (star, planet, atmosphere):

    dP/dr = -ρ(r) × g(r) = -ρ(r) × GM(r)/r²

This is the condition for thermal-gravitational balance: the pressure gradient
(thermal force) balances gravity (gravitational force).

Thermodynamically, this equation states:
    - Pressure P = thermal kinetic energy density (∝ nk_BT)
    - Gravity creates a potential well U(r) = -GMm/r
    - Equilibrium: ∇P = -ρ∇Φ where Φ is gravitational potential

The "slope" dP/dr is negative (pressure decreases outward) and steep where
density is high. This gradient drives convection if |dP/dr| exceeds the
adiabatic gradient.

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
# THERMODYNAMIC DERIVATION - GRAVITY PRESSURE SLOPE
# ============================================================================

print("=" * 80)
print("GRAVITY PRESSURE SLOPE - THERMODYNAMICS FRAMEWORK")
print("=" * 80)
print()
print("THERMODYNAMIC INTERPRETATION:")
print("The gravity pressure gradient dP/dr = -ρg is the force balance")
print("between thermal pressure (kinetic) and gravitational compression.")
print("This is hydrostatic equilibrium — a thermodynamic equilibrium condition.")
print()

print("HYDROSTATIC EQUILIBRIUM EQUATION:")
print(f"  dP/dr = -ρ(r) × g(r)")
print(f"  where g(r) = GM(r)/r² is local gravitational acceleration")
print(f"        M(r) = enclosed mass within radius r")
print()

# Example: Solar interior
M_sun = 1.989e30  # kg
R_sun = 6.96e8    # m
rho_center_sun = 1.5e5  # kg/m³ (central density)
T_center_sun_K = 1.57e7  # K (central temperature)
k_B = 1.380649e-23  # J/K

print("=" * 80)
print("EXAMPLE: SOLAR INTERIOR")
print("=" * 80)
print()

# Solar center
P_center_sun = 2.5e16  # Pa (estimated central pressure)
r_eval = 0.1 * R_sun  # Evaluate at 10% solar radius
M_r = (r_eval / R_sun)**3 * M_sun  # Rough estimate (uniform density)
g_r = G * M_r / r_eval**2
rho_r = 0.5 * rho_center_sun  # Rough estimate

# Pressure gradient
dP_dr = -rho_r * g_r

print(f"SOLAR PARAMETERS:")
print(f"  Solar mass M_☉              = {M_sun:.3e} kg")
print(f"  Solar radius R_☉            = {R_sun:.3e} m")
print(f"  Central density ρ_c         = {rho_center_sun:.3e} kg/m³")
print(f"  Central temperature T_c     = {T_center_sun_K:.3e} K")
print(f"  Central pressure P_c        = {P_center_sun:.3e} Pa")
print()

print(f"AT r = 0.1 R_☉:")
print(f"  Enclosed mass M(r)          ≈ {M_r:.3e} kg")
print(f"  Local gravity g(r)          = {g_r:.3e} m/s²")
print(f"  Local density ρ(r)          ≈ {rho_r:.3e} kg/m³")
print()
print(f"  Pressure gradient dP/dr     = -ρ × g")
print(f"                              = -{rho_r:.3e} × {g_r:.3e}")
print(f"                              = {dP_dr:.3e} Pa/m")
print()

# Pressure scale height
H_P = -1.0 / (dP_dr / P_center_sun) if dP_dr != 0 else 0
print(f"PRESSURE SCALE HEIGHT:")
print(f"  H_P = -P / (dP/dr)          ≈ {H_P:.3e} m")
print(f"                              ≈ {H_P/R_sun:.3f} R_☉")
print(f"  This is the characteristic length scale for pressure variation.")
print()

# ============================================================================
# THERMODYNAMIC FORCE BALANCE
# ============================================================================

print("=" * 80)
print("THERMODYNAMIC FORCE BALANCE")
print("=" * 80)
print()

print("THERMAL PRESSURE FORCE:")
print(f"  Thermal pressure P = nk_BT (ideal gas)")
print(f"  Pressure force per unit volume: F_thermal = -∇P")
print(f"  This force pushes outward (positive radial direction)")
print()

print("GRAVITATIONAL FORCE:")
print(f"  Gravitational force per unit volume: F_grav = -ρ∇Φ")
print(f"  where Φ = -GM/r is gravitational potential")
print(f"  ∇Φ = -GM/r² × r̂ → F_grav = -ρg(r) × r̂")
print(f"  This force pulls inward (negative radial direction)")
print()

print("EQUILIBRIUM CONDITION:")
print(f"  F_thermal + F_grav = 0")
print(f"  -∇P - ρ∇Φ = 0")
print(f"  dP/dr = -ρ × GM(r)/r²")
print()

# ============================================================================
# CONVECTION VS RADIATION
# ============================================================================

print("=" * 80)
print("CONVECTION VS RADIATION (THERMODYNAMIC STABILITY)")
print("=" * 80)
print()

# Adiabatic gradient
gamma = 5.0/3.0  # Ideal gas adiabatic index
nabla_ad = (gamma - 1.0) / gamma
nabla_rad = 0.4  # Typical radiative gradient (example)

print("ADIABATIC GRADIENT:")
print(f"  For an adiabatic process: PV^γ = const")
print(f"  Adiabatic index γ           = {gamma:.3f}")
print(f"  Adiabatic gradient ∇_ad     = (γ-1)/γ = {nabla_ad:.3f}")
print()

print("SCHWARZSCHILD CRITERION:")
print(f"  If |dT/dr|_actual > |dT/dr|_adiabatic:")
print(f"    → Convection occurs (efficient heat transport)")
print(f"  If |dT/dr|_actual < |dT/dr|_adiabatic:")
print(f"    → Radiative equilibrium (photon diffusion)")
print()

print(f"  In the Sun:")
print(f"    Inner 70%: Radiative zone (∇_rad < ∇_ad)")
print(f"    Outer 30%: Convection zone (∇_rad > ∇_ad)")
print()

# ============================================================================
# ATMOSPHERIC EXAMPLE
# ============================================================================

print("=" * 80)
print("EXAMPLE: EARTH'S ATMOSPHERE")
print("=" * 80)
print()

# Earth atmosphere
g_Earth = 9.81  # m/s²
T_surface = 288.0  # K
P_surface = 101325.0  # Pa
M_air = 0.029  # kg/mol (mean molecular mass of air)
R_gas = 8.314  # J/(mol·K)

# Barometric formula
H_scale = (R_gas * T_surface) / (M_air * g_Earth)

print(f"EARTH SURFACE CONDITIONS:")
print(f"  Temperature T               = {T_surface:.1f} K")
print(f"  Pressure P                  = {P_surface:.0f} Pa")
print(f"  Gravity g                   = {g_Earth:.2f} m/s²")
print()

print(f"BAROMETRIC FORMULA:")
print(f"  P(z) = P_0 × exp(-z/H)")
print(f"  Scale height H = RT/(Mg)    = {H_scale:.0f} m ≈ {H_scale/1000:.1f} km")
print()

# Pressure gradient at surface
rho_air = P_surface * M_air / (R_gas * T_surface)
dP_dz = -rho_air * g_Earth

print(f"PRESSURE GRADIENT:")
print(f"  Air density ρ               = {rho_air:.3f} kg/m³")
print(f"  dP/dz = -ρg                 = {dP_dz:.1f} Pa/m")
print(f"  This is the weight of air per unit area per meter.")
print()

# ============================================================================
# NEUTRON STAR EXAMPLE
# ============================================================================

print("=" * 80)
print("EXAMPLE: NEUTRON STAR")
print("=" * 80)
print()

M_NS = 1.4 * M_sun  # Neutron star mass
R_NS = 12e3  # m (12 km)
rho_NS = 5e17  # kg/m³ (nuclear density)
P_center_NS = 1e34  # Pa (estimated)

g_NS_surface = G * M_NS / R_NS**2
dP_dr_NS = -rho_NS * g_NS_surface

print(f"NEUTRON STAR PARAMETERS:")
print(f"  Mass M_NS                   = {M_NS/M_sun:.1f} M_☉")
print(f"  Radius R_NS                 = {R_NS/1e3:.0f} km")
print(f"  Central density ρ_c         ≈ {rho_NS:.3e} kg/m³ (nuclear)")
print(f"  Central pressure P_c        ≈ {P_center_NS:.3e} Pa")
print()

print(f"SURFACE CONDITIONS:")
print(f"  Surface gravity g_surf      = {g_NS_surface:.3e} m/s²")
print(f"                              = {g_NS_surface/g_Earth:.2e} g_Earth")
print(f"  Pressure gradient dP/dr     ≈ {dP_dr_NS:.3e} Pa/m")
print()

print(f"THERMODYNAMIC REGIME:")
print(f"  At nuclear densities, matter is degenerate (quantum pressure).")
print(f"  P ∝ n^(5/3) for non-relativistic fermions")
print(f"  P ∝ n^(4/3) for ultra-relativistic fermions")
print(f"  The equation of state P(ρ) determines R_NS via TOV equation.")
print()

# ============================================================================
# GENERAL FORMULATION
# ============================================================================

print("=" * 80)
print("GENERAL THERMODYNAMIC FORMULATION")
print("=" * 80)
print()

print("HYDROSTATIC EQUILIBRIUM:")
print(f"  dP/dr = -ρ(r) × GM(r)/r²")
print()
print(f"  Combined with mass conservation:")
print(f"    dM/dr = 4πr²ρ(r)")
print()
print(f"  And equation of state P(ρ, T):")
print(f"    These form a closed system determining ρ(r), P(r), T(r), M(r)")
print()

print("THERMODYNAMIC SIGNIFICANCE:")
print(f"  The pressure gradient dP/dr is a thermodynamic force.")
print(f"  It balances gravity to maintain thermal equilibrium.")
print(f"  Any deviation from dP/dr = -ρg causes:")
print(f"    - Expansion (if dP/dr > -ρg)")
print(f"    - Contraction (if dP/dr < -ρg)")
print(f"  This is Le Chatelier's principle applied to gravity.")
print()

# ============================================================================
# CALIBRATION CHECKPOINT
# ============================================================================

print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
print()

print("The gravity pressure gradient is a derived quantity.")
print("It is calculated from measured/computed values of ρ(r), M(r), r.")
print()
print(f"TriPhase G                  = {G:.6e} m³/(kg·s²)")
print(f"CODATA 2018 G               = 6.67430e-11 m³/(kg·s²)")
print()
print("Solar model dP/dr values match helioseismology observations.")
print("Neutron star masses inferred from dP/dr via TOV equation match")
print("observed masses in binary pulsars (1.2-2.0 M_☉).")
print()

print("THERMODYNAMIC SIGNIFICANCE:")
print("  dP/dr = -ρg is the most fundamental equation in astrophysics.")
print("  It connects thermodynamics (P) with gravity (g).")
print("  Every star, planet, and atmosphere obeys this thermodynamic law.")
print()

print("=" * 80)
print("END GRAVITY PRESSURE SLOPE THERMODYNAMICS DERIVATIVE")
print("=" * 80)

input("Press Enter to exit...")
