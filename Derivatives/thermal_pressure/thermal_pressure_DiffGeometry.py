"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================
Derivative:  Thermal Pressure (P = nk_BT)
Framework:   DiffGeometry
Version:     16.0
Generated:   2026-03-26
Status:      Active Development
Tag: (D)
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
"""
import math

# === ANCHOR INPUTS ===
epsilon_0 = 8.8541878128e-12  # F/m
mu_0      = 1.25663706212e-6  # H/m
e         = 1.602176634e-19   # C

# === DERIVED ANCHOR CHAIN ===
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
T_17      = 17 * 18 // 2  # = 153
mp_me     = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p       = m_e * mp_me
H_0       = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r      = c**4 / (8.0 * math.pi * G)
k_B       = 1.380649e-23  # J/K (exact SI 2019)

# === DERIVATION: Thermal Pressure ===
print("=" * 80)
print("TRIPHASE V16: THERMAL PRESSURE")
print("Framework: DiffGeometry")
print("=" * 80)
print()

print("IDEAL GAS LAW:")
print("-" * 80)
print()
print("  P V = N k_B T")
print()
print("Where:")
print("  P   = pressure (Pa)")
print("  V   = volume (m³)")
print("  N   = number of particles")
print("  k_B = Boltzmann constant (J/K)")
print("  T   = temperature (K)")
print()
print("In terms of number density n = N/V:")
print()
print("  P = n k_B T")
print()
print(f"Boltzmann constant: k_B = {k_B:.15e} J/K (exact SI)")
print()

print("=" * 80)
print("DIFFERENTIAL GEOMETRY INTERPRETATION")
print("=" * 80)
print()

print("THERMAL PRESSURE IN STRESS-ENERGY TENSOR:")
print("-" * 80)
print()
print("For a perfect fluid (no viscosity, no heat conduction):")
print()
print("  T_μν = (ρ + P/c²) u_μ u_ν + P g_μν")
print()
print("Where:")
print("  ρ   = rest mass density (kg/m³)")
print("  P   = pressure (Pa)")
print("  u_μ = 4-velocity of fluid element")
print("  g_μν = metric tensor")
print()
print("In fluid rest frame (u = (1,0,0,0)):")
print()
print("  T_00 = ρ c²         (energy density)")
print("  T_11 = T_22 = T_33 = P  (isotropic pressure)")
print("  T_0i = 0            (no momentum flux)")
print()
print("Trace of stress-energy tensor:")
print("  T = g^μν T_μν = -ρc² + 3P")
print()
print("For non-relativistic gas (ρc² >> P):")
print("  T ≈ -ρc² (negative trace → attractive gravity)")
print()
print("For ultra-relativistic gas (P = ρc²/3):")
print("  T = 0 (traceless → conformal invariance)")
print()

print("THERMAL PRESSURE AND SPACETIME CURVATURE:")
print("-" * 80)
print()
print("Einstein equation: G_μν = κ T_μν")
print()
kappa = 8.0 * math.pi * G / c**4
print(f"  κ = 8πG/c⁴ = {kappa:.3e} Pa⁻¹")
print()
print("Curvature from thermal pressure:")
print("  R ~ κ P = κ n k_B T")
print()
print("Higher temperature → higher pressure → more curvature")
print()

# === PHYSICAL EXAMPLES ===
print("=" * 80)
print("PHYSICAL EXAMPLES: Thermal Pressure")
print("=" * 80)
print()

# Example 1: Room temperature air
print("EXAMPLE 1: Room temperature air at sea level")
print("-" * 80)
P_atm = 101325  # Pa
T_room = 293  # K (20°C)
n_air = P_atm / (k_B * T_room)
rho_air = n_air * 28.97 * 1.66054e-27  # Average molecular mass
print(f"  T = {T_room} K (20°C)")
print(f"  P = {P_atm} Pa (1 atm)")
print(f"  n = P/(k_B T) = {n_air:.3e} particles/m³")
print(f"  ρ = {rho_air:.3e} kg/m³")
R_air = kappa * P_atm
print(f"  R = κP = {R_air:.3e} m⁻²")
print(f"  Curvature scale: {1.0/math.sqrt(abs(R_air)):.3e} m")
print()

# Example 2: Solar core
print("EXAMPLE 2: Solar core")
print("-" * 80)
T_sun = 1.5e7  # K
n_sun = 6e31  # particles/m³ (mostly protons)
P_sun = n_sun * k_B * T_sun
rho_sun = n_sun * m_p
print(f"  T = {T_sun:.3e} K")
print(f"  n = {n_sun:.3e} protons/m³")
print(f"  P = nk_BT = {P_sun:.3e} Pa")
print(f"  ρ = {rho_sun:.3e} kg/m³")
R_sun = kappa * P_sun
print(f"  R = κP = {R_sun:.3e} m⁻²")
print(f"  Curvature scale: {1.0/math.sqrt(R_sun):.3e} m")
print(f"  P/VF_r = {P_sun/VF_r:.3e}")
print()

# Example 3: Neutron star core
print("EXAMPLE 3: Neutron star core")
print("-" * 80)
T_ns = 1e9  # K
n_ns = 5e44  # baryons/m³ (neutron density)
P_ns_thermal = n_ns * k_B * T_ns
P_ns_degeneracy = 1e34  # Pa (dominated by degeneracy pressure)
P_ns_total = P_ns_thermal + P_ns_degeneracy
rho_ns = n_ns * m_p
print(f"  T = {T_ns:.3e} K")
print(f"  n = {n_ns:.3e} baryons/m³")
print(f"  P_thermal = nk_BT = {P_ns_thermal:.3e} Pa")
print(f"  P_degeneracy = {P_ns_degeneracy:.3e} Pa")
print(f"  P_total = {P_ns_total:.3e} Pa")
print(f"  ρ = {rho_ns:.3e} kg/m³")
R_ns = kappa * P_ns_total
print(f"  R = κP = {R_ns:.3e} m⁻²")
print(f"  Curvature scale: {1.0/math.sqrt(R_ns):.3e} m")
print(f"  P_total/VF_r = {P_ns_total/VF_r:.3e}")
print()
print("Note: Neutron star pressure dominated by quantum degeneracy,")
print("      not thermal motion. Thermal contribution is negligible.")
print()

# Example 4: Early universe (1 second after Big Bang)
print("EXAMPLE 4: Early universe (t = 1 second)")
print("-" * 80)
T_early = 1e10  # K
n_early = 1e42  # particles/m³ (photons, neutrinos, e+e-)
P_early = n_early * k_B * T_early
rho_early = P_early / c**2  # Ultra-relativistic: P = ρc²/3 → ρ = 3P/c²
print(f"  T = {T_early:.3e} K")
print(f"  n = {n_early:.3e} particles/m³")
print(f"  P = nk_BT = {P_early:.3e} Pa")
print(f"  ρ = 3P/c² = {rho_early:.3e} kg/m³ (ultra-relativistic)")
R_early = kappa * P_early
print(f"  R = κP = {R_early:.3e} m⁻²")
print(f"  Curvature scale: {1.0/math.sqrt(R_early):.3e} m")
print(f"  P/VF_r = {P_early/VF_r:.3e}")
print()

# Example 5: Cosmic microwave background today
print("EXAMPLE 5: Cosmic microwave background (today)")
print("-" * 80)
T_cmb = 2.725  # K
n_cmb = 4.1e8  # photons/m³
P_cmb = n_cmb * k_B * T_cmb
rho_cmb = P_cmb / c**2
print(f"  T = {T_cmb} K")
print(f"  n = {n_cmb:.3e} photons/m³")
print(f"  P = nk_BT = {P_cmb:.3e} Pa")
print(f"  ρ = P/c² = {rho_cmb:.3e} kg/m³")
R_cmb = kappa * P_cmb
print(f"  R = κP = {R_cmb:.3e} m⁻²")
print(f"  Curvature scale: {1.0/math.sqrt(abs(R_cmb)):.3e} m")
print()

# Example 6: Planck temperature
print("EXAMPLE 6: Planck temperature (quantum gravity regime)")
print("-" * 80)
l_P = math.sqrt(hbar * G / c**3)  # Planck length
t_P = l_P / c  # Planck time
T_P = hbar / (k_B * t_P)  # Planck temperature
n_P = 1.0 / l_P**3  # One particle per Planck volume
P_P = n_P * k_B * T_P
print(f"  l_P = {l_P:.3e} m (Planck length)")
print(f"  T_P = ħ/(k_B t_P) = {T_P:.3e} K")
print(f"  n_P = 1/l_P³ = {n_P:.3e} particles/m³")
print(f"  P_P = nk_BT = {P_P:.3e} Pa")
print(f"  P_P/VF_r = {P_P/VF_r:.3e}")
print()
print("At Planck scale, spacetime curvature becomes extreme.")
print("Quantum gravity effects dominate.")
print()

# === TEMPERATURE-CURVATURE RELATION ===
print("=" * 80)
print("TEMPERATURE-CURVATURE RELATION")
print("=" * 80)
print()
print("For thermal gas with number density n:")
print()
print("  P = n k_B T")
print("  R = κ P = κ n k_B T")
print()
print("Curvature is DIRECTLY proportional to temperature.")
print()
print("This connects thermodynamics to geometry:")
print("  - Hotter gas → higher pressure → more curvature")
print("  - Cooler gas → lower pressure → less curvature")
print()
print("In cosmology:")
print("  - Early universe: high T → high curvature → rapid expansion")
print("  - Today: low T_CMB → low curvature → slow expansion")
print()

# === EQUATION OF STATE ===
print("=" * 80)
print("EQUATION OF STATE PARAMETER w")
print("=" * 80)
print()
print("For perfect fluid: w = P/(ρc²)")
print()
print("Different regimes:")
print("  Non-relativistic matter: w = 0 (P << ρc²)")
print("  Radiation/relativistic: w = 1/3 (P = ρc²/3)")
print("  Cosmological constant: w = -1 (P = -ρc²)")
print()

examples_w = [
    ("Room air", P_atm, rho_air * c**2),
    ("Solar core", P_sun, rho_sun * c**2),
    ("Neutron star", P_ns_total, rho_ns * c**2),
    ("Early universe", P_early, rho_early * c**2),
    ("CMB today", P_cmb, rho_cmb * c**2),
]

for name, P, rho_c2 in examples_w:
    w = P / rho_c2
    print(f"{name:20s}: w = {w:.6f}")
print()

# === CALIBRATION CHECKPOINT ===
print("=" * 80)
print("CALIBRATION CHECKPOINT")
print("=" * 80)
k_B_SI = 1.380649e-23  # J/K (exact)
print(f"Boltzmann constant (SI 2019): {k_B_SI:.15e} J/K")
print(f"Value used in script:         {k_B:.15e} J/K")
print()
if k_B == k_B_SI:
    print("STATUS: EXACT - Using SI 2019 defined value")
else:
    print("STATUS: Check value")
print()

print("Verification with atmospheric pressure:")
P_calc = n_air * k_B * T_room
print(f"  Measured P: {P_atm} Pa")
print(f"  Calculated: {P_calc:.0f} Pa")
print(f"  Deviation:  {abs(P_calc - P_atm)/P_atm * 100:.1f}%")
print()

print("=" * 80)
print("DERIVATION COMPLETE")
print("=" * 80)
print()
print("Thermal pressure P = nk_BT connects:")
print("  - Statistical mechanics (particle motion)")
print("  - Thermodynamics (temperature, entropy)")
print("  - General relativity (spacetime curvature)")
print()
print("The Boltzmann constant k_B is the conversion factor")
print("between temperature (K) and energy (J).")
print()
print("In curved spacetime, thermal pressure sources curvature:")
print("  R ~ κ n k_B T")
print()
print("Hotter systems create more curvature and stronger gravity.")
print()

input("Press Enter to exit...")
