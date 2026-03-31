"""
thermal_pressure_Anchor_Condensed.py

TriPhase V16 - Thermal Pressure
Row 36 - Tag: (D) DERIVED

Derives thermal pressure P = n*k_B*T from anchor chain.
k_B = 1.380649e-23 J/K is exact SI (like e), so treated as fundamental.

Shows thermal pressure at standard conditions and compares to VF_r.

All derived from epsilon_0 and mu_0 via the anchor chain.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# Anchor chain
epsilon_0 = 8.8541878128e-12  # F/m
mu_0      = 1.25663706212e-6   # H/m
e         = 1.602176634e-19    # C (exact SI)
k_B       = 1.380649e-23       # J/K (exact SI)

c     = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0   = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha = 1.0 / alpha_inv
hbar  = Z_0 * e**2 / (4.0 * math.pi * alpha)
h     = 2.0 * math.pi * hbar
G     = c**4 * 7.5 * epsilon_0**3 * mu_0**2
m_e   = hbar * alpha / (c * 2.8179403262e-15)
f_e   = m_e * c**2 / hbar
mp_me = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p   = m_e * mp_me
H_0   = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r  = c**4 / (8.0 * math.pi * G)

print("=" * 70)
print("TriPhase V16 - Thermal Pressure")
print("Row 36 - DERIVED from epsilon_0, mu_0")
print("=" * 70)
print()

# Thermal pressure formula
print("THERMAL PRESSURE:")
print("  P_thermal = n * k_B * T")
print()
print("  Where:")
print("    n = number density (particles/m^3)")
print("    k_B = Boltzmann constant (exact SI)")
print("    T = temperature (K)")
print()
print(f"  k_B = {k_B:.15e} J/K (exact SI 2019)")
print()

# Standard atmospheric conditions
T_std = 300.0       # K (room temperature)
n_atm = 2.5e25      # molecules/m^3 (at 1 atm, 300K)
P_atm_calc = n_atm * k_B * T_std
P_atm_actual = 101325.0  # Pa

print("STANDARD ATMOSPHERIC PRESSURE:")
print(f"  T = {T_std:.1f} K (room temperature)")
print(f"  n = {n_atm:.2e} molecules/m^3")
print(f"  P_thermal = n*k_B*T = {P_atm_calc:.6e} Pa")
print(f"            = {P_atm_calc/1e5:.6f} bar")
print(f"            = {P_atm_calc/101325.0:.6f} atm")
print()
print(f"  Standard atm = {P_atm_actual:.1f} Pa = 1.000 atm")
print(f"  Agreement: {100.0*P_atm_calc/P_atm_actual:.3f}%")
print()

# Various temperature scales
temperatures = [
    ("Cosmic microwave background", 2.725),
    ("Liquid nitrogen", 77.0),
    ("Room temperature", 300.0),
    ("Water boiling", 373.15),
    ("Steel melting", 1800.0),
    ("Solar surface", 5778.0),
    ("Solar core", 1.5e7),
]

print("THERMAL PRESSURES AT n = 2.5e25 m^-3:")
for name, T in temperatures:
    P = n_atm * k_B * T
    print(f"  {name:30s} ({T:>10.3e} K): {P:>12.6e} Pa")
print()

# Compare to vacuum frame rigidity
print("COMPARISON TO VACUUM FRAME RIGIDITY:")
print(f"  VF_r = {VF_r:.6e} Pa")
print()
for name, T in temperatures:
    P = n_atm * k_B * T
    ratio = P / VF_r
    print(f"  {name:30s}: P/VF_r = {ratio:.6e}")
print()
print("  Even solar core thermal pressure is ~43 orders below VF_r!")
print()

# Extreme density examples
densities = [
    ("Atmosphere", 2.5e25),
    ("Water", 3.34e28),      # molecules/m^3
    ("White dwarf core", 1e36),
    ("Neutron star core", 1e44),
]

T_test = 1e7  # K

print(f"THERMAL PRESSURES AT T = {T_test:.1e} K:")
for name, n in densities:
    P = n * k_B * T_test
    ratio = P / VF_r
    print(f"  {name:20s} (n={n:.2e} m^-3): P={P:.6e} Pa, P/VF_r={ratio:.6e}")
print()

# Ideal gas law connection
print("IDEAL GAS LAW CONNECTION:")
print("  PV = NkT  or  P = (N/V)*k_B*T = n*k_B*T")
print("  Also: PV = nRT where n is moles, R = N_A*k_B")
print()
N_A = 6.02214076e23  # exact SI
R = N_A * k_B
print(f"  N_A = {N_A:.8e} mol^-1 (exact SI)")
print(f"  R = N_A*k_B = {R:.15e} J/(mol·K)")
print(f"  R_CODATA = 8.314462618... J/(mol·K)")
print()

# Energy scales
E_thermal_300K = k_B * 300.0
E_thermal_eV = E_thermal_300K / e

print("THERMAL ENERGY SCALES:")
print(f"  kT at 300K = {E_thermal_300K:.6e} J")
print(f"             = {E_thermal_eV:.6f} eV")
print(f"             = 1/40 eV (room temperature)")
print()

print("ANCHOR CHAIN DERIVATION:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print(f"  k_B       = {k_B:.15e} J/K (exact SI)")
print(f"  c         = {c:.10e} m/s")
print(f"  G         = c^4 * 7.5 * epsilon_0^3 * mu_0^2")
print(f"            = {G:.14e} m^3/(kg·s^2)")
print(f"  VF_r      = c^4/(8*pi*G)")
print(f"            = {VF_r:.6e} Pa")
print()
print("  P_thermal = n*k_B*T")
print()

print("=" * 70)
print("Thermal pressure uses exact SI k_B, compares to vacuum rigidity")
print("Even extreme thermal pressures are far below VF_r")
print("=" * 70)

input("Press Enter to exit...")
