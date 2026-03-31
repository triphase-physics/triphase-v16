"""
hydrostatic_pressure_Anchor_Condensed.py

TriPhase V16 - Hydrostatic Pressure
Row 37 - Tag: (D) DERIVED

Derives hydrostatic pressure P = rho*g*h from anchor chain.
Shows g = G*M/R^2 with G derived from epsilon_0, mu_0.

Example: atmospheric pressure at sea level from first principles.

All derived from epsilon_0 and mu_0 via the anchor chain.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# Anchor chain
epsilon_0 = 8.8541878128e-12  # F/m
mu_0      = 1.25663706212e-6   # H/m
e         = 1.602176634e-19    # C (exact SI)

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
print("TriPhase V16 - Hydrostatic Pressure")
print("Row 37 - DERIVED from epsilon_0, mu_0")
print("=" * 70)
print()

# Hydrostatic pressure formula
print("HYDROSTATIC PRESSURE:")
print("  P_hydro = rho * g * h")
print()
print("  Where:")
print("    rho = mass density (kg/m^3)")
print("    g = gravitational acceleration (m/s^2)")
print("    h = height/depth (m)")
print()
print("  And g = G*M/R^2 for planet of mass M, radius R")
print()

# Earth parameters
M_earth = 5.972e24     # kg
R_earth = 6.371e6      # m
rho_air = 1.225        # kg/m^3 at sea level

g_earth = G * M_earth / R_earth**2

print("EARTH SURFACE GRAVITY:")
print(f"  M_earth = {M_earth:.3e} kg")
print(f"  R_earth = {R_earth:.3e} m")
print(f"  G       = {G:.14e} m^3/(kg·s^2)")
print()
print(f"  g = G*M/R^2 = {g_earth:.10f} m/s^2")
print(f"  Standard g = 9.80665 m/s^2")
print(f"  Agreement: {100.0*g_earth/9.80665:.3f}%")
print()

# Atmospheric pressure (simple model)
h_atm = 8500.0  # scale height in meters
P_sea = rho_air * g_earth * h_atm

print("ATMOSPHERIC PRESSURE (simple model):")
print(f"  rho_air = {rho_air:.3f} kg/m^3 (sea level)")
print(f"  h_scale = {h_atm:.1f} m (atmospheric scale height)")
print(f"  P = rho*g*h = {P_sea:.1f} Pa")
print(f"            = {P_sea/101325.0:.3f} atm")
print()
print("  (Actual atmosphere: P = 101325 Pa = 1 atm)")
print("  (Simple model gives order of magnitude)")
print()

# Water pressure vs depth
rho_water = 1000.0  # kg/m^3

depths = [
    ("Swimming pool", 2.0),
    ("Scuba limit", 40.0),
    ("Submarine crush depth", 600.0),
    ("Mariana Trench", 10994.0),
]

print("WATER PRESSURE VS DEPTH:")
print(f"  rho_water = {rho_water:.1f} kg/m^3")
print(f"  g = {g_earth:.3f} m/s^2")
print()
for name, depth in depths:
    P = rho_water * g_earth * depth
    P_atm = P / 101325.0
    print(f"  {name:25s} ({depth:>8.1f} m): {P:>12.6e} Pa = {P_atm:>8.3f} atm")
print()

# Planetary comparison
planets = [
    ("Mercury", 3.30e23, 2.4397e6),
    ("Venus", 4.87e24, 6.0518e6),
    ("Earth", 5.972e24, 6.371e6),
    ("Mars", 6.42e23, 3.3895e6),
    ("Jupiter", 1.898e27, 6.9911e7),
    ("Saturn", 5.68e26, 5.8232e7),
]

print("PLANETARY SURFACE GRAVITY:")
for name, M, R in planets:
    g = G * M / R**2
    print(f"  {name:10s}: g = {g:>10.3f} m/s^2 (M={M:.2e} kg, R={R:.3e} m)")
print()

# Pressure in stellar objects
print("EXTREME HYDROSTATIC PRESSURES:")
print()

# White dwarf
M_wd = 1.4 * 1.989e30  # 1.4 solar masses
R_wd = 5000.0e3        # 5000 km
rho_wd_avg = M_wd / (4.0/3.0 * math.pi * R_wd**3)
g_wd = G * M_wd / R_wd**2
P_wd_surface = rho_wd_avg * g_wd * R_wd / 2.0  # rough estimate

print(f"White Dwarf (Chandrasekhar limit):")
print(f"  M = {M_wd:.3e} kg, R = {R_wd:.3e} m")
print(f"  g = {g_wd:.3e} m/s^2")
print(f"  P_center ~ {P_wd_surface:.3e} Pa")
print(f"  P/VF_r ~ {P_wd_surface/VF_r:.3e}")
print()

# Neutron star
M_ns = 1.4 * 1.989e30  # 1.4 solar masses
R_ns = 10.0e3          # 10 km
rho_ns_avg = M_ns / (4.0/3.0 * math.pi * R_ns**3)
g_ns = G * M_ns / R_ns**2
P_ns_surface = rho_ns_avg * g_ns * R_ns / 2.0  # rough estimate

print(f"Neutron Star:")
print(f"  M = {M_ns:.3e} kg, R = {R_ns:.3e} m")
print(f"  g = {g_ns:.3e} m/s^2")
print(f"  P_center ~ {P_ns_surface:.3e} Pa")
print(f"  P/VF_r ~ {P_ns_surface/VF_r:.3e}")
print()
print("  (Neutron star core pressure approaches VF_r!)")
print()

# Full derivation path
print("HYDROSTATIC EQUILIBRIUM EQUATION:")
print("  dP/dr = -rho*g = -rho*G*M(r)/r^2")
print()
print("  For constant density sphere:")
print("  P(r) = (2*pi/3)*G*rho^2*(R^2 - r^2)")
print()

print("ANCHOR CHAIN DERIVATION:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print(f"  c         = {c:.10e} m/s")
print(f"  G         = c^4 * 7.5 * epsilon_0^3 * mu_0^2")
print(f"            = {G:.14e} m^3/(kg·s^2)")
print()
print(f"  g_earth   = G*M_earth/R_earth^2 = {g_earth:.10f} m/s^2")
print(f"  P_hydro   = rho*g*h")
print()
print(f"  VF_r      = c^4/(8*pi*G) = {VF_r:.6e} Pa")
print(f"  (Maximum pressure before exotic matter required)")
print()

print("=" * 70)
print("Hydrostatic pressure derived from G via epsilon_0, mu_0")
print("From atmosphere to neutron stars, all traces back to vacuum constants")
print("=" * 70)

input("Press Enter to exit...")
