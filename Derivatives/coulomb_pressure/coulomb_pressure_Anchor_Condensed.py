"""
coulomb_pressure_Anchor_Condensed.py

TriPhase V16 - Coulomb Pressure
Row 38 - Tag: (D) DERIVED

Derives Coulomb pressure from electrostatic energy density.
P_Coulomb ~ e^2/(8*pi*epsilon_0*r^4) * n at atomic scales

All parameters (e, epsilon_0, hbar, m_e, alpha) from anchor chain.
Bohr radius a_0 = hbar/(m_e*c*alpha) emerges naturally.

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
print("TriPhase V16 - Coulomb Pressure (Electrostatic)")
print("Row 38 - DERIVED from epsilon_0, mu_0")
print("=" * 70)
print()

# Bohr radius (emerges from anchor chain)
a_0 = hbar / (m_e * c * alpha)
a_0_CODATA = 5.29177210903e-11  # m

print("BOHR RADIUS (from anchor chain):")
print(f"  a_0 = hbar/(m_e*c*alpha)")
print(f"      = {a_0:.14e} m")
print(f"  CODATA: {a_0_CODATA:.14e} m")
print(f"  Agreement: {100.0*a_0/a_0_CODATA:.6f}%")
print()

# Coulomb pressure formula
print("COULOMB PRESSURE:")
print("  Electrostatic energy density: u_E = epsilon_0*E^2/2")
print("  Electric field from point charge: E = e/(4*pi*epsilon_0*r^2)")
print("  Energy density: u_E = e^2/(32*pi^2*epsilon_0*r^4)")
print("  Pressure: P_Coulomb ~ u_E")
print()

# Pressure at Bohr radius
E_bohr = e / (4.0 * math.pi * epsilon_0 * a_0**2)
u_E_bohr = epsilon_0 * E_bohr**2 / 2.0
P_Coulomb_bohr = u_E_bohr

print("COULOMB PRESSURE AT BOHR RADIUS:")
print(f"  r = a_0 = {a_0:.6e} m")
print(f"  E(a_0) = e/(4*pi*epsilon_0*a_0^2) = {E_bohr:.6e} V/m")
print(f"  u_E = epsilon_0*E^2/2 = {u_E_bohr:.6e} J/m^3")
print(f"  P_Coulomb = {P_Coulomb_bohr:.6e} Pa")
print(f"            = {P_Coulomb_bohr/1e9:.3f} GPa")
print()

# Compare to vacuum frame rigidity
print("COMPARISON TO VACUUM FRAME RIGIDITY:")
print(f"  VF_r = {VF_r:.6e} Pa")
print(f"  P_Coulomb(a_0)/VF_r = {P_Coulomb_bohr/VF_r:.6e}")
print()
print("  Coulomb pressure at atomic scale is ~30 orders below VF_r")
print()

# Pressure at various scales
scales = [
    ("Nuclear scale", 1.0e-15),
    ("Bohr radius", a_0),
    ("Molecular scale", 1.0e-9),
    ("Micron scale", 1.0e-6),
]

print("COULOMB PRESSURE AT VARIOUS SCALES:")
for name, r in scales:
    E = e / (4.0 * math.pi * epsilon_0 * r**2)
    u_E = epsilon_0 * E**2 / 2.0
    P = u_E
    ratio = P / VF_r
    print(f"  {name:20s} (r={r:.2e} m): P={P:.6e} Pa, P/VF_r={ratio:.6e}")
print()

# Density-based Coulomb pressure
n_atom = 1.0 / a_0**3  # atomic number density (rough)
P_density = n_atom * e**2 / (4.0 * math.pi * epsilon_0 * a_0)

print("DENSITY-BASED COULOMB PRESSURE:")
print(f"  Atomic number density: n ~ 1/a_0^3 = {n_atom:.6e} m^-3")
print(f"  Coulomb energy per particle: E_c ~ e^2/(4*pi*epsilon_0*a_0)")
print(f"                                    = {e**2/(4.0*math.pi*epsilon_0*a_0):.6e} J")
print(f"  Pressure: P ~ n*E_c/a_0 = {P_density:.6e} Pa")
print()

# Electron degeneracy pressure (related but different)
print("RELATED: ELECTRON DEGENERACY PRESSURE:")
print("  (Not pure Coulomb, but relevant for dense matter)")
print(f"  P_deg ~ (hbar^2/m_e)*(n^(5/3))")
print(f"  At atomic density n ~ {n_atom:.2e} m^-3:")
P_deg = (hbar**2 / m_e) * n_atom**(5.0/3.0)
print(f"    P_deg ~ {P_deg:.6e} Pa = {P_deg/1e9:.3f} GPa")
print()

# White dwarf Coulomb pressure
n_wd = 1.0e36  # m^-3 (white dwarf density)
r_wd = (3.0/(4.0*math.pi*n_wd))**(1.0/3.0)
E_wd = e / (4.0 * math.pi * epsilon_0 * r_wd**2)
P_Coulomb_wd = epsilon_0 * E_wd**2 / 2.0

print("WHITE DWARF COULOMB PRESSURE:")
print(f"  Electron density: n ~ {n_wd:.2e} m^-3")
print(f"  Characteristic spacing: r ~ {r_wd:.6e} m")
print(f"  E(r) = {E_wd:.6e} V/m")
print(f"  P_Coulomb ~ {P_Coulomb_wd:.6e} Pa = {P_Coulomb_wd/1e9:.3e} GPa")
print(f"  P/VF_r = {P_Coulomb_wd/VF_r:.6e}")
print()

# Rydberg energy connection
E_Ryd = m_e * c**2 * alpha**2 / 2.0
E_Ryd_eV = E_Ryd / e

print("RYDBERG ENERGY CONNECTION:")
print(f"  E_Ryd = m_e*c^2*alpha^2/2 = {E_Ryd:.6e} J")
print(f"        = {E_Ryd_eV:.6f} eV")
print(f"  Also: E_Ryd = e^2/(8*pi*epsilon_0*a_0) = {e**2/(8.0*math.pi*epsilon_0*a_0):.6e} J")
print(f"  Pressure scale: P ~ E_Ryd/a_0^3 = {E_Ryd/a_0**3:.6e} Pa")
print()

print("ANCHOR CHAIN DERIVATION:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print(f"  e         = {e:.15e} C (exact SI)")
print(f"  c         = {c:.10e} m/s")
print(f"  alpha     = {alpha:.15e}")
print(f"  hbar      = {hbar:.15e} J·s")
print(f"  m_e       = {m_e:.15e} kg")
print(f"  a_0       = hbar/(m_e*c*alpha) = {a_0:.14e} m")
print()
print(f"  E(a_0)    = e/(4*pi*epsilon_0*a_0^2) = {E_bohr:.6e} V/m")
print(f"  P_Coulomb = epsilon_0*E^2/2 = {P_Coulomb_bohr:.6e} Pa")
print()
print(f"  VF_r      = {VF_r:.6e} Pa")
print()

print("=" * 70)
print("Coulomb pressure derived from e, epsilon_0 via anchor chain")
print("Atomic-scale electrostatic pressure from first principles")
print("=" * 70)

input("Press Enter to exit...")
