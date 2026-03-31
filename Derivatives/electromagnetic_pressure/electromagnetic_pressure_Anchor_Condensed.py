"""
electromagnetic_pressure_Anchor_Condensed.py

TriPhase V16 - Electromagnetic Pressure
Row 35 - Tag: (D) DERIVED

Derives electromagnetic pressure from vacuum constants:
P_EM = epsilon_0*E^2/2 + B^2/(2*mu_0)

Shows maximum possible EM pressure equals vacuum frame rigidity:
P_max = c^4/(8*pi*G) = VF_r

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
print("TriPhase V16 - Electromagnetic Pressure")
print("Row 35 - DERIVED from epsilon_0, mu_0")
print("=" * 70)
print()

# Electromagnetic pressure formula
print("ELECTROMAGNETIC PRESSURE:")
print("  P_EM = epsilon_0*E^2/2 + B^2/(2*mu_0)")
print()
print("  Electric contribution: u_E = epsilon_0*E^2/2")
print("  Magnetic contribution: u_B = B^2/(2*mu_0)")
print()

# Unit field examples
E_unit = 1.0  # V/m
B_unit = 1.0  # T

P_E_unit = epsilon_0 * E_unit**2 / 2.0
P_B_unit = B_unit**2 / (2.0 * mu_0)

print("UNIT FIELD PRESSURES:")
print(f"  E = {E_unit:.1f} V/m:")
print(f"    P_E = epsilon_0*E^2/2 = {P_E_unit:.6e} Pa")
print()
print(f"  B = {B_unit:.1f} T:")
print(f"    P_B = B^2/(2*mu_0) = {P_B_unit:.6e} Pa")
print()

# Strong field examples
E_strong = 1.0e11  # V/m (near atomic breakdown)
B_strong = 10.0    # T (strong lab magnet)

P_E_strong = epsilon_0 * E_strong**2 / 2.0
P_B_strong = B_strong**2 / (2.0 * mu_0)

print("STRONG FIELD PRESSURES:")
print(f"  E = {E_strong:.1e} V/m (atomic breakdown):")
print(f"    P_E = {P_E_strong:.6e} Pa = {P_E_strong/1e5:.3f} atm")
print()
print(f"  B = {B_strong:.1f} T (strong lab magnet):")
print(f"    P_B = {P_B_strong:.6e} Pa = {P_B_strong/1e5:.3f} atm")
print()

# Maximum EM pressure
print("MAXIMUM EM PRESSURE (VACUUM FRAME RIGIDITY):")
print(f"  VF_r = c^4/(8*pi*G)")
print(f"       = {VF_r:.6e} Pa")
print()
print("  This is the maximum pressure the vacuum frame can sustain.")
print("  Higher pressures would require exotic matter/negative energy.")
print()

# Field strengths at VF_r
E_max = math.sqrt(2.0 * VF_r / epsilon_0)
B_max = math.sqrt(2.0 * mu_0 * VF_r)

print("FIELD STRENGTHS AT VACUUM RIGIDITY:")
print(f"  E_max = sqrt(2*VF_r/epsilon_0) = {E_max:.6e} V/m")
print(f"  B_max = sqrt(2*mu_0*VF_r)      = {B_max:.6e} T")
print()
print("  (These are unphysical - would require energy density > VF_r)")
print()

# Practical field ratios
print("PRACTICAL FIELD RATIOS TO VF_r:")
print(f"  Strong lab E ({E_strong:.1e} V/m): P/VF_r = {P_E_strong/VF_r:.6e}")
print(f"  Strong lab B ({B_strong:.1f} T):      P/VF_r = {P_B_strong/VF_r:.6e}")
print()
print("  Even extreme lab fields are ~40 orders of magnitude below VF_r!")
print()

# Energy density equivalence
u_E = epsilon_0 * E_strong**2 / 2.0
u_B = B_strong**2 / (2.0 * mu_0)

print("ENERGY DENSITY = PRESSURE:")
print(f"  For E = {E_strong:.1e} V/m:")
print(f"    u_E = {u_E:.6e} J/m^3 = P_E")
print()
print(f"  For B = {B_strong:.1f} T:")
print(f"    u_B = {u_B:.6e} J/m^3 = P_B")
print()
print("  In relativistic units: energy density = pressure")
print()

# Poynting vector and radiation pressure
S = E_strong * B_strong / mu_0
P_rad = S / c

print("RADIATION PRESSURE:")
print(f"  For crossed fields E = {E_strong:.1e} V/m, B = {B_strong:.1f} T:")
print(f"    Poynting S = E*B/mu_0 = {S:.6e} W/m^2")
print(f"    P_rad = S/c = {P_rad:.6e} Pa")
print()

print("ANCHOR CHAIN DERIVATION:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print(f"  c         = {c:.10e} m/s")
print(f"  G         = c^4 * 7.5 * epsilon_0^3 * mu_0^2")
print(f"            = {G:.14e} m^3/(kg·s^2)")
print(f"  VF_r      = c^4/(8*pi*G)")
print(f"            = {VF_r:.6e} Pa")
print()
print("  P_EM = epsilon_0*E^2/2 + B^2/(2*mu_0)")
print(f"  P_max = VF_r = {VF_r:.6e} Pa")
print()

print("=" * 70)
print("Electromagnetic pressure derived from epsilon_0, mu_0")
print("Maximum EM pressure is vacuum frame rigidity")
print("=" * 70)

input("Press Enter to exit...")
