"""
TriPhase V16 Derivative - Proton Mass (ANCHOR_PRIMITIVE)
Framework: Anchor_Primitive
Tag: (D*)
DOI: 10.5281/zenodo.17855383

ANCHOR PRIMITIVE DERIVATION:
- ONLY inputs: epsilon_0 = 8.8541878128e-12 F/m, mu_0 = 1.25663706212e-6 H/m
- e = 1.602176634e-19 C (exact SI definition)
- c = 1/sqrt(epsilon_0 * mu_0)
- Z_0 = sqrt(mu_0/epsilon_0)
- alpha_inv = 137 + ln(137)/137, alpha = 1/alpha_inv
- hbar = Z_0 * e^2 / (4*pi*alpha)
- m_e from hbar, alpha, c chain
- m_p = m_e * mp_me_ratio
- mp_me_ratio = 2^2 * 3^3 * 17 * (1 + 5*alpha^2/pi) = 1836.15...

PURE ANCHOR CHAIN - No shortcuts, full derivation from electromagnetic constants.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
"""

import math

print("="*70)
print("TriPhase V16 - Proton Mass (m_p)")
print("Framework: Anchor_Primitive | Tag: (D*)")
print("="*70)
print()

# ============================================================================
# ANCHOR PRIMITIVE INPUTS (ONLY)
# ============================================================================
epsilon_0 = 8.8541878128e-12  # F/m (permittivity of free space)
mu_0 = 1.25663706212e-6       # H/m (permeability of free space)
e = 1.602176634e-19           # C (elementary charge, exact SI definition)

print("ANCHOR PRIMITIVE INPUTS:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print(f"  e         = {e:.12e} C (exact SI)")
print()

# ============================================================================
# DERIVED FUNDAMENTAL CONSTANTS
# ============================================================================
print("DERIVED CONSTANTS FROM ANCHOR CHAIN:")
print()

# Speed of light
c = 1.0 / math.sqrt(epsilon_0 * mu_0)
print(f"c = 1/sqrt(epsilon_0 * mu_0)")
print(f"  = {c:.10e} m/s")
print()

# Impedance of free space
Z_0 = math.sqrt(mu_0 / epsilon_0)
print(f"Z_0 = sqrt(mu_0/epsilon_0)")
print(f"    = {Z_0:.10f} Ohms")
print()

# Fine structure constant (TriPhase derivation)
alpha_inv = 137.0 + math.log(137.0)/137.0
alpha = 1.0 / alpha_inv
print(f"alpha_inv = 137 + ln(137)/137 = {alpha_inv:.10f}")
print(f"alpha     = 1/alpha_inv = {alpha:.12e}")
print()

# Reduced Planck constant
hbar = Z_0 * e**2 / (4.0 * math.pi * alpha)
print(f"hbar = Z_0 * e^2 / (4*pi*alpha)")
print(f"     = {hbar:.15e} J·s")
print()

# Compton wavelength of electron
lambda_C = 2.0 * math.pi / (alpha**2 * 137.0)
print(f"lambda_C = 2*pi / (alpha^2 * 137)")
print(f"         = {lambda_C:.15e} m")
print()

# Electron mass
m_e = hbar / (c * lambda_C)
print(f"m_e = hbar / (c * lambda_C)")
print(f"    = {m_e:.15e} kg")
print()

# ============================================================================
# PROTON MASS DERIVATION
# ============================================================================
print("="*70)
print("PROTON MASS DERIVATION:")
print("="*70)
print()

# Proton-to-electron mass ratio (TriPhase formula)
# mp_me_ratio = 2^2 * 3^3 * 17 * (1 + 5*alpha^2/pi)
base_factor = 2**2 * 3**3 * 17  # = 4 * 27 * 17 = 1836
correction_factor = 1.0 + 5.0 * alpha**2 / math.pi
mp_me_ratio = base_factor * correction_factor

print(f"mp_me_ratio = 2^2 * 3^3 * 17 * (1 + 5*alpha^2/pi)")
print(f"Base factor: 2^2 * 3^3 * 17 = {base_factor}")
print(f"Correction: (1 + 5*alpha^2/pi) = {correction_factor:.10f}")
print(f"mp_me_ratio = {mp_me_ratio:.10f}")
print()

# Proton mass
m_p = m_e * mp_me_ratio
print(f"m_p = m_e * mp_me_ratio")
print(f"    = {m_p:.15e} kg")
print()

# Convert to MeV/c^2
m_p_MeV = m_p * c**2 / (e * 1e6)
print(f"m_p = {m_p_MeV:.6f} MeV/c^2")
print()

# ============================================================================
# CODATA COMPARISON (CALIBRATION CHECKPOINT)
# ============================================================================
print("="*70)
print("CODATA COMPARISON (CALIBRATION CHECKPOINT):")
print("="*70)
m_p_CODATA = 1.67262192369e-27  # kg (CODATA 2018)
m_p_CODATA_MeV = 938.27208816    # MeV/c^2

print(f"TriPhase m_p:  {m_p:.15e} kg")
print(f"CODATA m_p:    {m_p_CODATA:.15e} kg")
print(f"Difference:    {abs(m_p - m_p_CODATA):.3e} kg")
print(f"Relative diff: {abs(m_p - m_p_CODATA)/m_p_CODATA * 100:.6f}%")
print()
print(f"TriPhase m_p:  {m_p_MeV:.6f} MeV/c^2")
print(f"CODATA m_p:    {m_p_CODATA_MeV:.6f} MeV/c^2")
print(f"Difference:    {abs(m_p_MeV - m_p_CODATA_MeV):.6f} MeV/c^2")
print()

print("="*70)
print("Derivation complete. Strong foundations built from anchor primitives.")
print("="*70)

input("Press Enter to exit...")
