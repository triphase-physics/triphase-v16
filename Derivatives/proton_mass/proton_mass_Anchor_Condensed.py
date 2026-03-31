"""
proton_mass_Anchor_Condensed.py
TriPhase V16 - Row 26 - Tag: (D*) DERIVED with discrete selection

Proton mass from electron mass and discrete geometric factor.
m_p = m_e * 2^2 * 3^3 * 17 * (1 + 5*alpha^2/pi)

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# Anchor constants
epsilon_0 = 8.8541878128e-12
mu_0      = 1.25663706212e-6
e         = 1.602176634e-19

# Chain derivation
c     = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0   = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha = 1.0 / alpha_inv
hbar  = Z_0 * e**2 / (4.0 * math.pi * alpha)
m_e   = hbar * alpha / (c * 2.8179403262e-15)

# Proton-to-electron mass ratio
mp_me = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)

# Proton mass
m_p = m_e * mp_me

# Convert to MeV
m_p_MeV = m_p * c**2 / (1.602176634e-19 * 1e6)

print("=" * 70)
print("PROTON MASS - Anchor Condensed")
print("=" * 70)
print("\nFormula:")
print("  m_p = m_e * 2^2 * 3^3 * 17 * (1 + 5*alpha^2/pi)")
print("  mp/me = 4 * 27 * 17 * (1 + 5*alpha^2/pi)")
print("\nComponents:")
print(f"  m_e           = {m_e:.12e} kg")
print(f"  alpha         = {alpha:.12f}")
print(f"  mp/me         = {mp_me:.12f}")
print("\nDerived:")
print(f"  m_p           = {m_p:.12e} kg")
print(f"  m_p           = {m_p_MeV:.6f} MeV/c^2")
print("\nCODATA 2022:")
print(f"  m_p           = 1.67262192369e-27 kg")
print(f"  m_p           = 938.272088 MeV/c^2")
print(f"  mp/me         = 1836.15267")
print("\nAgreement:")
print(f"  Mass ratio    = {abs(mp_me - 1836.15267) / 1836.15267 * 100:.6f}%")
print(f"  Mass (kg)     = {abs(m_p - 1.67262192369e-27) / 1.67262192369e-27 * 100:.6f}%")
print(f"  Mass (MeV)    = {abs(m_p_MeV - 938.272088) / 938.272088 * 100:.6f}%")
print("=" * 70)

input("\nPress Enter to exit...")
