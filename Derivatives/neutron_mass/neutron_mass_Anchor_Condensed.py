"""
neutron_mass_Anchor_Condensed.py
TriPhase V16 - Row 27 - Tag: (D*) DERIVED with discrete selection

Neutron mass from proton mass with alpha/(2*pi*T_17) correction.
m_n = m_p * (1 + alpha/(2*pi*T_17)) where T_17 = 153

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
mp_me = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p   = m_e * mp_me

# Neutron mass formula
T_17 = 153  # Discrete selection parameter
m_n = m_p * (1.0 + alpha / (2.0 * math.pi * T_17))

# Mass difference
Delta_m = m_n - m_p
Delta_m_MeV = Delta_m * c**2 / (1.602176634e-19 * 1e6)

# Convert to MeV
m_n_MeV = m_n * c**2 / (1.602176634e-19 * 1e6)
m_p_MeV = m_p * c**2 / (1.602176634e-19 * 1e6)

print("=" * 70)
print("NEUTRON MASS - Anchor Condensed")
print("=" * 70)
print("\nFormula:")
print("  m_n = m_p * (1 + alpha/(2*pi*T_17))")
print("  where T_17 = 153")
print("\nComponents:")
print(f"  m_p           = {m_p:.12e} kg")
print(f"  alpha         = {alpha:.12f}")
print(f"  T_17          = {T_17}")
print(f"  alpha/(2*pi*T_17) = {alpha / (2.0 * math.pi * T_17):.12e}")
print("\nDerived:")
print(f"  m_n           = {m_n:.12e} kg")
print(f"  m_n           = {m_n_MeV:.6f} MeV/c^2")
print(f"  Delta_m (n-p) = {Delta_m:.12e} kg")
print(f"  Delta_m (n-p) = {Delta_m_MeV:.6f} MeV/c^2")
print("\nCODATA 2022:")
print(f"  m_n           = 1.67492749804e-27 kg")
print(f"  m_n           = 939.565420 MeV/c^2")
print(f"  Delta_m (n-p) = 1.293 MeV/c^2")
print("\nAgreement:")
print(f"  Mass (kg)     = {abs(m_n - 1.67492749804e-27) / 1.67492749804e-27 * 100:.6f}%")
print(f"  Mass (MeV)    = {abs(m_n_MeV - 939.565420) / 939.565420 * 100:.6f}%")
print(f"  Delta_m       = {abs(Delta_m_MeV - 1.293) / 1.293 * 100:.6f}%")
print("=" * 70)

input("\nPress Enter to exit...")
