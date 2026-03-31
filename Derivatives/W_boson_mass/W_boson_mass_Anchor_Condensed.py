"""
W_boson_mass_Anchor_Condensed.py
TriPhase V16 - Row 28 - Tag: (D*H) DERIVED but hypothetical

W boson mass from proton mass and fine structure constant.
M_W = m_p * alpha^(-1) / 2

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

# W boson mass formula
M_W = m_p * alpha_inv / 2.0

# Convert to GeV
M_W_GeV = M_W * c**2 / (1.602176634e-19 * 1e9)

print("=" * 70)
print("W BOSON MASS - Anchor Condensed")
print("=" * 70)
print("\nFormula:")
print("  M_W = m_p * alpha^(-1) / 2")
print("  M_W = m_e * mp_me * alpha^(-1) / 2")
print("\nComponents:")
print(f"  m_p           = {m_p:.12e} kg")
print(f"  alpha         = {alpha:.12f}")
print(f"  alpha^(-1)    = {alpha_inv:.12f}")
print("\nDerived:")
print(f"  M_W           = {M_W:.12e} kg")
print(f"  M_W           = {M_W_GeV:.6f} GeV/c^2")
print("\nPDG 2022:")
print(f"  M_W           = 80.3692 +/- 0.0133 GeV/c^2")
print("\nAgreement:")
print(f"  Mass (GeV)    = {abs(M_W_GeV - 80.3692) / 80.3692 * 100:.6f}%")
print(f"  Difference    = {M_W_GeV - 80.3692:+.6f} GeV/c^2")
print("=" * 70)

input("\nPress Enter to exit...")
