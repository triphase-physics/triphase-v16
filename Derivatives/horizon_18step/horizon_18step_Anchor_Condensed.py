"""
horizon_18step_Anchor_Condensed.py
TriPhase V16 - Row 31 - Tag: (D*) DERIVED with discrete selection

Hubble horizon from 18-step alpha cascade.
H_0 = pi * sqrt(3) * f_e * alpha^18
d_H = c / H_0

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
f_e   = m_e * c**2 / hbar

# Hubble constant formula (18-step cascade)
H_0 = math.pi * math.sqrt(3.0) * f_e * alpha**18

# Hubble horizon
d_H = c / H_0

# Convert to various units
d_H_m = d_H
d_H_km = d_H / 1e3
d_H_Mpc = d_H / 3.0857e22  # 1 Mpc = 3.0857e22 m
d_H_Gly = d_H / 9.4607e24  # 1 Gly = 9.4607e24 m
d_H_Gpc = d_H / 3.0857e25  # 1 Gpc = 3.0857e25 m

# Hubble constant in km/s/Mpc
H_0_kmsMpc = H_0 * 3.0857e22 / 1e3

print("=" * 70)
print("HUBBLE HORIZON (18-STEP CASCADE) - Anchor Condensed")
print("=" * 70)
print("\nFormula:")
print("  H_0 = pi * sqrt(3) * f_e * alpha^18")
print("  d_H = c / H_0")
print("\nComponents:")
print(f"  f_e           = {f_e:.12e} Hz")
print(f"  alpha         = {alpha:.12f}")
print(f"  alpha^18      = {alpha**18:.12e}")
print(f"  pi*sqrt(3)    = {math.pi * math.sqrt(3.0):.12f}")
print("\nDerived:")
print(f"  H_0           = {H_0:.12e} s^-1")
print(f"  H_0           = {H_0_kmsMpc:.6f} km/s/Mpc")
print(f"  d_H           = {d_H:.12e} m")
print(f"  d_H           = {d_H_km:.12e} km")
print(f"  d_H           = {d_H_Mpc:.6f} Mpc")
print(f"  d_H           = {d_H_Gly:.6f} Gly (billion light-years)")
print(f"  d_H           = {d_H_Gpc:.6f} Gpc")
print("\n18-Step Interpretation:")
print("  The alpha^18 cascade represents 18 discrete coupling steps")
print("  from Planck scale to cosmological horizon scale.")
print("\nObservational Context:")
print("  Planck 2018: H_0 ~ 67.4 km/s/Mpc (CMB-based)")
print("  SH0ES 2022:  H_0 ~ 73.0 km/s/Mpc (local distance ladder)")
print("  Observable universe radius: ~46.5 Gly (comoving)")
print("  Hubble radius (horizon): ~14.5 Gly")
print("=" * 70)

input("\nPress Enter to exit...")
