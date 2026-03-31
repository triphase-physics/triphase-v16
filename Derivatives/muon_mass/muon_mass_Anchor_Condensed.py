"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================
Derivative:  Muon Mass
Framework:   Anchor_Condensed
Tag: (D*) DERIVED with discrete selection
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
"""
import math

# === ANCHOR INPUTS ===
epsilon_0 = 8.8541878128e-12
mu_0      = 1.25663706212e-6
e         = 1.602176634e-19

# === DERIVED CHAIN ===
c = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0 = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha = 1.0 / alpha_inv
hbar = Z_0 * e**2 / (4.0 * math.pi * alpha)
h = 2.0 * math.pi * hbar
G = c**4 * 7.5 * epsilon_0**3 * mu_0**2
m_e = hbar * alpha / (c * 2.8179403262e-15)
f_e = m_e * c**2 / hbar

# === DERIVATION ===
print("=" * 80)
print("ANCHOR CONDENSED DERIVATION: Muon Mass")
print("Framework: Anchor_Condensed | Tag: (D*) DERIVED with discrete selection")
print("=" * 80)
print()
print("The muon mass scales from electron mass via alpha and integer factors:")
print()

# Muon mass ratio
m_mu_m_e_ratio = 3.0 / (2.0 * alpha) * (1.0 + alpha / (2.0 * math.pi))
m_mu = m_e * m_mu_m_e_ratio

# Convert to MeV/c^2
m_mu_MeV = m_mu * c**2 / (e * 1e6)

print(f"INPUTS (Anchor):")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print(f"  e         = {e:.12e} C (exact SI)")
print()
print(f"DERIVED CHAIN:")
print(f"  c         = {c:.10e} m/s")
print(f"  Z_0       = {Z_0:.10f} ohm")
print(f"  alpha_inv = {alpha_inv:.10f}")
print(f"  alpha     = {alpha:.13e}")
print(f"  hbar      = {hbar:.13e} J·s")
print(f"  m_e       = {m_e:.13e} kg")
print()
print(f"DERIVATION:")
print(f"  m_mu/m_e = 3/(2*alpha) * (1 + alpha/(2*pi))")
print(f"           = {3.0/(2.0*alpha):.6f} * {1.0 + alpha/(2.0*math.pi):.10f}")
print(f"           = {m_mu_m_e_ratio:.10f}")
print()
print(f"  m_mu = m_e * (m_mu/m_e)")
print(f"       = {m_e:.13e} kg * {m_mu_m_e_ratio:.10f}")
print()
print(f"RESULT:")
print(f"  m_mu = {m_mu:.13e} kg")
print(f"  m_mu = {m_mu_MeV:.10f} MeV/c^2")
print()
print(f"CALIBRATION (CODATA 2022):")
print(f"  m_mu = 1.883531627e-28 kg")
print(f"  m_mu = 105.6583755 MeV/c^2")
print(f"  m_mu/m_e = 206.7682830")
print()
print(f"MATCH:")
print(f"  Mass ratio deviation = {abs(m_mu_m_e_ratio - 206.7682830):.3e}")
print()
print("=" * 80)

input("Press Enter to exit...")
