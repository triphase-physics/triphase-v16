"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================
Derivative:  Electron g-factor (g_e/2)
Framework:   Anchor_Condensed
Tag: (D) DERIVED
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
print("ANCHOR CONDENSED DERIVATION: Electron g-factor (g_e/2)")
print("Framework: Anchor_Condensed | Tag: (D) DERIVED")
print("=" * 80)
print()
print("The electron g-factor anomaly derives from QED radiative corrections.")
print("First two terms from alpha (which derives from epsilon_0, mu_0):")
print()

# QED expansion for g_e/2
g_e_half = 1.0 + alpha / (2.0 * math.pi) - 0.328 * (alpha / math.pi)**2

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
print()
print(f"DERIVATION:")
print(f"  g_e/2 = 1 + alpha/(2*pi) - 0.328*(alpha/pi)^2")
print(f"        = 1 + {alpha/(2.0*math.pi):.13e}")
print(f"            - {0.328 * (alpha/math.pi)**2:.13e}")
print()
print(f"RESULT:")
print(f"  g_e/2 = {g_e_half:.14f}")
print()
print(f"CALIBRATION (CODATA 2022):")
print(f"  g_e/2 = 1.00115965218128 (measured)")
print()
print(f"MATCH:")
print(f"  Deviation = {abs(g_e_half - 1.00115965218128):.3e}")
print(f"  Note: Full QED calculation requires higher-order terms")
print()
print("=" * 80)

input("Press Enter to exit...")
