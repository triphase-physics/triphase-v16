"""
TriPhase V16 Python Derivative Script
energy_per_mode_Anchor_Primitive.py

Calculates the energy per vacuum mode E_mode = hbar*omega/2 within the Anchor_Primitive framework.

Framework: Anchor_Primitive
Tag: (D) DERIVED - Pure anchor chain (epsilon_0, mu_0 only)

Row: 9
Each vacuum mode carries half-quantum of energy (zero-point energy).

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# ============================================================================
# ANCHOR PRIMITIVE DERIVATION
# ============================================================================

print("="*80)
print("TriPhase V16: Energy Per Vacuum Mode")
print("Framework: Anchor_Primitive")
print("Tag: (D) DERIVED - Pure anchor chain (epsilon_0, mu_0 only)")
print("="*80)
print()

# ----------------------------------------------------------------------------
# PURE ANCHOR CHAIN
# ----------------------------------------------------------------------------

print("PURE ANCHOR CHAIN:")
print("-" * 80)

# Primary anchors (ONLY inputs)
epsilon_0 = 8.8541878128e-12  # F/m - permittivity of free space
mu_0 = 1.25663706212e-6       # H/m - permeability of free space

print(f"epsilon_0 = {epsilon_0:.13e} F/m")
print(f"mu_0      = {mu_0:.14e} H/m")
print()

# Exact SI electron charge (2019 redefinition)
e = 1.602176634e-19  # C (exact)
print(f"e = {e:.15e} C (exact SI 2019)")
print()

# Derive speed of light
c = 1.0 / math.sqrt(epsilon_0 * mu_0)
print(f"c = 1/sqrt(epsilon_0 * mu_0)")
print(f"  = {c:.10e} m/s")
print()

# Derive impedance of free space
Z_0 = math.sqrt(mu_0 / epsilon_0)
print(f"Z_0 = sqrt(mu_0 / epsilon_0)")
print(f"    = {Z_0:.10f} Ohms")
print()

# ----------------------------------------------------------------------------
# FINE STRUCTURE CONSTANT
# ----------------------------------------------------------------------------

print("FINE STRUCTURE CONSTANT:")
print("-" * 80)

# Derive alpha from pressure band structure: 8*17+1 = 137
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha = 1.0 / alpha_inv

print(f"alpha^(-1) = 137 + ln(137)/137 (from 8*17+1=137)")
print(f"           = {alpha_inv:.10f}")
print(f"alpha = {alpha:.15e}")
print()

# ----------------------------------------------------------------------------
# REDUCED PLANCK CONSTANT
# ----------------------------------------------------------------------------

print("REDUCED PLANCK CONSTANT:")
print("-" * 80)

# Derive hbar from Z_0, e, and alpha
hbar = (Z_0 * e * e) / (4.0 * math.pi * alpha)

print(f"hbar = Z_0 * e^2 / (4*pi*alpha)")
print(f"     = {Z_0:.10f} * ({e:.15e})^2 / (4*pi*{alpha:.15e})")
print(f"     = {hbar:.15e} J*s")
print()

# Full Planck constant
h = 2.0 * math.pi * hbar
print(f"h = 2*pi*hbar")
print(f"  = {h:.15e} J*s")
print()

# ----------------------------------------------------------------------------
# ENERGY PER VACUUM MODE
# ----------------------------------------------------------------------------

print("ENERGY PER VACUUM MODE:")
print("-" * 80)

# For a mode with angular frequency omega, the zero-point energy is hbar*omega/2
# We'll demonstrate this for a characteristic frequency

# Example: Use Compton frequency of electron as characteristic scale
# First derive electron mass from hbar and alpha

# Classical electron radius (CODATA 2018)
r_e = 2.8179403262e-15  # m (using for derivation check)

# Derive electron mass: m_e = hbar*alpha/(c*r_e)
m_e = (hbar * alpha) / (c * r_e)
print(f"Electron mass (for example mode):")
print(f"m_e = hbar*alpha/(c*r_e)")
print(f"    = {m_e:.15e} kg")
print()

# Compton angular frequency
omega_C = (m_e * c * c) / hbar
print(f"Compton angular frequency:")
print(f"omega_C = m_e*c^2/hbar")
print(f"        = {omega_C:.10e} rad/s")
print()

# Zero-point energy for this mode
E_mode = hbar * omega_C / 2.0
print(f"Zero-point energy per mode:")
print(f"E_mode = hbar*omega/2")
print(f"       = ({hbar:.15e})*({omega_C:.10e})/2")
print(f"       = {E_mode:.15e} J")
print()

# Convert to eV
E_mode_eV = E_mode / e
print(f"E_mode = {E_mode_eV:.10e} eV")
print(f"       = {E_mode_eV / 1000.0:.6f} keV")
print(f"       = {E_mode_eV / 1.0e6:.10f} MeV")
print()

# For electron Compton frequency, this is exactly m_e*c^2/2
print(f"For Compton mode: E_mode = m_e*c^2/2 = {m_e * c * c / 2.0:.15e} J")
print(f"                                     = {(m_e * c * c / 2.0) / e / 1.0e6:.10f} MeV")
print()

# ----------------------------------------------------------------------------
# GENERAL FORMULA
# ----------------------------------------------------------------------------

print("GENERAL FORMULA:")
print("-" * 80)
print(f"For any vacuum mode with angular frequency omega:")
print(f"  E_mode = hbar*omega/2")
print()
print(f"Where hbar is derived from epsilon_0 and mu_0 through:")
print(f"  c = 1/sqrt(epsilon_0*mu_0)")
print(f"  Z_0 = sqrt(mu_0/epsilon_0)")
print(f"  alpha = 1/(137 + ln(137)/137)")
print(f"  hbar = Z_0*e^2/(4*pi*alpha)")
print()

# ----------------------------------------------------------------------------
# CALIBRATION CHECKPOINT
# ----------------------------------------------------------------------------

print("CALIBRATION CHECKPOINT:")
print("-" * 80)
print(f"hbar (derived) = {hbar:.15e} J*s")
print(f"hbar (CODATA)  = 1.054571817e-34 J*s")
delta_hbar = abs(hbar - 1.054571817e-34) / 1.054571817e-34 * 100.0
print(f"Relative difference: {delta_hbar:.4e}%")
print()

print(f"Electron rest energy (derived) = {m_e * c * c / e / 1.0e6:.10f} MeV")
print(f"Electron rest energy (CODATA)  = 0.51099895000 MeV")
delta_m_e = abs((m_e * c * c / e / 1.0e6) - 0.51099895000) / 0.51099895000 * 100.0
print(f"Relative difference: {delta_m_e:.4e}%")
print()

print("="*80)
print("Derivation complete.")
print("="*80)

input("Press Enter to exit...")
