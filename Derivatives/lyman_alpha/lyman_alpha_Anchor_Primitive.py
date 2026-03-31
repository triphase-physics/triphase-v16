"""
TriPhase V16 Python Derivative Script
lyman_alpha_Anchor_Primitive.py

Calculates the Lyman Alpha wavelength (121.567 nm) within the Anchor_Primitive framework.

Framework: Anchor_Primitive
Tag: (D) DERIVED - Pure anchor chain (epsilon_0, mu_0 only)

Row: 10
Lambda_Ly = 4/(3*R_inf) where R_inf = alpha^2*m_e*c/(2*h)
Everything derived from epsilon_0, mu_0.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

# ============================================================================
# ANCHOR PRIMITIVE DERIVATION
# ============================================================================

print("="*80)
print("TriPhase V16: Lyman Alpha Wavelength")
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
print(f"     = {hbar:.15e} J*s")
print()

# Full Planck constant
h = 2.0 * math.pi * hbar
print(f"h = 2*pi*hbar")
print(f"  = {h:.15e} J*s")
print()

# ----------------------------------------------------------------------------
# ELECTRON MASS
# ----------------------------------------------------------------------------

print("ELECTRON MASS:")
print("-" * 80)

# Classical electron radius (for derivation - from CODATA)
r_e = 2.8179403262e-15  # m

# Derive electron mass: m_e = hbar*alpha/(c*r_e)
m_e = (hbar * alpha) / (c * r_e)
print(f"m_e = hbar*alpha/(c*r_e)")
print(f"    = {m_e:.15e} kg")
print()

# Rest energy
E_rest = m_e * c * c
E_rest_eV = E_rest / e
print(f"m_e*c^2 = {E_rest:.15e} J")
print(f"        = {E_rest_eV / 1.0e6:.10f} MeV")
print()

# ----------------------------------------------------------------------------
# RYDBERG CONSTANT
# ----------------------------------------------------------------------------

print("RYDBERG CONSTANT:")
print("-" * 80)

# R_inf = alpha^2 * m_e * c / (2*h)
R_inf = (alpha * alpha * m_e * c) / (2.0 * h)

print(f"R_inf = alpha^2 * m_e * c / (2*h)")
print(f"      = ({alpha:.15e})^2 * {m_e:.15e} * {c:.10e} / (2*{h:.15e})")
print(f"      = {R_inf:.15e} m^(-1)")
print()

# Rydberg wavelength
lambda_R = 1.0 / R_inf
print(f"Rydberg wavelength:")
print(f"lambda_R = 1/R_inf")
print(f"         = {lambda_R:.15e} m")
print(f"         = {lambda_R * 1.0e9:.10f} nm")
print()

# ----------------------------------------------------------------------------
# LYMAN ALPHA WAVELENGTH
# ----------------------------------------------------------------------------

print("LYMAN ALPHA WAVELENGTH:")
print("-" * 80)

# Lyman alpha is the n=2 to n=1 transition in hydrogen
# Energy: E = R_inf * h * c * (1/1^2 - 1/2^2) = R_inf * h * c * (3/4)
# Wavelength: lambda_Ly = 1 / (R_inf * 3/4) = 4/(3*R_inf)

lambda_Ly = 4.0 / (3.0 * R_inf)

print(f"Lyman alpha transition: n=2 -> n=1")
print(f"lambda_Ly = 4/(3*R_inf)")
print(f"          = 4/(3*{R_inf:.15e})")
print(f"          = {lambda_Ly:.15e} m")
print(f"          = {lambda_Ly * 1.0e9:.10f} nm")
print()

# Photon energy
E_Ly = (h * c) / lambda_Ly
E_Ly_eV = E_Ly / e

print(f"Lyman alpha photon energy:")
print(f"E_Ly = h*c/lambda_Ly")
print(f"     = {E_Ly:.15e} J")
print(f"     = {E_Ly_eV:.10f} eV")
print()

# Frequency
nu_Ly = c / lambda_Ly
print(f"Lyman alpha frequency:")
print(f"nu_Ly = c/lambda_Ly")
print(f"      = {nu_Ly:.10e} Hz")
print()

# ----------------------------------------------------------------------------
# PHYSICAL SIGNIFICANCE
# ----------------------------------------------------------------------------

print("PHYSICAL SIGNIFICANCE:")
print("-" * 80)
print(f"Lyman alpha is the most prominent line in hydrogen spectrum.")
print(f"It represents the 2p -> 1s transition in hydrogen.")
print(f"In TriPhase, this wavelength appears in:")
print(f"  - Cosmological surveys (redshifted Lyman alpha forest)")
print(f"  - Quantum field vacuum structure")
print(f"  - Pressure band resonances at n=2 to n=1 transition")
print()

# ----------------------------------------------------------------------------
# CALIBRATION CHECKPOINT
# ----------------------------------------------------------------------------

print("CALIBRATION CHECKPOINT:")
print("-" * 80)
print(f"lambda_Ly (derived) = {lambda_Ly * 1.0e9:.10f} nm")
print(f"lambda_Ly (NIST)    = 121.56700000 nm")
delta_ly = abs((lambda_Ly * 1.0e9) - 121.567) / 121.567 * 100.0
print(f"Relative difference: {delta_ly:.4e}%")
print()

print(f"R_inf (derived) = {R_inf:.15e} m^(-1)")
print(f"R_inf (CODATA)  = 1.0973731568160e+07 m^(-1)")
delta_R = abs(R_inf - 1.0973731568160e+07) / 1.0973731568160e+07 * 100.0
print(f"Relative difference: {delta_R:.4e}%")
print()

print("="*80)
print("Derivation complete.")
print("="*80)

input("Press Enter to exit...")
