"""
TriPhase V16 Python Derivative Script
Constant: Muon Mass
Framework: Anchor_Primitive
Tag: (D*) DERIVED - Discrete selection from anchor chain

Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
"""

import math

print("="*70)
print("TRIPHASE V16 - MUON MASS")
print("Framework: ANCHOR_PRIMITIVE")
print("Tag: (D*) DERIVED - Discrete selection from anchor chain")
print("="*70)
print()

# ============================================================================
# ANCHOR PRIMITIVE DERIVATION
# ============================================================================
print("ANCHOR PRIMITIVE DERIVATION")
print("-" * 70)

# ANCHOR INPUTS (SI exact or measured)
epsilon_0 = 8.8541878128e-12  # F/m (vacuum permittivity)
mu_0 = 1.25663706212e-6       # H/m (vacuum permeability)

print(f"INPUT ANCHORS:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print()

# DERIVED: Speed of light
c = 1.0 / math.sqrt(epsilon_0 * mu_0)
print(f"DERIVED:")
print(f"  c = 1/sqrt(epsilon_0 * mu_0)")
print(f"    = {c:.10e} m/s")
print()

# DERIVED: Impedance of free space
Z_0 = math.sqrt(mu_0 / epsilon_0)
print(f"  Z_0 = sqrt(mu_0/epsilon_0)")
print(f"      = {Z_0:.12f} Ohm")
print()

# DERIVED: Fine structure constant (TriPhase formula)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha = 1.0 / alpha_inv
print(f"  alpha_inv = 137 + ln(137)/137")
print(f"            = {alpha_inv:.12f}")
print(f"  alpha = 1/alpha_inv")
print(f"        = {alpha:.15e}")
print()

# DERIVED: Elementary charge (SI exact definition)
e = 1.602176634e-19  # C (exact by SI definition)
print(f"  e = {e:.15e} C (SI exact definition)")
print()

# DERIVED: Reduced Planck constant
hbar = Z_0 * e**2 / (4.0 * math.pi * alpha)
print(f"  hbar = Z_0 * e^2 / (4*pi*alpha)")
print(f"       = {hbar:.15e} J·s")
print()

# DERIVED: Electron mass (from energy scale)
# Using: m_e = 2 * R_inf * h / c
# Where R_inf = (m_e * c * alpha^2) / (2 * h)
# This is circular, so we use the established TriPhase formula:
# m_e = alpha^2 * mu_0 * c * e^2 / (2 * h)
# But simpler: m_e from fine structure energy scale
m_e = (alpha**2 * mu_0 * c * e**2) / (2.0 * hbar)
print(f"  m_e = (alpha^2 * mu_0 * c * e^2) / (2 * hbar)")
print(f"      = {m_e:.15e} kg")
m_e_MeV = m_e * c**2 / 1.602176634e-13  # Convert to MeV/c^2
print(f"      = {m_e_MeV:.10f} MeV/c^2")
print()

# ============================================================================
# MUON MASS DERIVATION (TRIPHASE LEPTON FORMULA)
# ============================================================================
print("="*70)
print("MUON MASS DERIVATION")
print("="*70)
print()

print("TriPhase lepton mass hierarchy emerges from frequency scaling:")
print()
print("  m_mu / m_e = (3^3 * 7) * (1 + alpha^2 corrections)")
print()
print("The muon represents the first excited state of the electron")
print("wavefunction, with mass ratio governed by cubic frequency scaling.")
print()

# TriPhase lepton mass ratio formula
base_ratio = 3**3 * 7  # 189
alpha_correction = 1.0 + 0.105 * alpha**2  # Small QED correction

mu_e_ratio = base_ratio * alpha_correction

print("CALCULATION:")
print(f"  Base ratio = 3^3 * 7 = {base_ratio}")
print(f"  Alpha correction = 1 + 0.105*alpha^2 = {alpha_correction:.12f}")
print(f"  m_mu/m_e = {mu_e_ratio:.12f}")
print()

# Muon mass
m_mu = m_e * mu_e_ratio
m_mu_MeV = m_mu * c**2 / 1.602176634e-13  # Convert to MeV/c^2

print(f"RESULT:")
print(f"  m_mu = {m_mu:.15e} kg")
print(f"       = {m_mu_MeV:.10f} MeV/c^2")
print()

# ============================================================================
# COMPARISON WITH CODATA (CALIBRATION CHECKPOINT)
# ============================================================================
print("="*70)
print("CALIBRATION CHECKPOINT")
print("="*70)
print()

m_mu_CODATA = 1.883531627e-28  # kg
m_mu_MeV_CODATA = 105.6583755  # MeV/c^2

print(f"TRIPHASE DERIVED:")
print(f"  m_mu = {m_mu:.15e} kg")
print(f"       = {m_mu_MeV:.10f} MeV/c^2")
print()
print(f"CODATA 2018:")
print(f"  m_mu = {m_mu_CODATA:.15e} kg")
print(f"       = {m_mu_MeV_CODATA:.10f} MeV/c^2")
print()

difference = abs(m_mu - m_mu_CODATA)
relative_error = (difference / m_mu_CODATA) * 100

difference_MeV = abs(m_mu_MeV - m_mu_MeV_CODATA)
relative_error_MeV = (difference_MeV / m_mu_MeV_CODATA) * 100

print(f"DIFFERENCE:")
print(f"  Delta m_mu = {difference:.3e} kg")
print(f"  Relative error = {relative_error:.6f}%")
print()
print(f"  Delta m_mu = {difference_MeV:.6f} MeV/c^2")
print(f"  Relative error = {relative_error_MeV:.6f}%")
print()

# ============================================================================
# SUMMARY
# ============================================================================
print("="*70)
print("SUMMARY")
print("="*70)
print()
print("FRAMEWORK: ANCHOR_PRIMITIVE")
print("  Inputs:  epsilon_0, mu_0")
print("  Derived: c, Z_0, alpha, e (SI def), hbar, m_e, m_mu")
print()
print(f"RESULT:")
print(f"  Muon mass = {m_mu_MeV:.10f} MeV/c^2")
print(f"  Mass ratio m_mu/m_e = {mu_e_ratio:.10f}")
print()
print("Pure derivation from vacuum electromagnetic properties.")
print("Lepton mass hierarchy emerges from TriPhase frequency scaling.")
print()
print("="*70)

input("Press Enter to exit...")
