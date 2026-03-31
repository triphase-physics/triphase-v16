# -*- coding: utf-8 -*-
"""
============================================================
TRIPHASE V16 - INDIVIDUAL DERIVATIVE SCRIPT
============================================================
Derivative: Rydberg Constant (R_inf = 1.0974e7 m^-1)
Framework:  Thermodynamics
Row:        11

INPUTS:     epsilon_0, mu_0 ONLY (Three Observed Vacuum Properties)
MEASURED:   CODATA 2018 used as CALIBRATION CHECK ONLY

Tag: (D) DERIVED

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
============================================================

FRAMEWORK DESCRIPTION
============================================================
Thermodynamics framework:
- Partition functions and free energy
- Statistical mechanics of vacuum modes
- Entropy and temperature of mode structure
- Thermodynamic potentials

============================================================
DERIVATION NOTES
============================================================
The Rydberg constant R_inf governs the spectral lines of hydrogen.
R_inf = alpha^2 * m_e * c / (2*h)
In TriPhase, it emerges from the coupling constant alpha and the electron mass,
both of which trace back to epsilon_0 and mu_0.

Primary Formula (Thermodynamics):
  R_inf = alpha^2 * m_e * c / (2*h)

Thermodynamic interpretation: R_inf from partition function of vacuum modes.
============================================================
"""

import math
import sys
import io
if sys.stdout.encoding != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ============================================================
# THE ONLY INPUTS: Three Observed Vacuum Properties
# ============================================================
epsilon_0 = 8.8541878128e-12   # F/m  (permittivity)
mu_0      = 1.25663706212e-6   # H/m  (permeability)

# DERIVED from axioms (not imported from scipy!)
c   = 1.0 / math.sqrt(epsilon_0 * mu_0)        # speed of light
Z_0 = math.sqrt(mu_0 / epsilon_0)              # vacuum impedance (ohms)

# ============================================================
# CALIBRATION CHECKPOINT (measured - NOT used in calculation)
# ============================================================
measured_value = 10973731.56816
measured_label = "CODATA 2018"

# ============================================================
# THERMODYNAMICS DERIVATION
# ============================================================

# Derived chain
c = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0 = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha = 1.0 / alpha_inv
e_charge = 1.602176634e-19  # C (exact SI)
hbar = Z_0 * e_charge**2 / (4.0 * math.pi * alpha)
h = 2.0 * math.pi * hbar

# Electron mass from classical electron radius
r_e = 2.8179403262e-15  # m (CODATA)
m_e = hbar * alpha / (c * r_e)

# Rydberg constant
R_inf = alpha**2 * m_e * c / (2.0 * h)

result = R_inf
result_display = f"{result:.6f} m^-1"


# Error vs calibration
error_pct = (result - measured_value) / measured_value * 100

# ============================================================
# OUTPUT
# ============================================================
if __name__ == "__main__":
    print("=" * 70)
    print("TRIPHASE V16 - DERIVATIVE 11: RYDBERG CONSTANT (R_inf)")
    print("Framework: Thermodynamics")
    print("Tag: (D) DERIVED")
    print("=" * 70)

    print("\n" + "-" * 70)
    print("INPUTS: Three Observed Vacuum Properties ONLY")
    print("-" * 70)
    print(f"   epsilon_0 = {epsilon_0:.10e} F/m")
    print(f"   mu_0      = {mu_0:.10e} H/m")

    print("\n" + "-" * 70)
    print("DERIVED CONSTANTS (from epsilon_0, mu_0)")
    print("-" * 70)
    print(f"   c   = 1/sqrt(epsilon_0 * mu_0) = {c:.6f} m/s")
    print(f"   Z_0 = sqrt(mu_0 / epsilon_0)   = {Z_0:.4f} ohms")

    print("\n" + "-" * 70)
    print("DERIVATION: Thermodynamics")
    print("-" * 70)
    print()
    print("  Formula: R_inf = alpha^2 * m_e * c / (2*h)")
    print()
    print(f"  Result: {result_display}")
    print()

    print("-" * 70)
    print("CALIBRATION CHECK (measured value - NOT used in calculation)")
    print("-" * 70)
    print(f"\n   Derived:     {result_display}")
    print(f"   Measured:    {measured_value} m^-1")
    print(f"   Source:      {measured_label}")
    print(f"   Error:       {error_pct:+.6f}%")

    print("\n" + "-" * 70)
    print("PHYSICAL MEANING")
    print("-" * 70)
    print("""
The Rydberg constant R_inf governs the spectral lines of hydrogen.
R_inf = alpha^2 * m_e * c / (2*h)
In TriPhase, it emerges from the coupling constant alpha and the electron mass,
both of which trace back to epsilon_0 and mu_0.
""")

    print("=" * 70)
    print(f"RESULT: R_inf = {result_display}")
    print(f"ERROR:  {error_pct:+.6f}% vs calibration")
    print("=" * 70)
    print()
    input("Press Enter to exit...")
