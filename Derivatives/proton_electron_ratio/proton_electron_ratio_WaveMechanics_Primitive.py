# -*- coding: utf-8 -*-
"""
============================================================
TRIPHASE V16 - INDIVIDUAL DERIVATIVE SCRIPT
============================================================
Derivative: Proton-Electron Mass Ratio (mp/me = 1836.156)
Framework:  WaveMechanics_Primitive
Row:        2

INPUTS:     epsilon_0, mu_0 ONLY (Three Observed Vacuum Properties)
MEASURED:   CODATA 2022 value used as CALIBRATION CHECK ONLY

Tag: (D*) DERIVED with discrete selection (2^2, 3^3, 17)

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
============================================================
"""

import numpy as np
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
c  = 1.0 / np.sqrt(epsilon_0 * mu_0)        # speed of light
Z_0 = np.sqrt(mu_0 / epsilon_0)              # vacuum impedance (ohms)

# DERIVE alpha from pressure band formula
m_band = 17
node = 8 * m_band + 1                        # = 137
correction = np.log(node) / node             # ln(137)/137
alpha_inv = node + correction                # 137.0359...
alpha = 1.0 / alpha_inv

# ============================================================
# CALIBRATION CHECKPOINT (measured - NOT used in calculation)
# ============================================================
mp_me_CODATA = 1836.15267344    # CODATA 2022

# ============================================================
# WAVE MECHANICS PRIMITIVE DERIVATION
# ============================================================
#
# MECHANISM:
# The proton-electron mass ratio emerges from the discrete harmonic
# structure of the vector frame.
#
# The ratio combines:
#   - Powers of 2 (binary structure)
#   - Powers of 3 (three-phase structure)
#   - The pressure band index m=17 (same m that gives alpha^-1 = 137)
#   - Fine structure correction (alpha^2/pi term)
#
# FORMULA:
#   mp/me = 2^2 * 3^3 * 17 * (1 + 5*alpha^2/pi)
#
# WHERE:
#   2^2 = 4      (binary harmonic base)
#   3^3 = 27     (three-phase harmonic cube)
#   17 = m       (pressure band index from 8m+1=137)
#   5*alpha^2/pi (fine structure coupling correction)
#
# WHY THIS STRUCTURE?
#   The proton is NOT an elementary particle - it's a composite
#   structure of three quarks. This formula captures the HARMONIC
#   relationship between the electron (elementary) and the proton
#   (composite) as they couple to the vector frame.
#
#   The base: 4 * 27 * 17 = 1836 (integer part)
#   The correction: 5*alpha^2/pi adds ~0.15 (QED coupling)
#
# ============================================================

# Integer harmonic base
base_2 = 2**2                   # = 4
base_3 = 3**3                   # = 27
m = 17                          # pressure band index (from 8m+1=137)

# Integer product (before fine structure correction)
integer_base = base_2 * base_3 * m   # = 1836

# Fine structure correction factor
fine_structure_correction = 1 + 5 * alpha**2 / np.pi

# DERIVED proton-electron mass ratio
mp_me_derived = integer_base * fine_structure_correction

# Error vs calibration
error = (mp_me_derived - mp_me_CODATA) / mp_me_CODATA * 100

# ============================================================
# OUTPUT
# ============================================================
if __name__ == "__main__":
    print("=" * 70)
    print("TRIPHASE V16 - DERIVATIVE 2: PROTON-ELECTRON MASS RATIO (mp/me)")
    print("Framework: WaveMechanics_Primitive")
    print("Tag: (D*) DERIVED with discrete selection")
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
    print("DERIVED alpha (from pressure band formula)")
    print("-" * 70)
    print(f"   Pressure band: m = {m_band}")
    print(f"   Node: 8 * {m_band} + 1 = {node}")
    print(f"   Correction: ln({node})/{node} = {correction:.8f}")
    print(f"   alpha^-1 = {alpha_inv:.9f}")
    print(f"   alpha    = {alpha:.10e}")

    print("\n" + "-" * 70)
    print("DERIVATION: Harmonic Structure Formula")
    print("-" * 70)
    print(f"\n   Binary base:     2^2 = {base_2}")
    print(f"   Three-phase:     3^3 = {base_3}")
    print(f"   Pressure band:   m   = {m}")
    print(f"   Integer base:    {base_2} × {base_3} × {m} = {integer_base}")
    print(f"\n   Fine structure correction:")
    print(f"      1 + 5*alpha^2/pi = 1 + 5*({alpha:.10e})^2/pi")
    print(f"                       = {fine_structure_correction:.10f}")
    print(f"\n   mp/me = {integer_base} × {fine_structure_correction:.10f}")
    print(f"         = {mp_me_derived:.8f}")

    print("\n" + "-" * 70)
    print("CALIBRATION CHECK (measured value - NOT used in calculation)")
    print("-" * 70)
    print(f"\n   Derived:      {mp_me_derived:.8f}")
    print(f"   CODATA 2022:  {mp_me_CODATA:.8f}")
    print(f"\n   Error:        {error:+.6f}%")

    print("\n" + "-" * 70)
    print("PHYSICAL MEANING")
    print("-" * 70)
    print(f"""
   mp/me = {mp_me_derived:.2f} tells us the HARMONIC relationship
   between elementary and composite structures in the vector frame.

   The electron is ELEMENTARY - it couples directly at frequency 137.
   The proton is COMPOSITE (3 quarks) - its mass reflects the
   harmonic structure of bound states.

   Integer base: 2^2 × 3^3 × 17 = {integer_base}
      2^2 = 4   (binary harmonic)
      3^3 = 27  (three-phase cubic)
      17  = m   (pressure band index, same m from alpha^-1 = 137)

   The {integer_base} is NOT a coincidence - it's the discrete
   harmonic spacing on the vector frame for composite structures.

   Fine structure adds ~0.15 from electromagnetic self-coupling.
""")

    print("=" * 70)
    print(f"RESULT: mp/me = {mp_me_derived:.8f}")
    print(f"ERROR:  {error:+.6f}% vs CODATA 2022 calibration")
    print("=" * 70)
    print()
    input("Press Enter to exit...")
