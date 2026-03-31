# -*- coding: utf-8 -*-
"""
============================================================
TRIPHASE V16 - INDIVIDUAL DERIVATIVE SCRIPT
============================================================
Derivative: Speed of Light (c)
Framework:  WaveMechanics_Primitive
Row:        5

INPUTS:     epsilon_0, mu_0 ONLY (Three Observed Vacuum Properties)
MEASURED:   SI 2019 definition (exact) used as CALIBRATION CHECK ONLY

Tag: (D) DERIVED - Direct from the two axioms

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

# ============================================================
# CALIBRATION CHECKPOINT (SI definition - NOT used in calculation)
# ============================================================
c_SI = 299792458   # m/s (exact by SI definition since 1983)

# ============================================================
# WAVE MECHANICS PRIMITIVE DERIVATION
# ============================================================
#
# MECHANISM:
# The speed of light is the PROPAGATION RATE of electromagnetic
# waves through the vector frame.
#
# The vector frame has TWO intrinsic properties:
#    epsilon_0 (electric permittivity) - how easily E-fields form
#    mu_0      (magnetic permeability)  - how easily B-fields form
#
# These two properties FULLY DETERMINE the wave speed.
#
# FORMULA:
#   c = 1 / sqrt(epsilon_0 * mu_0)
#
# WHY THIS STRUCTURE?
#   An electromagnetic wave is a coupled oscillation of E and B fields.
#   The wave equation for such oscillations in free space is:
#
#      ∂²E/∂t² = (1/(epsilon_0 * mu_0)) * ∂²E/∂x²
#
#   This is a standard wave equation with speed v = 1/sqrt(epsilon_0 * mu_0).
#
# PHYSICAL INTERPRETATION:
#   epsilon_0 * mu_0 is the INERTIA of the electromagnetic field.
#   - epsilon_0: electric inertia (resistance to E-field changes)
#   - mu_0:      magnetic inertia (resistance to B-field changes)
#
#   The product epsilon_0 * mu_0 has units of (s/m)^2.
#   Taking the reciprocal square root gives m/s - a velocity!
#
#   Lower inertia → faster wave propagation
#   Higher inertia → slower wave propagation
#
# THIS IS THE MOST FUNDAMENTAL DERIVATION IN PHYSICS:
#   The speed of light is NOT an independent constant.
#   It EMERGES from the structure of the electromagnetic field.
#
# ============================================================

# DERIVED speed of light (direct from axioms)
c_derived = 1.0 / np.sqrt(epsilon_0 * mu_0)

# Also derive vacuum impedance (for context)
Z_0 = np.sqrt(mu_0 / epsilon_0)

# Error vs SI definition
error = (c_derived - c_SI) / c_SI * 100

# ============================================================
# OUTPUT
# ============================================================
if __name__ == "__main__":
    print("=" * 70)
    print("TRIPHASE V16 - DERIVATIVE 5: SPEED OF LIGHT (c)")
    print("Framework: WaveMechanics_Primitive")
    print("Tag: (D) DERIVED - Direct from the two axioms")
    print("=" * 70)

    print("\n" + "-" * 70)
    print("INPUTS: Three Observed Vacuum Properties ONLY")
    print("-" * 70)
    print(f"   epsilon_0 = {epsilon_0:.10e} F/m")
    print(f"   mu_0      = {mu_0:.10e} H/m")

    print("\n" + "-" * 70)
    print("DERIVATION: Wave Propagation in the Vector Frame")
    print("-" * 70)
    print(f"\n   Formula: c = 1 / sqrt(epsilon_0 * mu_0)")
    print(f"\n   Step 1: Calculate the product")
    print(f"      epsilon_0 * mu_0 = {epsilon_0:.6e} × {mu_0:.6e}")
    print(f"                       = {epsilon_0 * mu_0:.6e} s^2/m^2")
    print(f"\n   Step 2: Take the square root")
    print(f"      sqrt(epsilon_0 * mu_0) = {np.sqrt(epsilon_0 * mu_0):.6e} s/m")
    print(f"\n   Step 3: Take the reciprocal")
    print(f"      c = 1 / {np.sqrt(epsilon_0 * mu_0):.6e}")
    print(f"        = {c_derived:.0f} m/s")

    print("\n" + "-" * 70)
    print("ALSO DERIVED: Vacuum Impedance")
    print("-" * 70)
    print(f"\n   Z_0 = sqrt(mu_0 / epsilon_0)")
    print(f"       = sqrt({mu_0:.6e} / {epsilon_0:.6e})")
    print(f"       = {Z_0:.4f} ohms")
    print(f"\n   Note: Z_0 = mu_0 * c = 1 / (epsilon_0 * c)")
    print(f"         This is the impedance of free space to EM waves.")

    print("\n" + "-" * 70)
    print("CALIBRATION CHECK (SI definition - NOT used in calculation)")
    print("-" * 70)
    print(f"\n   Derived:     {c_derived:.0f} m/s")
    print(f"   SI 2019:     {c_SI} m/s (exact by definition)")
    print(f"\n   Error:       {error:.10f}%")
    print(f"\n   NOTE: The SI definition of the meter is BASED ON c.")
    print(f"   Since 1983, the meter is defined as the distance light")
    print(f"   travels in 1/299792458 seconds. This makes c EXACT.")

    print("\n" + "-" * 70)
    print("PHYSICAL MEANING")
    print("-" * 70)
    print(f"""
   c = {c_derived:.0f} m/s is the WAVE SPEED of the electromagnetic
   vector frame.

   This is NOT an arbitrary constant - it's DERIVED from the
   fundamental properties of the vacuum:
      epsilon_0: electric permittivity (field formation inertia)
      mu_0:      magnetic permeability (field flow inertia)

   The formula c = 1/sqrt(epsilon_0 * mu_0) shows that the speed
   of light is determined by the INERTIA of the electromagnetic field.

   DIMENSIONAL ANALYSIS:
      [epsilon_0] = F/m = C^2 s^2 / (kg m^3)
      [mu_0]      = H/m = kg m / C^2
      [epsilon_0 * mu_0] = (C^2 s^2 / (kg m^3)) × (kg m / C^2)
                         = s^2 / m^2
      [1/sqrt(epsilon_0 * mu_0)] = m/s  ✓

   Lower electromagnetic inertia → faster light propagation
   Higher electromagnetic inertia → slower light propagation

   In our universe: epsilon_0 * mu_0 = {epsilon_0 * mu_0:.3e} s^2/m^2
   This gives: c = {c_derived:.0f} m/s

   ALL electromagnetic phenomena (radio, light, X-rays, gamma rays)
   propagate at this speed in vacuum.

   This is the MAXIMUM INFORMATION TRANSFER RATE in the universe.
""")

    print("=" * 70)
    print(f"RESULT: c = {c_derived:.0f} m/s")
    print(f"ERROR:  {error:.10f}% vs SI 2019 definition")
    print("=" * 70)
    print()
    input("Press Enter to exit...")
