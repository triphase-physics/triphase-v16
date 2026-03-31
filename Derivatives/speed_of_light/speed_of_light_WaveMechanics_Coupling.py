# -*- coding: utf-8 -*-
"""
============================================================
TRIPHASE V16 - INDIVIDUAL DERIVATIVE SCRIPT
============================================================
Derivative: Speed of Light (c = 299792458 m/s)
Framework:  WaveMechanics_Coupling
Row:        5

INPUTS:     epsilon_0, mu_0 ONLY (Three Observed Vacuum Properties)
MEASURED:   SI 2019 exact definition - NOT used in calculation

Tag: (D) DERIVED - Pure epsilon_0, mu_0 derivation

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

# ============================================================
# CALIBRATION CHECKPOINT (SI exact - NOT used in calculation)
# ============================================================
c_SI = 299792458.0    # m/s (SI 2019 exact definition)

# ============================================================
# WAVE MECHANICS COUPLING DERIVATION
# ============================================================
#
# COUPLING CHAIN:
# The speed of light is the COUPLING VELOCITY - the maximum
# rate at which signals can propagate through the vector frame.
#
# c is NOT arbitrary - it's determined by the IMPEDANCE STRUCTURE
# of the electromagnetic vector frame.
#
# FUNDAMENTAL FORMULA:
#   c = 1 / sqrt(epsilon_0 * mu_0)
#
# This is the most fundamental equation in physics. It says:
# The propagation velocity is determined by the IMPEDANCE
# PROPERTIES (epsilon_0, mu_0) of the medium.
#
# COUPLING INTERPRETATION:
#
#   epsilon_0: Electric field storage capacity (F/m)
#   mu_0:      Magnetic field storage capacity (H/m)
#   c:         Maximum coupling velocity through this medium
#
# WHY c IS THE MAXIMUM VELOCITY:
#
#   When a disturbance propagates, it must couple energy between
#   electric and magnetic field modes. The rate is limited by
#   how fast the frame can store and release energy.
#
#   epsilon_0 * mu_0 determines the "stiffness" of the frame.
#   c = 1/sqrt(epsilon_0 * mu_0) is the wave velocity in this
#   "stiff" medium.
#
# IMPEDANCE RELATIONS:
#
#   Z_0 = sqrt(mu_0 / epsilon_0) = 376.73 ohms
#   c = 1 / (Z_0 * epsilon_0) = Z_0 / mu_0
#
# These relations show c is the coupling velocity that matches
# the impedance Z_0.
#
# DIMENSIONAL ANALYSIS:
#   [epsilon_0] = F/m = C^2 s^2 / (kg m^3)
#   [mu_0]      = H/m = kg m / C^2
#
#   [epsilon_0 * mu_0] = [C^2 s^2 / (kg m^3)] * [kg m / C^2]
#                      = s^2 / m^2
#
#   [1 / sqrt(epsilon_0 * mu_0)] = [m^2 / s^2]^(1/2)
#                                 = m/s  ✓
#
# COUPLING HIERARCHY:
#   c is the REFERENCE VELOCITY for all coupling strengths:
#   - EM coupling: alpha = e^2 * Z_0 / (4*pi*hbar) at velocity c
#   - Gravity coupling: G = c^4 * 7.5 * epsilon_0^3 * mu_0^2
#   - All interactions propagate at maximum velocity c
#
# ============================================================

# Alternative impedance forms
c_from_impedance_1 = 1.0 / (Z_0 * epsilon_0)   # c = 1 / (Z_0 * epsilon_0)
c_from_impedance_2 = Z_0 / mu_0                 # c = Z_0 / mu_0

# Verify consistency
error_1 = (c_from_impedance_1 - c) / c * 100
error_2 = (c_from_impedance_2 - c) / c * 100

# Error vs SI definition
error_SI = (c - c_SI) / c_SI * 100

# ============================================================
# OUTPUT
# ============================================================
if __name__ == "__main__":
    print("=" * 70)
    print("TRIPHASE V16 - DERIVATIVE 5: SPEED OF LIGHT (c)")
    print("Framework: WaveMechanics_Coupling")
    print("Tag: (D) DERIVED - Pure epsilon_0, mu_0 derivation")
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
    print("COUPLING CHAIN: Maximum Signal Velocity")
    print("-" * 70)
    print(f"""
   STEP 1: Vector frame impedance properties
      epsilon_0 = {epsilon_0:.6e} F/m (electric storage)
      mu_0      = {mu_0:.6e} H/m (magnetic storage)

   STEP 2: Frame "stiffness"
      epsilon_0 * mu_0 = {epsilon_0 * mu_0:.6e} s^2/m^2
      This determines how fast disturbances propagate

   STEP 3: Coupling velocity
      c = 1 / sqrt(epsilon_0 * mu_0)
        = 1 / sqrt({epsilon_0 * mu_0:.6e})
        = {c:.6f} m/s

   STEP 4: Impedance matching
      Z_0 = sqrt(mu_0 / epsilon_0) = {Z_0:.6f} ohms
      c = 1 / (Z_0 * epsilon_0) = {c_from_impedance_1:.6f} m/s
      c = Z_0 / mu_0 = {c_from_impedance_2:.6f} m/s

   COUPLING INTERPRETATION:
      c is the MAXIMUM VELOCITY at which signals can couple
      energy between electric and magnetic field modes.

      The vector frame has impedance Z_0 = {Z_0:.2f} ohms.
      Waves propagating at velocity c match this impedance.

      WHY c IS CONSTANT:
      epsilon_0 and mu_0 are properties of the vector frame.
      They don't change with observer motion.
      Therefore c = 1/sqrt(epsilon_0 * mu_0) is the SAME
      for all observers - this is special relativity.
""")

    print("\n" + "-" * 70)
    print("DERIVATION: Three Equivalent Forms")
    print("-" * 70)
    print(f"\n   Form 1: Direct from epsilon_0, mu_0")
    print(f"      c = 1 / sqrt(epsilon_0 * mu_0)")
    print(f"        = 1 / sqrt({epsilon_0:.6e} × {mu_0:.6e})")
    print(f"        = {c:.6f} m/s")
    print(f"\n   Form 2: From impedance (electric)")
    print(f"      c = 1 / (Z_0 * epsilon_0)")
    print(f"        = 1 / ({Z_0:.6f} × {epsilon_0:.6e})")
    print(f"        = {c_from_impedance_1:.6f} m/s")
    print(f"      Error: {error_1:.2e}%")
    print(f"\n   Form 3: From impedance (magnetic)")
    print(f"      c = Z_0 / mu_0")
    print(f"        = {Z_0:.6f} / {mu_0:.6e}")
    print(f"        = {c_from_impedance_2:.6f} m/s")
    print(f"      Error: {error_2:.2e}%")

    print("\n" + "-" * 70)
    print("CALIBRATION CHECK (SI exact - NOT used in calculation)")
    print("-" * 70)
    print(f"\n   Derived:      {c:.6f} m/s")
    print(f"   SI 2019:      {c_SI:.6f} m/s (exact definition)")
    print(f"\n   Error:        {error_SI:+.2e}%")

    print("\n" + "-" * 70)
    print("PHYSICAL MEANING")
    print("-" * 70)
    print(f"""
   c = {c:.0f} m/s is the COUPLING VELOCITY - the maximum
   rate at which signals can propagate through spacetime.

   FUNDAMENTAL SIGNIFICANCE:

   1. MAXIMUM VELOCITY
      No signal can travel faster than c.
      Why? Because c is determined by the IMPEDANCE of the
      vector frame itself (epsilon_0, mu_0).
      Going faster would require infinite energy to overcome
      frame impedance.

   2. IMPEDANCE MATCHING
      Waves traveling at c perfectly match the frame impedance
      Z_0 = {Z_0:.2f} ohms. This is why light propagates without
      dispersion - it's impedance-matched to the medium.

   3. RELATIVITY
      c is the same for all observers because epsilon_0 and mu_0
      are properties of the vector frame, not the observer.
      This is the foundation of special relativity.

   4. COUPLING HIERARCHY
      All interactions propagate at velocity c:
      - EM coupling: alpha at velocity c
      - Gravity coupling: G uses c^4 scaling
      - Strong/weak: also propagate at c (massless bosons)

   DIMENSIONAL STRUCTURE:
      [c] = m/s (velocity)
      [epsilon_0 * mu_0] = s^2/m^2 (inverse velocity squared)
      [1/sqrt(epsilon_0 * mu_0)] = m/s ✓

   The speed of light is NOT a "speed limit" arbitrarily imposed.
   It's the NATURAL VELOCITY of wave propagation in a medium
   with impedance properties (epsilon_0, mu_0).

   Just as sound has a maximum speed in air (set by air's
   compressibility and density), light has a maximum speed in
   spacetime (set by spacetime's epsilon_0 and mu_0).

   The difference: spacetime is UNIVERSAL, so c is universal.
""")

    print("=" * 70)
    print(f"RESULT: c = {c:.6f} m/s")
    print(f"             = {c:.6e} m/s")
    print(f"ERROR:  {error_SI:+.2e}% vs SI 2019 exact definition")
    print("=" * 70)
    print()
    input("Press Enter to exit...")
