# -*- coding: utf-8 -*-
"""
============================================================
TRIPHASE V16 - INDIVIDUAL DERIVATIVE SCRIPT
============================================================
Derivative: Vector Frame Properties (6 fundamental constants)
Framework:  WaveMechanics_Primitive
Row:        6

INPUTS:     epsilon_0, mu_0 ONLY (Three Observed Vacuum Properties)
MEASURED:   CODATA 2022 values used as CALIBRATION CHECK ONLY

Tag: (D) DERIVED - All six from epsilon_0, mu_0

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
# CALIBRATION CHECKPOINT (measured - NOT used in calculation)
# ============================================================
# CODATA 2022 values for comparison
c_SI = 299792458              # m/s (exact by SI definition)
Z_0_CODATA = 376.730313668    # ohms

# ============================================================
# WAVE MECHANICS PRIMITIVE DERIVATION
# ============================================================
#
# MECHANISM:
# The "vector frame" is the electromagnetic vacuum structure that
# pervades all of space. It has SIX fundamental properties that
# can be arranged into three complementary pairs:
#
# PRIMARY PROPERTIES (directly observed):
#    1. epsilon_0: Electric permittivity (8.854e-12 F/m)
#    2. mu_0:      Magnetic permeability (1.257e-6 H/m)
#    3. Z_0:       Vacuum impedance (376.7 ohms)
#
# RECIPROCAL PROPERTIES (derived):
#    4. vf_0 = 1/epsilon_0: Imperivity (1.129e11 m/F)
#    5. rho_0 = 1/mu_0:     Reluctivity (7.958e5 m/H)
#    6. Y_0 = 1/Z_0:        Admittance (2.654e-3 siemens)
#
# RELATIONSHIPS:
#    Z_0 = sqrt(mu_0 / epsilon_0)     (impedance from permeabilities)
#    c = 1 / sqrt(epsilon_0 * mu_0)   (wave speed from inertia)
#    Z_0 = mu_0 * c = 1/(epsilon_0 * c)  (impedance-speed relations)
#
# THESE SIX PROPERTIES FULLY DEFINE THE VECTOR FRAME.
#
# All electromagnetic phenomena emerge from interactions with
# this six-dimensional parameter space.
#
# WHY SIX PROPERTIES?
#    The vector frame has THREE fundamental aspects:
#       - Electric field structure (epsilon_0, vf_0)
#       - Magnetic field structure (mu_0, rho_0)
#       - Wave coupling structure (Z_0, Y_0)
#
#    Each aspect has a PROPERTY and its RECIPROCAL.
#    3 aspects × 2 properties = 6 total properties
#
# DIMENSIONAL CLOSURE:
#    These six properties form a closed dimensional algebra:
#       [epsilon_0] = F/m     [vf_0] = m/F
#       [mu_0]      = H/m     [rho_0] = m/H
#       [Z_0]       = ohm     [Y_0] = siemens = 1/ohm
#
#    Every electromagnetic quantity can be expressed as products
#    and ratios of these six properties.
#
# ============================================================

# DERIVED PRIMARY PROPERTIES
Z_0 = np.sqrt(mu_0 / epsilon_0)              # vacuum impedance
c   = 1.0 / np.sqrt(epsilon_0 * mu_0)        # speed of light

# DERIVED RECIPROCAL PROPERTIES
vf_0  = 1.0 / epsilon_0   # imperivity (reciprocal permittivity)
rho_0 = 1.0 / mu_0        # reluctivity (reciprocal permeability)
Y_0   = 1.0 / Z_0         # admittance (reciprocal impedance)

# Verify closure relationships
check_Z0_from_c = mu_0 * c
check_Z0_from_epsilon = 1.0 / (epsilon_0 * c)

# Errors vs calibration
error_c = (c - c_SI) / c_SI * 100
error_Z0 = (Z_0 - Z_0_CODATA) / Z_0_CODATA * 100

# ============================================================
# OUTPUT
# ============================================================
if __name__ == "__main__":
    print("=" * 70)
    print("TRIPHASE V16 - DERIVATIVE 6: VECTOR FRAME PROPERTIES")
    print("Framework: WaveMechanics_Primitive")
    print("Tag: (D) DERIVED - All six from epsilon_0, mu_0")
    print("=" * 70)

    print("\n" + "-" * 70)
    print("INPUTS: Three Observed Vacuum Properties")
    print("-" * 70)
    print(f"   1. epsilon_0 = {epsilon_0:.10e} F/m  (electric permittivity)")
    print(f"   2. mu_0      = {mu_0:.10e} H/m  (magnetic permeability)")

    print("\n" + "-" * 70)
    print("DERIVED PRIMARY PROPERTIES")
    print("-" * 70)
    print(f"\n   3. Z_0 = sqrt(mu_0 / epsilon_0)")
    print(f"        = sqrt({mu_0:.6e} / {epsilon_0:.6e})")
    print(f"        = {Z_0:.9f} ohms  (vacuum impedance)")
    print(f"\n   Also derived: c = 1/sqrt(epsilon_0 * mu_0)")
    print(f"                   = {c:.0f} m/s  (speed of light)")

    print("\n" + "-" * 70)
    print("DERIVED RECIPROCAL PROPERTIES")
    print("-" * 70)
    print(f"\n   4. vf_0 = 1/epsilon_0")
    print(f"         = 1/{epsilon_0:.6e}")
    print(f"         = {vf_0:.6e} m/F  (imperivity)")
    print(f"\n   5. rho_0 = 1/mu_0")
    print(f"          = 1/{mu_0:.6e}")
    print(f"          = {rho_0:.6e} m/H  (reluctivity)")
    print(f"\n   6. Y_0 = 1/Z_0")
    print(f"        = 1/{Z_0:.9f}")
    print(f"        = {Y_0:.10f} S  (admittance)")

    print("\n" + "-" * 70)
    print("CLOSURE RELATIONSHIPS (verification)")
    print("-" * 70)
    print(f"\n   Z_0 = mu_0 * c")
    print(f"       = {mu_0:.6e} × {c:.6e}")
    print(f"       = {check_Z0_from_c:.9f} ohms  ✓")
    print(f"\n   Z_0 = 1 / (epsilon_0 * c)")
    print(f"       = 1 / ({epsilon_0:.6e} × {c:.6e})")
    print(f"       = {check_Z0_from_epsilon:.9f} ohms  ✓")
    print(f"\n   epsilon_0 * vf_0 = {epsilon_0:.6e} × {vf_0:.6e} = {epsilon_0 * vf_0:.1f}  ✓")
    print(f"   mu_0 * rho_0     = {mu_0:.6e} × {rho_0:.6e} = {mu_0 * rho_0:.1f}  ✓")
    print(f"   Z_0 * Y_0        = {Z_0:.6f} × {Y_0:.6f} = {Z_0 * Y_0:.1f}  ✓")

    print("\n" + "-" * 70)
    print("CALIBRATION CHECK (measured values - NOT used in calculation)")
    print("-" * 70)
    print(f"\n   Derived c:   {c:.0f} m/s")
    print(f"   SI 2019:     {c_SI} m/s (exact)")
    print(f"   Error:       {error_c:.10f}%")
    print(f"\n   Derived Z_0: {Z_0:.9f} ohms")
    print(f"   CODATA 2022: {Z_0_CODATA:.9f} ohms")
    print(f"   Error:       {error_Z0:.10f}%")

    print("\n" + "-" * 70)
    print("PHYSICAL MEANING")
    print("-" * 70)
    print(f"""
   The VECTOR FRAME is the electromagnetic vacuum structure.
   It is defined by SIX fundamental properties:

   THREE PRIMARY PROPERTIES (directly observed):
      epsilon_0 = {epsilon_0:.3e} F/m  (how easily E-fields form)
      mu_0      = {mu_0:.3e} H/m  (how easily B-fields form)
      Z_0       = {Z_0:.1f} ohms      (wave coupling strength)

   THREE RECIPROCAL PROPERTIES (derived):
      vf_0  = {vf_0:.3e} m/F  (electric field resistance)
      rho_0 = {rho_0:.3e} m/H  (magnetic field resistance)
      Y_0   = {Y_0:.3e} S    (wave coupling conductance)

   These six properties form a COMPLETE BASIS for all electromagnetic
   phenomena. Every EM quantity can be expressed as products and
   ratios of these six.

   WHY THREE PAIRS?
      The vacuum has three fundamental aspects:
         1. Electric structure (epsilon_0 ↔ vf_0)
         2. Magnetic structure (mu_0 ↔ rho_0)
         3. Wave coupling (Z_0 ↔ Y_0)

      Each aspect has a property and its reciprocal.
      Together they form a 6-dimensional parameter space.

   RELATIONSHIPS:
      Z_0 = sqrt(mu_0/epsilon_0) = mu_0*c = 1/(epsilon_0*c)
      c = 1/sqrt(epsilon_0*mu_0)

   These closure relations show the six properties are NOT
   independent - they're interconnected through the geometry
   of the electromagnetic field.

   ALL OF PHYSICS emerges from these six numbers.
""")

    print("\n" + "-" * 70)
    print("SUMMARY TABLE")
    print("-" * 70)
    print("\n   PROPERTY        VALUE              UNITS    NAME")
    print("   " + "-" * 62)
    print(f"   epsilon_0    {epsilon_0:12.6e}   F/m      Electric permittivity")
    print(f"   mu_0         {mu_0:12.6e}   H/m      Magnetic permeability")
    print(f"   Z_0          {Z_0:12.6f}   ohm      Vacuum impedance")
    print(f"   vf_0         {vf_0:12.6e}   m/F      Imperivity (1/epsilon_0)")
    print(f"   rho_0        {rho_0:12.6e}   m/H      Reluctivity (1/mu_0)")
    print(f"   Y_0          {Y_0:12.6e}   S        Admittance (1/Z_0)")
    print(f"   c            {c:12.0f}   m/s      Speed of light")

    print("\n" + "=" * 70)
    print("RESULT: Vector frame fully characterized by 6 properties")
    print("        All derived from epsilon_0 and mu_0")
    print("=" * 70)
    print()
    input("Press Enter to exit...")
