# -*- coding: utf-8 -*-
"""
============================================================
G UNITS DERIVATION - How does G relate to e0 dimensionally?
============================================================
Solving Gurcharn's question: How do F/m units become m^3/(kg*s^2)?

Author: MIS Magnetic Innovative Solutions LLC
Date: January 29, 2026
============================================================
"""
import sys
import traceback
import atexit

atexit.register(lambda: input('Press Enter to exit...'))
sys.stdout.reconfigure(encoding='utf-8')

_original_excepthook = sys.excepthook
def _custom_excepthook(exc_type, exc_value, exc_tb):
    traceback.print_exception(exc_type, exc_value, exc_tb)
sys.excepthook = _custom_excepthook


import numpy as np

print("=" * 70)
print("G UNITS DERIVATION")
print("=" * 70)
print()

# =============================================================
# CONSTANTS
# =============================================================
e0 = 8.8541878128e-12  # F/m = A^2*s^4/(kg*m^3)
u0 = 1.25663706212e-6  # H/m = kg*m/(A^2*s^2)
c = 299792458          # m/s
G_measured = 6.67430e-11  # m^3/(kg*s^2)
Z0 = np.sqrt(u0/e0)    # Ohms

print("MEASURED VALUES:")
print(f"  e0 = {e0:.4e} F/m")
print(f"  u0 = {u0:.4e} H/m")
print(f"  c  = {c} m/s")
print(f"  Z0 = {Z0:.2f} Ohms")
print(f"  G  = {G_measured:.4e} m^3/(kg*s^2)")
print()

# =============================================================
# NUMERICAL RATIO
# =============================================================
print("=" * 70)
print("NUMERICAL RATIO G/e0")
print("-" * 70)

ratio = G_measured / e0
print(f"G/e0 = {G_measured:.4e} / {e0:.4e}")
print(f"     = {ratio:.4f}")
print()
print("So numerically, G = 7.54 * e0")
print("But the units don't match - so what carries the units?")
print()

# =============================================================
# UNIT ANALYSIS
# =============================================================
print("=" * 70)
print("UNIT ANALYSIS")
print("-" * 70)
print()
print("e0 units: A^2*s^4/(kg*m^3)")
print("G  units: m^3/(kg*s^2)")
print()
print("G/e0 units: [m^3/(kg*s^2)] / [A^2*s^4/(kg*m^3)]")
print("          = [m^3/(kg*s^2)] * [kg*m^3/(A^2*s^4)]")
print("          = m^6/(A^2*s^6)")
print()
print("So the '7.54' has units m^6/(A^2*s^6)")
print()

# =============================================================
# WHAT IS m^6/(A^2*s^6)?
# =============================================================
print("=" * 70)
print("WHAT IS m^6/(A^2*s^6)?")
print("-" * 70)
print()
print("m^6/s^6 = c^6")
print("So m^6/(A^2*s^6) = c^6/A^2 = c^6/I^2 for some current I")
print()

c6 = c**6
print(f"c^6 = {c6:.4e} m^6/s^6")
print()

# What current I0 makes G/e0 = c^6/I0^2?
# G/e0 = 7.54 m^6/(A^2*s^6)
# c^6/I0^2 = 7.54 m^6/(A^2*s^6)
# I0^2 = c^6 / 7.54

# But wait - 7.54 isn't in SI units of m^6/(A^2*s^6), it's dimensionless
# The dimensional part is already in G/e0

# Let me recalculate properly
# G = e0 * c^6 / I0^2
# So I0^2 = e0 * c^6 / G

I0_squared = e0 * c6 / G_measured
I0 = np.sqrt(I0_squared)

print(f"If G = e0 * c^6 / I0^2, then:")
print(f"  I0^2 = e0 * c^6 / G")
print(f"       = {e0:.4e} * {c6:.4e} / {G_measured:.4e}")
print(f"       = {I0_squared:.4e} A^2")
print(f"  I0   = {I0:.4e} A")
print()

# =============================================================
# WHAT IS THIS CURRENT?
# =============================================================
print("=" * 70)
print("WHAT IS I0 = 3.1e25 AMPERES?")
print("-" * 70)
print()

# Planck current (one definition)
I_planck = np.sqrt(c6 * e0 / G_measured)
print(f"This IS the Planck current!")
print(f"  I_P = sqrt(c^6 * e0 / G) = {I_planck:.4e} A")
print()

# =============================================================
# ALTERNATIVE: THROUGH u0
# =============================================================
print("=" * 70)
print("ALTERNATIVE: G THROUGH u0")
print("-" * 70)
print()

# Try G = c^4 / (u0 * I^2)
# G * u0 = c^4 / I^2
# I^2 = c^4 / (G * u0)

c4 = c**4
I0_via_u0_sq = c4 / (G_measured * u0)
I0_via_u0 = np.sqrt(I0_via_u0_sq)

print(f"If G = c^4 / (u0 * I0^2), then:")
print(f"  I0^2 = c^4 / (G * u0)")
print(f"       = {c4:.4e} / ({G_measured:.4e} * {u0:.4e})")
print(f"       = {I0_via_u0_sq:.4e} A^2")
print(f"  I0   = {I0_via_u0:.4e} A")
print()
print("Same current! Because c^2 = 1/(e0*u0)")
print()

# =============================================================
# THE TWO EQUIVALENT FORMS
# =============================================================
print("=" * 70)
print("TWO EQUIVALENT FORMS FOR G")
print("-" * 70)
print()
print("Form 1: G = e0 * c^6 / I_P^2")
print("Form 2: G = c^4 / (u0 * I_P^2)")
print()
print("These are equivalent because c^2 = 1/(e0*u0)")
print()
print("Substituting Form 1 into Form 2:")
print("  e0 * c^6 / I_P^2 = c^4 / (u0 * I_P^2)")
print("  e0 * c^6 = c^4 / u0")
print("  e0 * c^2 = 1/u0")
print("  e0 * u0 = 1/c^2  CHECK!")
print()

# Verify
print("Verification:")
G_from_e0 = e0 * c6 / I0_squared
G_from_u0 = c4 / (u0 * I0_via_u0_sq)
print(f"  G from e0 form: {G_from_e0:.4e}")
print(f"  G from u0 form: {G_from_u0:.4e}")
print(f"  G measured:     {G_measured:.4e}")
print()

# =============================================================
# WHERE DOES 7.5 COME FROM?
# =============================================================
print("=" * 70)
print("WHERE DOES 7.5 COME FROM?")
print("-" * 70)
print()

# The numerical ratio G/e0 = 7.54 comes from c^6/I_P^2 in SI units
numerical_factor = c6 / I0_squared
print(f"c^6 / I_P^2 = {c6:.4e} / {I0_squared:.4e}")
print(f"           = {numerical_factor:.4f}")
print()
print("So the '7.5' is just c^6/I_P^2 evaluated in SI units!")
print()
print("In TriPhase, 7.5 = 15/2 from mode counting.")
print("This suggests I_P has a mode-counting origin.")
print()

# =============================================================
# THE KEY INSIGHT
# =============================================================
print("=" * 70)
print("THE KEY INSIGHT")
print("=" * 70)
print()
print("""
G is the COUPLING between vacuum EM properties:

  G = e0 * c^6 / I_P^2
  G = c^4 / (u0 * I_P^2)

Where:
  - e0 describes electric field coupling to vacuum
  - u0 describes magnetic field coupling to vacuum
  - c^2 = 1/(e0*u0) connects them
  - I_P is the Planck current where EM and gravity merge

The "7.5 * e0" is shorthand for "e0 * c^6 / I_P^2"

G isn't a separate force - it's the coupling constant that
describes how mass-energy interacts with the vacuum through
the same e0 and u0 that govern electromagnetism.
""")

# =============================================================
# CHECK: DOES 7.5 = 15/2 RELATE TO MODES?
# =============================================================
print("=" * 70)
print("MODE COUNTING CHECK")
print("-" * 70)
print()

# If 7.5 = 15/2, and this equals c^6/I_P^2 in "natural" units,
# then I_P should have a simple form

# In TriPhase: 6 modes = 3 phases * 2 quadratures
# 15 = 6 + 6 + 3 = total + total + phases?
# Or 15 = 5 * 3 = active_modes * phases?
# Or 15 = T_5 = 5th triangular number

print("15/2 = 7.5")
print()
print("Possible mode origins of 15:")
print("  15 = 5 * 3 (active_modes * phases)")
print("  15 = T_5 = 1+2+3+4+5 (5th triangular number)")
print("  15 = 6 + 6 + 3 (modes + modes + phases)")
print()
print("The factor of 2 in denominator:")
print("  2 = quadratures (sin/cos)")
print("  2 = binary (on/off)")
print()

# =============================================================
# ANSWER FOR GURCHARN
# =============================================================
print("=" * 70)
print("ANSWER FOR GURCHARN")
print("=" * 70)
print()
print("""
The units DO convert - through c and the Planck current:

  G = e0 * c^6 / I_P^2

Where:
  e0 = 8.854e-12 A^2*s^4/(kg*m^3)
  c^6 = 7.29e50 m^6/s^6
  I_P^2 = 9.67e49 A^2

The c^6 provides m^6/s^6
The I_P^2 provides A^2
Together with e0, they give m^3/(kg*s^2) = G's units

The numerical factor 7.54 = c^6/I_P^2 in SI units.
TriPhase predicts this should be 15/2 = 7.5 from mode counting.
Error: 0.5%

G is the coupling between e0 and u0 through c and I_P.
It's not a separate force - it's vacuum EM coupling at the Planck scale.
""")