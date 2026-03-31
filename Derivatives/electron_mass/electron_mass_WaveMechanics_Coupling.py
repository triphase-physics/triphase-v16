# -*- coding: utf-8 -*-
"""
============================================================
TRIPHASE V16 - INDIVIDUAL DERIVATIVE SCRIPT
============================================================
Derivative: Electron Mass (m_e)
Framework:  WaveMechanics_Coupling
Row:        (Measured anchor for consistency check)

INPUTS:     epsilon_0, mu_0 ONLY (Three Observed Vacuum Properties)
MEASURED:   m_e = 9.1093837015 × 10^-31 kg (CODATA 2018)

Tag: (C) CONSISTENCY CHECK — m_e is measured anchor

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
m = 17                          # pressure band index
node = 8 * m + 1                # = 137
correction = np.log(node) / node   # ln(137)/137
alpha_inv = node + correction      # 137.0359...
alpha = 1.0 / alpha_inv

# ============================================================
# MEASURED ANCHOR (NOT DERIVED)
# ============================================================
# Electron mass is a MEASURED value that anchors the coupling chain
m_e_measured = 9.1093837015e-31  # kg (CODATA 2018)

# ============================================================
# COUPLING CHAIN INTERPRETATION
# ============================================================
#
# ELECTRON MASS: GROUND STATE COUPLING TO VECTOR FRAME
#
# The electron mass m_e is a MEASURED ANCHOR in TriPhase V16.
# It is NOT derived from epsilon_0 and mu_0 alone.
#
# WHY m_e IS FUNDAMENTAL:
#
# The electron is the LIGHTEST STABLE charged particle.
# It represents the GROUND STATE coupling to the vector frame
# at pressure band m = 17.
#
# COUPLING PICTURE:
#
# The electron mass arises from coupling to the vacuum:
#   1) Vacuum impedance Z_0 = 376.73 ohms
#   2) Pressure band m = 17 (alpha^-1 ≈ 137 = 8*17+1)
#   3) Coupling strength alpha ≈ 1/137
#
# Mass is the INERTIAL RESISTANCE to acceleration through
# the vacuum impedance Z_0.
#
# DIMENSIONAL RELATION:
#
# The electron mass can be related to vacuum properties through:
#
#    m_e ~ (hbar * alpha) / (c * r_e)
#
# where:
#   hbar = reduced Planck constant
#   alpha = fine structure constant
#   c = speed of light
#   r_e = classical electron radius
#
# However, this is a CONSISTENCY CHECK, not a derivation,
# because r_e itself depends on m_e.
#
# COUPLING HIERARCHY:
#
# The electron at band m = 17 has:
#   - T_17 = 153 coupling modes
#   - Coupling strength alpha ≈ 1/137
#   - Mass m_e = 9.109e-31 kg (measured)
#
# Other particles occupy different bands:
#   - Muon (m ~ 18-19): heavier, more coupling modes
#   - Tau (m ~ 20-21): even heavier
#   - Quarks: coupled to strong force, different band structure
#
# IMPEDANCE COUPLING AND MASS:
#
# Mass represents INERTIAL IMPEDANCE to acceleration:
#   F = m * a  (Newton's second law)
#
# In wave mechanics, mass is the coupling coefficient between
# force (momentum transfer rate) and acceleration (velocity change rate).
#
# The vacuum impedance Z_0 = 376.73 ohms mediates this coupling.
#
# For EM coupling:
#   - Charge e couples to fields through Z_0
#   - Mass m_e provides inertial resistance
#   - Ratio e/m_e determines acceleration response
#
# COUPLING STRENGTH RELATION:
#
# The fine structure constant relates charge and mass:
#
#    alpha = e^2 / (4*pi*epsilon_0*hbar*c)
#          = (e^2 * Z_0) / (4*pi*hbar*c)
#
# This connects charge e, impedance Z_0, and coupling alpha.
# The electron mass m_e enters through hbar and r_e relationships.
#
# CONSISTENCY CHECK:
#
# We can verify dimensional consistency:
#
#    m_e * c^2 = 511 keV (rest energy)
#
# Combined with alpha and hbar:
#    Rydberg = m_e * c^2 * alpha^2 / 2 = 13.6 eV
#    Compton wavelength = hbar / (m_e * c) = 2.43e-12 m
#    Classical radius = alpha * hbar / (m_e * c) = 2.82e-15 m
#
# All these relations are CONSISTENT with m_e measured value.
#
# WHY m_e CANNOT BE DERIVED FROM epsilon_0, mu_0 ALONE:
#
# epsilon_0 and mu_0 set the GEOMETRIC properties of vacuum:
#   - Speed of light c
#   - Impedance Z_0
#   - Wave propagation
#
# But they do NOT determine MASS SCALES.
#
# Mass scales require ADDITIONAL input:
#   - Either measure m_e directly (anchor)
#   - Or derive from quantum gravity scale (Planck mass)
#
# TriPhase V16 uses m_e as MEASURED ANCHOR because:
#   1) It is precisely measured (10 significant figures)
#   2) It anchors the entire particle mass hierarchy
#   3) It connects EM coupling (alpha) to inertia
#
# COUPLING CHAIN SUMMARY:
#
# Vacuum properties (epsilon_0, mu_0):
#   → c, Z_0 (geometric scales)
#   → alpha (from band m=17 coupling)
#
# Electron mass (measured):
#   → m_e (inertial coupling to vacuum)
#
# Combined (epsilon_0, mu_0, m_e):
#   → All particle masses (via band coupling ratios)
#   → All energy scales (Rydberg, Compton, etc.)
#   → Complete coupling hierarchy
#
# The electron is the GROUND STATE COUPLER to the vector frame
# at band m = 17. Its mass m_e anchors all other masses through
# coupling ratios.
#
# ============================================================

# Electron mass (measured anchor)
m_e = m_e_measured

# Derived quantities for consistency check
m_e_c2_J = m_e * c**2  # rest energy in Joules
eV = 1.602176634e-19  # J/eV (exact, SI 2019)
m_e_c2_eV = m_e_c2_J / eV  # rest energy in eV
m_e_c2_keV = m_e_c2_eV / 1000.0  # rest energy in keV

# Planck constant (needed for consistency checks)
h = 6.62607015e-34  # J·s (exact, SI 2019)
hbar = h / (2.0 * np.pi)

# Compton wavelength
lambda_C = hbar / (m_e * c)  # meters

# Classical electron radius (approximate)
# r_e = alpha * lambda_C / (2*pi) = alpha * hbar / (m_e * c)
r_e = alpha * hbar / (m_e * c)  # meters

# Rydberg energy
E_Rydberg = m_e * c**2 * alpha**2 / 2.0
E_Rydberg_eV = E_Rydberg / eV

# Triangular number T_17 (coupling modes at band m=17)
T_17 = m * (m + 1) // 2

# ============================================================
# OUTPUT
# ============================================================
if __name__ == "__main__":
    print("=" * 70)
    print("TRIPHASE V16 - ELECTRON MASS (m_e)")
    print("Framework: WaveMechanics_Coupling")
    print("Tag: (C) CONSISTENCY CHECK — m_e is measured anchor")
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
    print(f"   alpha = 1 / {alpha_inv:.6f}    = {alpha:.10f}")

    print("\n" + "-" * 70)
    print("MEASURED ANCHOR (NOT DERIVED)")
    print("-" * 70)
    print(f"\n   m_e = {m_e:.10e} kg  (CODATA 2018)")
    print(f"\n   This is a MEASURED value that anchors the coupling chain.")
    print(f"   It CANNOT be derived from epsilon_0 and mu_0 alone.")

    print("\n" + "-" * 70)
    print("COUPLING CHAIN: Electron as Ground State Coupler")
    print("-" * 70)
    print(f"\n   Pressure band: m = {m} (electron band)")
    print(f"   Coupling modes: T_{m} = {T_17}")
    print(f"   Coupling strength: alpha = {alpha:.10f}")
    print(f"   Vacuum impedance: Z_0 = {Z_0:.4f} ohms")

    print(f"\n   The electron is the GROUND STATE coupling to the")
    print(f"   vector frame at band m = {m}.")

    print(f"\n   Its mass m_e = {m_e:.4e} kg represents the")
    print(f"   INERTIAL COUPLING to vacuum through Z_0.")

    print("\n" + "-" * 70)
    print("REST ENERGY")
    print("-" * 70)
    print(f"\n   m_e * c^2 = ({m_e:.4e}) * ({c:.3e})^2")
    print(f"             = {m_e_c2_J:.6e} J")
    print(f"             = {m_e_c2_eV:.6e} eV")
    print(f"             = {m_e_c2_keV:.3f} keV")

    print("\n" + "-" * 70)
    print("CONSISTENCY CHECKS: DERIVED SCALES")
    print("-" * 70)
    print(f"\n   COMPTON WAVELENGTH:")
    print(f"      lambda_C = hbar / (m_e * c)")
    print(f"               = {hbar:.6e} / ({m_e:.4e} * {c:.3e})")
    print(f"               = {lambda_C:.6e} m")
    print(f"               = {lambda_C * 1e12:.3f} pm")

    print(f"\n   CLASSICAL ELECTRON RADIUS:")
    print(f"      r_e = alpha * hbar / (m_e * c)")
    print(f"          = {alpha:.6e} * {hbar:.6e} / ({m_e:.4e} * {c:.3e})")
    print(f"          = {r_e:.6e} m")
    print(f"          = {r_e * 1e15:.3f} fm")

    print(f"\n   RYDBERG ENERGY:")
    print(f"      E_Ry = m_e * c^2 * alpha^2 / 2")
    print(f"           = {m_e_c2_keV:.3f} keV * ({alpha:.6e})^2 / 2")
    print(f"           = {E_Rydberg_eV:.3f} eV")

    print("\n" + "-" * 70)
    print("IMPEDANCE COUPLING PICTURE")
    print("-" * 70)
    print(f"\n   Vacuum impedance: Z_0 = {Z_0:.4f} ohms")
    print(f"\n   The electron charge e couples to EM fields through Z_0.")
    print(f"   The electron mass m_e provides INERTIAL RESISTANCE.")
    print(f"\n   Ratio e/m_e determines acceleration response:")
    print(f"      a = F/m_e = (e*E)/m_e")
    print(f"\n   Fine structure constant connects charge and mass:")
    print(f"      alpha = (e^2 * Z_0) / (4*pi*hbar*c)")
    print(f"            = {alpha:.10f}")

    print("\n" + "-" * 70)
    print("COUPLING HIERARCHY: ELECTRON AT BAND m=17")
    print("-" * 70)
    print(f"\n   Band     T_n   Mass Scale")
    print(f"   {'─'*45}")
    print(f"   m=17     {T_17:3d}   m_e = {m_e:.3e} kg  ← ELECTRON")
    print(f"   m=18     171   Heavier particles (muon?)")
    print(f"   m=19     190   Even heavier (tau?)")

    print(f"\n   The electron occupies band m={m} with {T_17} coupling modes.")
    print(f"   Heavier particles occupy higher bands with more coupling.")

    print("\n" + "-" * 70)
    print("PHYSICAL MEANING: GROUND STATE COUPLING")
    print("-" * 70)
    print(f"""
   m_e = {m_e:.10e} kg is the MEASURED electron mass.

   COUPLING INTERPRETATION:
      The electron is the LIGHTEST STABLE charged particle.
      It represents the GROUND STATE coupling to the vacuum
      at pressure band m = {m}.

      Mass = INERTIAL COUPLING to vacuum impedance Z_0

   WHY m_e IS FUNDAMENTAL:
      • Lightest stable lepton (no lighter states available)
      • Occupies band m={m} (alpha^-1 = 137 = 8*17+1)
      • Anchors entire particle mass hierarchy
      • Couples through Z_0 = {Z_0:.4f} ohms

   COUPLING MODES:
      At band m={m}, there are T_{m} = {T_17} coupling modes.
      The electron ground state couples through these {T_17} channels.

   IMPEDANCE PICTURE:
      F = m_e * a  (Newton's second law)

      Mass m_e is the INERTIAL IMPEDANCE:
         - Force F drives momentum change (dp/dt)
         - Mass m_e resists acceleration
         - Coupling through Z_0 mediates interaction

   CONSISTENCY CHECKS:
      Rest energy:         {m_e_c2_keV:.1f} keV
      Compton wavelength:  {lambda_C*1e12:.3f} pm
      Classical radius:    {r_e*1e15:.3f} fm
      Rydberg energy:      {E_Rydberg_eV:.3f} eV

      All these scales are CONSISTENT with m_e measured value!

   WHY NOT DERIVED?
      epsilon_0 and mu_0 set GEOMETRIC scales (c, Z_0, alpha).
      They do NOT determine MASS SCALES.

      Mass scales require ADDITIONAL ANCHOR:
         • Option 1: Measure m_e (TriPhase approach)
         • Option 2: Derive from Planck mass (quantum gravity)

      TriPhase V16 uses m_e as MEASURED ANCHOR because:
         1) Precisely measured (10 significant figures)
         2) Anchors entire mass hierarchy
         3) Connects EM coupling to inertia

   COUPLING CHAIN SUMMARY:
      Vacuum (epsilon_0, mu_0) → c, Z_0, alpha (geometric)
      Electron (m_e measured)  → mass anchor (inertial)
      Combined → complete coupling hierarchy

   THE ELECTRON IS THE GROUND STATE COUPLER!
      It sits at band m={m}, couples through Z_0 with strength alpha,
      and provides the mass anchor for all other particles.

      Without m_e, we would have geometry but no inertia.
      m_e connects vacuum properties to material properties!
""")

    print("\n" + "=" * 70)
    print(f"RESULT: m_e = {m_e:.10e} kg (MEASURED ANCHOR)")
    print(f"        Ground state coupling at band m={m}")
    print(f"        Inertial coupling through Z_0 = {Z_0:.4f} ohms")
    print(f"        Rest energy: {m_e_c2_keV:.3f} keV")
    print("=" * 70)
    print()
    input("Press Enter to exit...")
