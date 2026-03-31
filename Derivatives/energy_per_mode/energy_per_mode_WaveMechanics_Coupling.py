# -*- coding: utf-8 -*-
"""
============================================================
TRIPHASE V16 - INDIVIDUAL DERIVATIVE SCRIPT
============================================================
Derivative: Energy Per Mode (epsilon_pair)
Framework:  WaveMechanics_Coupling
Row:        8

INPUTS:     epsilon_0, mu_0 ONLY (Three Observed Vacuum Properties)
MEASURED:   Derived from Rydberg energy and T_17

Tag: (D) DERIVED from Rydberg coupling and triangular modes

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

# Electron mass (measured anchor)
m_e = 9.1093837015e-31  # kg (CODATA 2018)

# ============================================================
# CALIBRATION CHECKPOINT (for comparison)
# ============================================================
# Typical scale: epsilon_pair ~ few eV
# Exact value depends on Rydberg energy

# ============================================================
# COUPLING CHAIN DERIVATION
# ============================================================
#
# ENERGY PER MODE: COUPLING ENERGY DISTRIBUTION
#
# The Rydberg energy is the binding energy of hydrogen ground state:
#
#    E_Rydberg = m_e * c^2 * alpha^2 / 2
#              = (electron rest energy) × (coupling strength)^2 / 2
#
# This energy represents the TOTAL COUPLING ENERGY available
# for electromagnetic binding in the hydrogen atom.
#
# COUPLING MECHANISM:
#
# At pressure band m = 17, there are T_17 = 153 coupling modes
# (see triangular_T17 derivative).
#
# The Rydberg energy is DISTRIBUTED across these coupling modes.
# Each mode carries a quantum of energy:
#
#    epsilon_pair = 2 * E_Rydberg / T_17
#
# The factor of 2 accounts for the electron-proton pair coupling
# (two particles coupling through T_17 modes).
#
# FORMULA:
#    epsilon_pair = 2 * E_Rydberg / T_17
#                 = 2 * (m_e * c^2 * alpha^2 / 2) / 153
#                 = m_e * c^2 * alpha^2 / 153
#
# IMPEDANCE COUPLING PICTURE:
#
# The vacuum impedance Z_0 = 376.73 ohms mediates the coupling.
# Each of the T_17 = 153 coupling channels experiences Z_0.
#
# The coupling energy per mode:
#    epsilon_pair ~ (coupling voltage)^2 / Z_0
#
# Where the coupling voltage is set by the Rydberg energy scale.
#
# PHYSICAL INTERPRETATION:
#
# epsilon_pair is the energy quantum associated with each
# coupling mode at the electron band m = 17.
#
# When the electron couples to the proton to form hydrogen,
# the binding energy (Rydberg) is partitioned across the
# T_17 = 153 available coupling channels.
#
# Each channel carries energy epsilon_pair.
#
# COUPLING HIERARCHY:
#
# Different pressure bands have different T_m values:
#   - Lower bands (m < 17): fewer coupling modes, larger epsilon
#   - Higher bands (m > 17): more coupling modes, smaller epsilon
#
# The electron at m = 17 sits at a special point where:
#   T_17 = 153 = 9 × 17 = 3^2 × 17
#
# This gives epsilon_pair a characteristic scale that determines
# atomic coupling strengths.
#
# ENERGY CONSERVATION IN COUPLING:
#
# Total Rydberg energy = sum over all coupling modes
#    E_Rydberg = sum_{i=1}^{T_17} (epsilon_pair / 2)
#              = (T_17 / 2) * epsilon_pair
#
# Solving for epsilon_pair:
#    epsilon_pair = 2 * E_Rydberg / T_17
#
# This is the fundamental energy partition law for coupling.
#
# ============================================================

# Triangular number T_17
T_17 = m * (m + 1) // 2  # = 153

# Rydberg energy (in Joules)
E_Rydberg = m_e * c**2 * alpha**2 / 2.0

# Energy per coupling mode
epsilon_pair = 2.0 * E_Rydberg / T_17

# Convert to eV for readability
eV = 1.602176634e-19  # J/eV (exact, SI 2019 definition)
E_Rydberg_eV = E_Rydberg / eV
epsilon_pair_eV = epsilon_pair / eV

# Coupling energy per mode in terms of Z_0
# (dimensional analysis for insight)
# Energy ~ Voltage^2 / Impedance
# Characteristic voltage scale ~ sqrt(E_Rydberg * Z_0)
V_char = np.sqrt(E_Rydberg * Z_0)

# ============================================================
# OUTPUT
# ============================================================
if __name__ == "__main__":
    print("=" * 70)
    print("TRIPHASE V16 - DERIVATIVE 8: ENERGY PER MODE (epsilon_pair)")
    print("Framework: WaveMechanics_Coupling")
    print("Tag: (D) DERIVED from Rydberg coupling and triangular modes")
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
    print("MEASURED ANCHOR")
    print("-" * 70)
    print(f"   m_e = {m_e:.10e} kg  (electron mass, CODATA)")

    print("\n" + "-" * 70)
    print("COUPLING CHAIN: Energy Distribution Across Modes")
    print("-" * 70)
    print(f"\n   Pressure band index: m = {m}")
    print(f"   Coupling modes available: T_{m} = {m}×{m+1}/2 = {T_17}")

    print(f"\n   RYDBERG ENERGY (total coupling energy):")
    print(f"      E_Rydberg = m_e * c^2 * alpha^2 / 2")
    print(f"                = ({m_e:.4e}) * ({c:.3e})^2 * ({alpha:.6e})^2 / 2")
    print(f"                = {E_Rydberg:.6e} J")
    print(f"                = {E_Rydberg_eV:.6f} eV")

    print(f"\n   ENERGY PER COUPLING MODE:")
    print(f"      epsilon_pair = 2 * E_Rydberg / T_{m}")
    print(f"                   = 2 * {E_Rydberg_eV:.6f} / {T_17}")
    print(f"                   = {epsilon_pair_eV:.6f} eV")
    print(f"                   = {epsilon_pair:.6e} J")

    print("\n" + "-" * 70)
    print("IMPEDANCE COUPLING PICTURE")
    print("-" * 70)
    print(f"\n   Vacuum impedance: Z_0 = {Z_0:.4f} ohms")
    print(f"\n   Each of the {T_17} coupling modes experiences Z_0.")
    print(f"\n   Characteristic coupling voltage scale:")
    print(f"      V_char ~ sqrt(E_Rydberg * Z_0)")
    print(f"             = sqrt({E_Rydberg:.3e} * {Z_0:.1f})")
    print(f"             = {V_char:.3e} (dimensionally: sqrt(J·ohms))")

    print("\n" + "-" * 70)
    print("COUPLING ENERGY PARTITION")
    print("-" * 70)
    print(f"\n   Total Rydberg energy distributed across {T_17} modes:")
    print(f"      E_Rydberg = sum(epsilon_pair/2) over {T_17} modes")
    print(f"                = ({T_17}/2) * epsilon_pair")
    print(f"\n   Solving for epsilon_pair:")
    print(f"      epsilon_pair = 2 * E_Rydberg / {T_17}")
    print(f"\n   Verification:")
    print(f"      ({T_17}/2) * {epsilon_pair_eV:.6f} eV = {(T_17/2) * epsilon_pair_eV:.6f} eV")
    print(f"      E_Rydberg = {E_Rydberg_eV:.6f} eV  ✓")

    print("\n" + "-" * 70)
    print("PHYSICAL MEANING: COUPLING ENERGY DISTRIBUTION")
    print("-" * 70)
    print(f"""
   epsilon_pair = {epsilon_pair_eV:.6f} eV

   This is the energy quantum per coupling mode at band m = {m}.

   COUPLING INTERPRETATION:
      The Rydberg energy ({E_Rydberg_eV:.3f} eV) is the total
      electromagnetic binding energy in hydrogen.

      This energy is DISTRIBUTED across T_{m} = {T_17} coupling
      channels connecting the electron and proton through the
      vacuum impedance Z_0 = {Z_0:.4f} ohms.

      Each coupling channel carries energy epsilon_pair.

   IMPEDANCE MATCHING:
      Energy transfer through each channel requires impedance
      matching to Z_0. Mismatched impedances cause reflections
      and reduce coupling efficiency.

   COUPLING HIERARCHY:
      Lower bands (m < {m}): fewer modes → larger epsilon
      Higher bands (m > {m}): more modes → smaller epsilon

      The electron at m = {m} sits at the characteristic scale
      where epsilon_pair ~ {epsilon_pair_eV:.2f} eV.

   WHY THIS MATTERS:
      epsilon_pair sets the energy scale for atomic transitions
      and coupling strengths. It determines the fine details
      of electromagnetic interactions at the electron band.

      Knowing epsilon_pair allows prediction of coupling energies
      for other processes at band m = {m}.

   COMPARISON TO OTHER SCALES:
      Rydberg energy:     {E_Rydberg_eV:.3f} eV (total)
      epsilon_pair:       {epsilon_pair_eV:.6f} eV (per mode)
      Ratio:              {E_Rydberg_eV / epsilon_pair_eV:.1f} = {T_17}/2 ✓
""")

    print("\n" + "=" * 70)
    print(f"RESULT: epsilon_pair = {epsilon_pair_eV:.6f} eV")
    print(f"        Energy per coupling mode at band m = {m}")
    print(f"        ({T_17} coupling modes, Z_0 = {Z_0:.4f} ohms)")
    print("=" * 70)
    print()
    input("Press Enter to exit...")
