# -*- coding: utf-8 -*-
"""
TriPhase V16 Individual Derivation
============================================================
Constant: Strange Quark Mass (Steady State)
Symbol: m_s
Row: 21
Framework: WaveMechanics_Primitive

Derived: 95.1 MeV
Measured: 93.5 +8/-3.4 MeV (PDG MS-bar at 2 GeV)
Source: m_d * 20 * 57/56 (drain rule)

Tag: (D*) Drain rule from down quark

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
============================================================
"""

import numpy as np
import sys
import io
if sys.stdout.encoding != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')


# ============================================================
# FUNDAMENTAL VACUUM CONSTANTS (ONLY INPUTS)
# ============================================================
epsilon_0 = 8.8541878128e-12  # F/m - Electric permittivity
mu_0 = 1.25663706212e-6       # H/m - Magnetic permeability

# ============================================================
# SI-DEFINED EXACT CONSTANTS
# ============================================================
h = 6.62607015e-34   # J*s (exact by SI definition 2019)
hbar = h / (2.0 * np.pi)
e = 1.602176634e-19  # C (exact by SI definition 2019)
eV = 1.602176634e-19 # J (exact by SI definition 2019)

# ============================================================
# ELECTRON MASS ANCHOR
# ============================================================
m_e = 9.1093837015e-31  # kg (measured anchor)

# ============================================================
# DERIVE c AND Z_0
# ============================================================
c = 1.0 / np.sqrt(epsilon_0 * mu_0)  # Speed of light
Z_0 = np.sqrt(mu_0 / epsilon_0)      # Vacuum impedance

# ============================================================
# DERIVE ALPHA (FINE STRUCTURE CONSTANT)
# ============================================================
# TriPhase formula: m=17, node=8*17+1=137
m = 17
node = 8 * m + 1  # = 137
correction = np.log(node) / node
alpha_inv = node + correction
alpha = 1.0 / alpha_inv

# ============================================================
# QUARK MASSES (FROM PREVIOUS ROWS)
# ============================================================
m_u = 2.20   # MeV (Row 19)
m_d = 4.675  # MeV (Row 20)

# ============================================================
# MEASURED CALIBRATION CHECKPOINT
# ============================================================
# PDG 2024: m_s = 93.5 +8/-3.4 MeV (MS-bar at 2 GeV)
m_s_PDG = 93.5  # MeV (central value)
m_s_PDG_upper = 93.5 + 8.0
m_s_PDG_lower = 93.5 - 3.4

# ============================================================
# THE FOUR WAVE PRIMITIVES
# ============================================================

# 1. FREQUENCY
#    f_s = f_d * 20 * 57/56 (drain rule step)

# 2. WAVELENGTH
#    lambda_s = lambda_d * 56/(20*57) (much shorter)

# 3. AMPLITUDE
#    m_s = m_d * 20 * 57/56 = 95.1 MeV

# 4. PHASE
#    Standing wave confined by QCD, second generation

# ============================================================
# TRIPHASE WAVE MECHANICS DERIVATION
# ============================================================
#
# MECHANISM:
# The strange quark is the second-generation partner of the down
# quark. The mass step uses the TriPhase "drain rule":
#
# FORMULA (D* - Drain Rule):
#   m_s = m_d * 20 * (57/56)
#
# WHERE:
#   20    = Generation step factor
#   57/56 = Drain rule correction
#
# WHY THIS FORMULA?
# The drain rule 57/56 appears throughout TriPhase physics:
# - Saturn's hexagon at 57° latitude
# - CMB acoustic peaks
# - Pressure band stability
#
# The factor 20 represents the generation step from first to
# second generation for down-type quarks.
#
# ============================================================

# Convert electron mass to MeV
m_e_MeV = m_e * c**2 / eV / 1e6

# TriPhase derivation
generation_step = 20
drain_numerator = 57
drain_denominator = 56
drain_ratio = drain_numerator / drain_denominator

m_s_triphase = m_d * generation_step * drain_ratio

# Frequency
f_d = m_d * 1e6 * eV / h
f_s = f_d * generation_step * drain_ratio

# Wavelength
lambda_d = h / (m_d * 1e6 * eV / c)
lambda_s = lambda_d / (generation_step * drain_ratio)

# Error calculation
error_percent = (m_s_triphase - m_s_PDG) / m_s_PDG * 100

# ============================================================
# OUTPUT
# ============================================================
if __name__ == "__main__":
    print("=" * 70)
    print("TRIPHASE V16 -- ROW 21: STRANGE QUARK MASS")
    print("Framework: WaveMechanics_Primitive")
    print("Tag: (D*) Drain rule from down quark")
    print("=" * 70)

    print("\n" + "-" * 70)
    print("STEP 1: DERIVE VACUUM PROPERTIES")
    print("-" * 70)
    print(f"\n   INPUT: epsilon_0 = {epsilon_0:.13e} F/m")
    print(f"   INPUT: mu_0      = {mu_0:.11e} H/m")
    print(f"\n   c   = 1/sqrt(epsilon_0 * mu_0) = {c:.10e} m/s")
    print(f"   Z_0 = sqrt(mu_0/epsilon_0)      = {Z_0:.10f} Ohm")

    print("\n" + "-" * 70)
    print("STEP 2: DERIVE ALPHA (FINE STRUCTURE CONSTANT)")
    print("-" * 70)
    print(f"\n   m = {m}")
    print(f"   node = 8*m + 1 = 8*{m} + 1 = {node}")
    print(f"   correction = ln({node})/{node} = {correction:.10f}")
    print(f"   alpha_inv = {node} + {correction:.10f} = {alpha_inv:.10f}")
    print(f"   alpha = 1/alpha_inv = {alpha:.15f}")

    print("\n" + "-" * 70)
    print("STEP 3: QUARK MASSES (FROM PREVIOUS ROWS)")
    print("-" * 70)
    print(f"\n   m_u = {m_u:.2f} MeV (Row 19 - base)")
    print(f"   m_d = {m_d:.3f} MeV (Row 20 - isospin partner)")

    print("\n" + "-" * 70)
    print("FOUR WAVE PRIMITIVES")
    print("-" * 70)

    print("\n1. FREQUENCY:")
    print(f"   f_d = {f_d:.6e} Hz")
    print(f"   f_s = f_d * 20 * {drain_ratio:.6f} = {f_s:.6e} Hz")

    print("\n2. WAVELENGTH:")
    print(f"   lambda_d = {lambda_d:.6e} m")
    print(f"   lambda_s = lambda_d / (20 * {drain_ratio:.6f}) = {lambda_s:.6e} m")

    print("\n3. AMPLITUDE:")
    print(f"   m_s = {m_s_triphase:.1f} MeV")

    print("\n4. PHASE:")
    print("   Standing wave confined by QCD")
    print("   Second generation, down-type quark")

    print("\n" + "-" * 70)
    print("STEP 4: DERIVE STRANGE QUARK MASS")
    print("-" * 70)

    print("\n   Formula: m_s = m_d * 20 * (57/56)")
    print(f"\n   WHERE:")
    print(f"   - m_d = {m_d:.3f} MeV (down quark mass)")
    print(f"   - 20 = generation step factor")
    print(f"   - 57/56 = {drain_ratio:.10f} (drain rule)")
    print(f"\n   m_s = {m_d:.3f} MeV * 20 * {drain_ratio:.10f}")
    print(f"       = {m_s_triphase:.3f} MeV")

    print("\n   THE DRAIN RULE (57/56):")
    print("   - Appears in Saturn's hexagon (57° latitude)")
    print("   - CMB acoustic peak ratios")
    print("   - Pressure band stability")
    print("   - Universal TriPhase correction factor")

    print("\n" + "-" * 70)
    print("CALIBRATION CHECKPOINT")
    print("-" * 70)
    print(f"\n   TriPhase m_s: {m_s_triphase:.1f} MeV (steady state)")
    print(f"   PDG m_s:      {m_s_PDG:.1f} MeV (MS-bar at 2 GeV)")
    print(f"   PDG range:    {m_s_PDG_lower:.1f} - {m_s_PDG_upper:.1f} MeV")
    print(f"   Error: {error_percent:+.2f}%")

    print(f"\n   Measured ratio: m_s/m_d = {m_s_PDG}/{m_d:.2f} = {m_s_PDG/m_d:.2f}")
    print(f"   TriPhase ratio: 20 * 57/56 = {generation_step * drain_ratio:.2f}")

    print("\n" + "-" * 70)
    print("PHYSICAL MEANING")
    print("-" * 70)
    print("\n   - Strange quark = 2nd generation of down quark")
    print("   - Mass step: 20 * (57/56) from down quark")
    print("   - Drain rule 57/56 is universal TriPhase factor")
    print("\n   QUARK GENERATIONS:")
    print("   1st generation (light):")
    print(f"     Up-type:   u = {m_u:.2f} MeV")
    print(f"     Down-type: d = {m_d:.2f} MeV")
    print("   2nd generation (strange):")
    print(f"     Up-type:   c = 1267 MeV (charm)")
    print(f"     Down-type: s = {m_s_triphase:.1f} MeV (strange)")
    print("   3rd generation (heavy):")
    print("     Up-type:   t = 172.5 GeV (top)")
    print("     Down-type: b = 4.18 GeV (bottom)")

    print("\n   WHY 'STRANGE'?")
    print("   - Discovered in 1947 in cosmic rays")
    print("   - Called 'strange' because it lived longer than expected")
    print("   - Strangeness conservation in strong/EM interactions")
    print("   - Only weak force can change strangeness")

    print("\n   FOUND IN:")
    print("   - Kaons: K+ (us-bar), K- (su-bar)")
    print("   - Lambda: uds (strangeness -1)")
    print("   - Sigma baryons: various sss, ssu, ssd")
    print("   - Omega-: sss (strangeness -3)")

    print("\n   QUARK MASS LADDER (down-type):")
    print(f"   - Down:    {m_d:.2f} MeV (1st gen)")
    print(f"   - Strange: {m_s_triphase:.1f} MeV = down * 20 * 57/56")
    print("   - Bottom:  4180 MeV = down * 880 * 57/56")

    print("\n" + "=" * 70)
    print("STATUS: DRAIN RULE (D*)")
    print(f"        m_s = {m_s_triphase:.1f} MeV")
    print(f"        m_s/m_d = 20 * 57/56 = {generation_step * drain_ratio:.3f}")
    print(f"        Error: {error_percent:+.2f}%")
    print("        Drain rule 57/56")
    print("=" * 70)
    input("Press Enter to exit...")
