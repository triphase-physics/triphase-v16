# -*- coding: utf-8 -*-
"""
TriPhase V16 Individual Derivation
============================================================
Constant: Down Quark Mass (Steady State)
Symbol: m_d
Row: 20
Framework: WaveMechanics_Primitive

Derived: 4.675 MeV
Measured: 4.67 +0.48/-0.17 MeV (PDG MS-bar at 2 GeV)
Source: m_u * 17/8 (isospin step from up quark)

Tag: (D*) Via 17 from alpha

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
# UP QUARK MASS (FROM ROW 19)
# ============================================================
m_u = 2.20  # MeV (steady state confined mass)

# ============================================================
# MEASURED CALIBRATION CHECKPOINT
# ============================================================
# PDG 2024: m_d = 4.67 +0.48/-0.17 MeV (MS-bar at 2 GeV)
m_d_PDG = 4.67  # MeV (central value)
m_d_PDG_upper = 4.67 + 0.48
m_d_PDG_lower = 4.67 - 0.17

# ============================================================
# THE FOUR WAVE PRIMITIVES
# ============================================================

# 1. FREQUENCY
#    f_d = f_u * (17/8) (isospin step)

# 2. WAVELENGTH
#    lambda_d = lambda_u * (8/17) (shorter than up)

# 3. AMPLITUDE
#    m_d = m_u * (17/8) = 4.675 MeV

# 4. PHASE
#    Standing wave confined by QCD, isospin partner of up

# ============================================================
# TRIPHASE WAVE MECHANICS DERIVATION
# ============================================================
#
# MECHANISM:
# The down quark is the isospin partner of the up quark.
# Their mass ratio is determined by the TriPhase isospin step:
#
# FORMULA (D* - Via 17 from alpha):
#   m_d = m_u * (17/8)
#
# WHERE:
#   17 = from alpha (8*17+1 = 137)
#   8  = 2^3 (generator cube)
#   17/8 = 2.125 (isospin mass step)
#
# WHY THIS RATIO?
# The number 17 appears in alpha_inv = 137 = 8*17 + 1.
# The isospin doublet (u,d) has mass ratio 17/8, where
# 8 = 2^3 represents the three-dimensional generator structure.
#
# This is NOT arbitrary - it's built into the wave mechanics
# of the vacuum via alpha.
#
# ============================================================

# Convert electron mass to MeV
m_e_MeV = m_e * c**2 / eV / 1e6

# TriPhase derivation
isospin_numerator = 17  # From alpha: 137 = 8*17 + 1
isospin_denominator = 8  # 2^3 = generator cube
isospin_ratio = isospin_numerator / isospin_denominator

m_d_triphase = m_u * isospin_ratio

# Frequency
f_u = m_u * 1e6 * eV / h
f_d = f_u * isospin_ratio

# Wavelength
lambda_u = h / (m_u * 1e6 * eV / c)
lambda_d = lambda_u / isospin_ratio

# Error calculation
error_percent = (m_d_triphase - m_d_PDG) / m_d_PDG * 100

# ============================================================
# OUTPUT
# ============================================================
if __name__ == "__main__":
    print("=" * 70)
    print("TRIPHASE V16 -- ROW 20: DOWN QUARK MASS")
    print("Framework: WaveMechanics_Primitive")
    print("Tag: (D*) Via 17 from alpha")
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
    print(f"\n   NOTE: 137 = 8*17 + 1")
    print(f"         The 17 appears in isospin ratio!")

    print("\n" + "-" * 70)
    print("STEP 3: UP QUARK MASS (FROM ROW 19)")
    print("-" * 70)
    print(f"\n   m_u = {m_u:.2f} MeV (base of quark ladder)")

    print("\n" + "-" * 70)
    print("FOUR WAVE PRIMITIVES")
    print("-" * 70)

    print("\n1. FREQUENCY:")
    print(f"   f_u = {f_u:.6e} Hz")
    print(f"   f_d = f_u * {isospin_ratio:.3f} = {f_d:.6e} Hz")

    print("\n2. WAVELENGTH:")
    print(f"   lambda_u = {lambda_u:.6e} m")
    print(f"   lambda_d = lambda_u / {isospin_ratio:.3f} = {lambda_d:.6e} m")

    print("\n3. AMPLITUDE:")
    print(f"   m_d = {m_d_triphase:.3f} MeV")

    print("\n4. PHASE:")
    print("   Standing wave confined by QCD")
    print("   Isospin partner of up quark")

    print("\n" + "-" * 70)
    print("STEP 4: DERIVE DOWN QUARK MASS")
    print("-" * 70)

    print("\n   Formula: m_d = m_u * (17/8)")
    print(f"\n   WHERE:")
    print(f"   - 17 comes from alpha: 137 = 8*17 + 1")
    print(f"   - 8 = 2^3 (generator cube)")
    print(f"   - 17/8 = {isospin_ratio:.6f} (isospin step)")
    print(f"\n   m_d = {m_u:.2f} MeV * {isospin_ratio:.6f}")
    print(f"       = {m_d_triphase:.3f} MeV")

    print("\n" + "-" * 70)
    print("CALIBRATION CHECKPOINT")
    print("-" * 70)
    print(f"\n   TriPhase m_d: {m_d_triphase:.3f} MeV (steady state)")
    print(f"   PDG m_d:      {m_d_PDG:.2f} MeV (MS-bar at 2 GeV)")
    print(f"   PDG range:    {m_d_PDG_lower:.2f} - {m_d_PDG_upper:.2f} MeV")
    print(f"   Error: {error_percent:+.3f}%")

    print(f"\n   Measured ratio: m_d/m_u = {m_d_PDG}/{m_u:.2f} = {m_d_PDG/m_u:.4f}")
    print(f"   TriPhase ratio: 17/8 = {isospin_ratio:.4f}")
    print(f"   Agreement: {(1 - abs(m_d_PDG/m_u - isospin_ratio)/isospin_ratio)*100:.2f}%")

    print("\n" + "-" * 70)
    print("PHYSICAL MEANING")
    print("-" * 70)
    print("\n   - Down quark = isospin partner of up quark")
    print("   - Mass ratio 17/8 from vacuum alpha structure")
    print("   - 17 appears in: alpha^-1 = 137 = 8*17 + 1")
    print("   - 8 = 2^3 = three-dimensional generator")
    print("\n   ISOSPIN DOUBLET:")
    print(f"   - Up quark:   charge +2/3, mass {m_u:.2f} MeV")
    print(f"   - Down quark: charge -1/3, mass {m_d_triphase:.2f} MeV")
    print(f"   - Mass ratio: {isospin_ratio:.4f} (NOT arbitrary!)")
    print("\n   WHY 17/8?")
    print("   - Alpha contains 17: (8*17 + 1)^-1")
    print("   - Isospin is SU(2) symmetry")
    print("   - 8 = 2^3 connects to generator structure")
    print("   - Result: mass splitting determined by alpha")

    print("\n   FOUND IN:")
    print("   - Proton: uud (2 up, 1 down)")
    print("   - Neutron: udd (1 up, 2 down)")
    print("   - Pions: pi+ (ud-bar), pi- (du-bar), pi0 (mix)")

    print("\n   QUARK MASS LADDER (1st generation):")
    print(f"   - Up:   {m_u:.2f} MeV (base)")
    print(f"   - Down: {m_d_triphase:.2f} MeV = up * 17/8")
    print("   Next: strange = down * 20 * 57/56")

    print("\n" + "=" * 70)
    print("STATUS: VIA 17 FROM ALPHA (D*)")
    print(f"        m_d = {m_d_triphase:.3f} MeV")
    print(f"        m_d/m_u = 17/8 = {isospin_ratio:.4f}")
    print(f"        Error: {error_percent:+.3f}%")
    print("        Isospin ratio from alpha")
    print("=" * 70)
    input("Press Enter to exit...")
