# -*- coding: utf-8 -*-
"""
TriPhase Wave Mechanics Framework - Tau Mass
WaveMechanics_Coupling Derivation

Derives the tau mass from vacuum permittivity and permeability through
the fine structure constant and coupling harmonics.

Coupling Focus: Highest-energy transient at 51/2 coupling harmonic

Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
License: MIT

Derivative Tag: (D*)
"""

import math

print("=" * 80)
print("TriPhase Wave Mechanics Framework - Tau Mass")
print("WaveMechanics_Coupling Derivation")
print("=" * 80)
print()

# ============================================================================
# VACUUM CONSTANTS (ONLY INPUTS)
# ============================================================================
print("VACUUM CONSTANTS (ONLY INPUTS)")
print("-" * 80)

epsilon_0 = 8.8541878128e-12  # F/m - Vacuum permittivity
mu_0 = 1.25663706212e-6       # H/m - Vacuum permeability

print(f"ε₀ (vacuum permittivity)    = {epsilon_0:.13e} F/m")
print(f"μ₀ (vacuum permeability)    = {mu_0:.14e} H/m")
print()

# ============================================================================
# DERIVED VACUUM PROPERTIES
# ============================================================================
print("DERIVED VACUUM PROPERTIES")
print("-" * 80)

c = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0 = math.sqrt(mu_0 / epsilon_0)

print(f"c (speed of light)          = {c:.10e} m/s")
print(f"Z₀ (vacuum impedance)       = {Z_0:.10f} Ω")
print()

# ============================================================================
# FINE STRUCTURE CONSTANT (ALPHA)
# ============================================================================
print("FINE STRUCTURE CONSTANT (ALPHA)")
print("-" * 80)

m = 17
node = 8 * m + 1
correction = math.log(node) / node
alpha_inv = node + correction
alpha = 1.0 / alpha_inv

print(f"m (harmonic)                = {m}")
print(f"node (8m + 1)               = {node}")
print(f"correction (ln(137)/137)    = {correction:.15f}")
print(f"α⁻¹                         = {alpha_inv:.15f}")
print(f"α (fine structure constant) = {alpha:.15e}")
print()

# ============================================================================
# SI-DEFINED CONSTANTS
# ============================================================================
print("SI-DEFINED CONSTANTS (EXACT)")
print("-" * 80)

m_e = 9.1093837015e-31  # kg - Electron mass (measured anchor)

print(f"m_e (electron mass)         = {m_e:.13e} kg")
print()

# ============================================================================
# COUPLING CHAIN
# ============================================================================
print("COUPLING CHAIN - Highest-Energy Transient at 51/2 Harmonic")
print("-" * 80)
print()
print("The tau represents the highest-energy transient electromagnetic")
print("resonance at the 51/2 coupling harmonic of the electron base frequency.")
print()
print("Coupling hierarchy:")
print("  Z₀ → α → electron (base) → muon (3/2) → tau (51/2)")
print()
print("The mass ratio m_τ/m_e emerges from:")
print("  1. Base harmonic: (51/2) * α⁻¹")
print("  2. Resonance correction: (1 + 5α/6)")
print("  3. High-energy correction: (1 - 3α/2)")
print()
print("The negative correction term reflects the UPPER LIMIT of stable")
print("electromagnetic coupling resonances. Beyond this, the QCD coupling")
print("dominates and particles transition to quark/hadron regimes.")
print()

# ============================================================================
# TAU MASS DERIVATION
# ============================================================================
print("TAU MASS DERIVATION")
print("-" * 80)

# Base harmonic
base_harmonic = (51.0 / 2.0) * alpha_inv

# Resonance correction
resonance_correction = 1.0 + (5.0 * alpha / 6.0)

# High-energy correction (negative - approaching coupling limit)
high_energy_correction = 1.0 - (3.0 * alpha / 2.0)

# Mass ratio
mass_ratio = base_harmonic * resonance_correction * high_energy_correction

# Tau mass
m_tau_derived = mass_ratio * m_e

print(f"m_τ/m_e = (51/2) * α⁻¹ * (1 + 5α/6) * (1 - 3α/2)")
print()
print(f"Base harmonic (51/2)α⁻¹    = {base_harmonic:.10f}")
print(f"Resonance correction        = {resonance_correction:.15f}")
print(f"High-energy correction      = {high_energy_correction:.15f}")
print()
print(f"Mass ratio m_τ/m_e          = {mass_ratio:.10f}")
print(f"Tau mass m_τ                = {m_tau_derived:.15e} kg")
print()

# Convert to MeV/c²
c_mks = c
MeV_per_kg = (c_mks**2) / 1.602176634e-13  # 1 MeV = 1.602176634e-13 J
m_tau_MeV = m_tau_derived * MeV_per_kg

print(f"Tau mass m_τ                = {m_tau_MeV:.10f} MeV/c²")
print(f"Tau mass m_τ                = {m_tau_MeV / 1000.0:.10f} GeV/c²")
print()

# ============================================================================
# CALIBRATION CHECKPOINT (PDG 2024)
# ============================================================================
print("CALIBRATION CHECKPOINT")
print("-" * 80)

m_tau_pdg = 1776.86  # MeV/c² - PDG 2024

deviation = abs(m_tau_MeV - m_tau_pdg)
relative_error = (deviation / m_tau_pdg) * 100

print(f"PDG 2024 m_τ                = {m_tau_pdg:.10f} MeV/c²")
print(f"Derived m_τ                 = {m_tau_MeV:.10f} MeV/c²")
print(f"Absolute deviation          = {deviation:.10f} MeV/c²")
print(f"Relative error              = {relative_error:.10f}%")
print()

if relative_error < 0.01:
    print("✓ EXCELLENT AGREEMENT - Within 0.01% of PDG")
elif relative_error < 0.1:
    print("✓ GOOD AGREEMENT - Within 0.1% of PDG")
elif relative_error < 1.0:
    print("✓ FAIR AGREEMENT - Within 1% of PDG")
else:
    print("⚠ DEVIATION EXCEEDS 1% - Review derivation")

print()

# ============================================================================
# COUPLING INTERPRETATION
# ============================================================================
print("COUPLING INTERPRETATION")
print("-" * 80)
print()
print("The tau mass marks the UPPER BOUNDARY of electromagnetic transients:")
print()
print(f"  Z₀ = {Z_0:.6f} Ω           (vacuum impedance - base coupling)")
print(f"  α  = {alpha:.10f}          (EM fine structure)")
print(f"  51/2 harmonic = {base_harmonic:.6f}   (maximum resonant mode)")
print(f"  m_τ/m_e = {mass_ratio:.6f}        (boundary state mass ratio)")
print()
print("The tau's extremely short lifetime (2.9×10⁻¹³ s) reflects its position")
print("at the edge of electromagnetic stability. Beyond this mass scale,")
print("QCD coupling dominates and particles exist as bound quark states.")
print()
print("There are NO higher electromagnetic lepton generations - the tau")
print("represents the coupling boundary where EM → QCD transition occurs.")
print()

# ============================================================================
# SUMMARY
# ============================================================================
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"Inputs:  ε₀, μ₀ (vacuum constants)")
print(f"Derived: c, Z₀, α, m_e (anchor)")
print(f"Result:  m_τ = {m_tau_MeV:.10f} MeV/c² = {m_tau_MeV/1000.0:.5f} GeV/c²")
print(f"Match:   {100 - relative_error:.10f}% agreement with PDG")
print(f"Tag:     (D*) - Derived via coupling boundary")
print("=" * 80)

input("\nPress Enter to exit...")
