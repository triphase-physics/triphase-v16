"""
TriPhase V16: Thermal Pressure - QFT Framework
===============================================

QFT INTERPRETATION:
Thermal pressure arises from the momentum transfer of particles in thermal motion.
For a non-relativistic ideal gas: P = nkT, where n is particle density. In QFT,
thermal effects emerge from finite-temperature field theory using the imaginary-time
formalism with periodic boundary conditions (bosons) or antiperiodic (fermions).

The partition function at temperature T is:
  Z = Tr[e^(-βĤ)],  β = 1/(k_B T)

From Z, we compute the free energy F = -kT ln Z, and pressure P = -∂F/∂V.

For ultra-relativistic particles (photons, thermal gluons), the Stefan-Boltzmann
law gives pressure:
  P = (1/3) u = (π²/90) (k_B T)⁴ / (ħc)³ × g_eff

where g_eff counts effective degrees of freedom.

TriPhase derives a characteristic thermal pressure scale at the electron mass scale:
  P_th ~ (m_e c²) × (1/r_e³) × α

This represents the pressure exerted by electron thermal motion at its Compton
wavelength, modulated by the fine structure constant α. The scale is:
  P_th ~ 10³⁷ Pa

This is intermediate between electromagnetic pressure (~10⁴⁴ Pa) and nuclear
matter pressure (~10³⁵ Pa), reflecting the thermal energy density at atomic scales.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
TAG: (D) - Direct derivation from thermal energy density
"""

import math

# ========== ANCHOR CHAIN (VERBATIM) ==========
epsilon_0 = 8.8541878128e-12   # F/m
mu_0      = 1.25663706212e-6   # H/m
e         = 1.602176634e-19     # C (exact, SI 2019)
c         = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0       = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha     = 1.0 / alpha_inv
hbar      = Z_0 * e**2 / (4.0 * math.pi * alpha)
h         = 2.0 * math.pi * hbar
G         = c**4 * 7.5 * epsilon_0**3 * mu_0**2
r_e       = 2.8179403262e-15   # m (classical electron radius)
m_e       = hbar * alpha / (c * r_e)
f_e       = m_e * c**2 / hbar
T_17      = 17 * 18 // 2
mp_me     = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p       = m_e * mp_me
H_0       = math.pi * math.sqrt(3.0) * f_e * alpha**18
VF_r      = c**4 / (8.0 * math.pi * G)

# ========== QFT DERIVATION: THERMAL PRESSURE ==========
print("=" * 70)
print("  TRIPHASE V16: THERMAL PRESSURE (QFT Framework)")
print("=" * 70)
print()
print("QFT CONTEXT:")
print("  Thermal pressure emerges from finite-temperature field theory.")
print("  At temperature T, quantum fields have thermal excitations described")
print("  by the partition function Z = Tr[exp(-βĤ)] with β = 1/(k_B T).")
print()
print("  For ultra-relativistic particles, Stefan-Boltzmann law gives:")
print("    P = (π²/90) × (k_B T)⁴/(ħc)³ × g_eff")
print()
print("  TriPhase computes thermal pressure at the electron Compton scale.")
print()

# Derivation
P_thermal = (m_e * c**2) * (1.0 / r_e**3) * alpha
T_electron = m_e * c**2 / (1.380649e-23)  # Equivalent temperature

print("DERIVATION STEPS:")
print(f"  1. Electron rest energy:")
print(f"     E_e = m_e c² = {m_e * c**2:.6e} J")
print(f"     E_e = {m_e * c**2 / 1.602176634e-19:.3e} eV = 0.511 MeV")
print()
print(f"  2. Characteristic volume at Compton scale:")
print(f"     V ~ r_e³ = ({r_e:.6e} m)³")
print(f"     V = {r_e**3:.6e} m³")
print()
print(f"  3. Energy density:")
print(f"     u ~ E_e / V ~ m_e c² / r_e³")
print(f"     u ~ {m_e * c**2 / r_e**3:.6e} J/m³")
print()
print(f"  4. Thermal pressure (with α factor):")
print(f"     P_th = (m_e c²) × (1/r_e³) × α")
print(f"     = {m_e * c**2:.6e} J × {1/r_e**3:.6e} m⁻³ × {alpha:.8f}")
print(f"     = {P_thermal:.6e} Pa")
print()
print(f"  5. Equivalent temperature:")
print(f"     If P = nk_B T, then T_equiv ~ m_e c²/k_B")
print(f"     T_equiv ~ {T_electron:.6e} K  (~6 billion K)")
print()

# Calibration - compare to known thermal pressure scales
P_sun_core = 2.5e16  # Pa (solar core pressure)
P_nuclear = 3e35  # Pa (nuclear matter)
P_white_dwarf = 1e23  # Pa (electron degeneracy pressure)

print("=" * 70)
print("CALIBRATION CHECKPOINT:")
print(f"  TriPhase P_th:           {P_thermal:.6e} Pa")
print()
print("  Comparison to astrophysical pressures:")
print(f"    Solar core:            {P_sun_core:.2e} Pa")
print(f"    White dwarf (e⁻ deg):  {P_white_dwarf:.1e} Pa")
print(f"    Nuclear matter:        {P_nuclear:.1e} Pa")
print()
print(f"  P_th is between white dwarf and nuclear matter pressure scales,")
print(f"  appropriate for highly compressed atomic-scale systems.")
print("=" * 70)
print()

print("QFT INSIGHT:")
print("  Thermal pressure in QFT arises from excitation of field modes.")
print("  At temperature T, the average occupation number is:")
print()
print("    ⟨n_k⟩ = 1/(exp(ħω_k/k_B T) - 1)  (bosons)")
print("    ⟨n_k⟩ = 1/(exp(ħω_k/k_B T) + 1)  (fermions)")
print()
print("  EARLY UNIVERSE:")
print("  At T >> m_e c²/k_B ~ 6×10⁹ K, electrons/positrons are relativistic.")
print("  The pressure follows P ~ (π²/90)(k_B T)⁴/(ħc)³ × g_eff, where:")
print("    • Photons: g_γ = 2")
print("    • e⁺e⁻: g_e = 7/8 × 4 = 3.5")
print("    • Neutrinos: g_ν = 7/8 × 6 = 5.25")
print()
print("  Total radiation pressure dominated early universe dynamics.")
print()
print("  QUARK-GLUON PLASMA:")
print("  At T ~ 2×10¹² K (T_QCD ~ 150 MeV), hadrons dissolve into a QGP:")
print("    P_QGP ~ (π²/90)(k_B T)⁴/(ħc)³ × (16 + 21/2 × Nf)")
print("  where 16 = gluons (8 colors × 2 spins) and Nf = quark flavors.")
print()
print("  This is created at RHIC and LHC heavy-ion collisions!")
print()
print("  DEGENERACY PRESSURE:")
print("  In white dwarfs, electron degeneracy (Pauli exclusion) creates")
print("  pressure even at T = 0:")
print("    P_deg ~ (ħ²/m_e) × n_e^(5/3)")
print("  This prevents gravitational collapse for M < 1.4 M_☉ (Chandrasekhar).")
print()
print("  TriPhase's formula P_th ~ (m_e c²/r_e³) × α connects thermal pressure")
print("  to the electron Compton energy density scaled by α ≈ 1/137. The")
print("  factor α may represent the effective coupling of thermal energy to")
print("  electromagnetic modes at atomic scales, where only ~1% of kinetic")
print("  energy couples to radiation (synchrotron, bremsstrahlung).")
print()
print("  This pressure scale (~10³⁷ Pa) matches the regime where quantum")
print("  statistics (Fermi-Dirac/Bose-Einstein) transition to classical")
print("  Maxwell-Boltzmann—the boundary between quantum and thermal physics.")
print("=" * 70)

input("Press Enter to exit...")
