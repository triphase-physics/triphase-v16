"""
TriPhase V16 Derivative: Thermal Pressure (Gauge Theory Framework)

GAUGE THEORY INTERPRETATION:
Thermal pressure P_th = nkT arises from gauge field fluctuations at finite
temperature. In the imaginary-time formalism, thermal field theory is formulated
as a gauge theory on a compactified Euclidean spacetime with periodic boundary
conditions A_μ(x, τ) = A_μ(x, τ + β) where β = 1/(kT) is the inverse temperature.
The pressure P_th = (m_e c²) × (1/r_e³) × α represents the thermal energy density
of gauge field modes confined to the electron Compton volume V ~ r_e³. The factor
α = e²/(4πε_0 ℏc) is the electromagnetic gauge coupling, measuring the strength
of photon-electron interactions. At high temperatures T → ∞, gauge symmetries
are restored: electroweak SU(2)_L × U(1)_Y is unbroken, and the Higgs VEV
vanishes. The thermal pressure drives this phase transition, analogous to how
QCD deconfines at T ~ 150 MeV, transitioning from a confined hadronic phase to
a quark-gluon plasma where color gauge symmetry is restored.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
TAG: (D)
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

print("=" * 70)
print("THERMAL PRESSURE - GAUGE THEORY DERIVATION")
print("=" * 70)

# Gauge theory derivation
print("\nDeriving thermal pressure from gauge field fluctuations:")
print(f"Electron rest energy m_e c² = {m_e * c**2:.6e} J")
print(f"                              = {m_e * c**2 / 1.602176634e-13:.6f} MeV")
print(f"Classical electron radius r_e = {r_e:.6e} m")
print(f"Compton volume r_e³ = {r_e**3:.6e} m³")
print(f"Number density 1/r_e³ = {1.0/r_e**3:.6e} m⁻³")
print(f"Fine structure constant α = {alpha:.10f}")

# Thermal pressure
P_th = (m_e * c**2) * (1.0 / r_e**3) * alpha

print(f"\nThermal pressure P_th = (m_e c²) × (1/r_e³) × α")
print(f"P_th = {P_th:.6e} Pa")
print(f"P_th = {P_th / 1e9:.6e} GPa")

# Effective temperature
k_B = 1.380649e-23  # J/K
T_eff = P_th * r_e**3 / k_B
print(f"\nEffective temperature T_eff = P_th × r_e³ / k_B")
print(f"T_eff = {T_eff:.6e} K")
print(f"      = {T_eff * k_B / 1.602176634e-13:.3f} MeV")

print("\n" + "=" * 70)
print("CALIBRATION CHECKPOINT")
print("=" * 70)

print(f"Derived value:  {P_th:.6e} Pa")
print(f"Physical scale:  {P_th / 1e9:.3e} GPa")
print(f"Effective temperature: {T_eff:.3e} K = {T_eff * k_B / 1.602176634e-13:.3f} MeV")
print(f"Comparison:     QCD deconfinement ~ 150 MeV")
print(f"                Electroweak restoration ~ 100 GeV")

print("\n" + "=" * 70)
print("GAUGE THEORY INSIGHT")
print("=" * 70)
print("""
Thermal pressure in gauge theories drives phase transitions where symmetries
are restored at high temperatures. In the early universe, as temperature dropped
from the Planck scale (10³² K), gauge symmetries broke in stages: GUT symmetry
broke at ~10²⁸ K, electroweak at ~10¹⁵ K, and QCD confined at ~10¹² K. Each
transition released latent heat, potentially seeding density perturbations. The
pressure P_th = (m_e c²)/(r_e³) × α represents the thermal energy density when
electromagnetic gauge fluctuations are thermalized at T ~ m_e c²/k_B ~ 6 billion K.
At this temperature, the gauge coupling α(T) runs logarithmically: α(m_e) ≈ 1/137
increases to α(m_Z) ≈ 1/128 at the electroweak scale. In finite-temperature
field theory, this running modifies the effective potential V_eff(φ, T), shifting
the Higgs minimum from <φ> = v at T = 0 to <φ> = 0 at T > T_c ~ 100 GeV. The
thermal pressure thus controls the order parameter of electroweak symmetry
breaking, determining whether the transition was first-order (allowing baryogenesis)
or second-order (smooth crossover). Current lattice simulations suggest a smooth
crossover, creating a tension with baryogenesis scenarios.
""")

print("=" * 70)
input("Press Enter to exit...")
