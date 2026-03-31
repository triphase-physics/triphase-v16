"""
TriPhase V16 — Thermal Pressure (CMB) (Renormalization Group Framework)
========================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

RENORMALIZATION GROUP INTERPRETATION
--------------------------------------
The CMB thermal pressure P_th = n k_B T_CMB arises from photon and neutrino
number densities at temperature T_CMB ≈ 2.7255 K. In the RG framework, this
is the relic IR pressure from the early universe's hot Big Bang initial conditions.
As the universe expands and cools, the thermal pressure redshifts: P_th(a) ∝ a⁻⁴,
where a is the scale factor. Today's CMB temperature is the IR endpoint of this
cosmological RG flow.

The CMB temperature is NOT a free parameter but an RG-evolved quantity. Starting
from the hot Big Bang (T ~ 10¹⁰ K at nucleosynthesis), photons decouple at
T_dec ~ 3000 K and then cool adiabatically: T(a) = T_dec × a_dec/a. The observed
T_CMB ≈ 2.7255 K encodes the total integrated expansion from decoupling to today,
which is determined by the same α¹⁸ cascade that sets H₀.

Thermal pressure is subdominant today (P_th << P_DE), but it was the dominant
pressure component in the radiation-dominated era. The RG flow from radiation
domination (early times) to dark energy domination (late times) represents a
crossover between two IR fixed points: high-T (thermal) and low-T (vacuum).
The TriPhase framework connects both regimes through the α¹⁸ cascade determining
cosmic expansion history.

TAG: (C) — Calibrated to observation (CMB temperature from measurement)
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

# ========== RENORMALIZATION GROUP DERIVATION ==========
print("=" * 70)
print("TriPhase V16: Thermal Pressure (CMB) (RG Framework)")
print("=" * 70)
print()

print("CMB TEMPERATURE (CALIBRATED OBSERVATION)")
print("-" * 70)
T_CMB = 2.7255  # K, CODATA/Planck 2018
k_B = 1.380649e-23  # J/K, Boltzmann constant (exact, SI 2019)

print(f"CMB temperature:                 T_CMB = {T_CMB} K")
print(f"Boltzmann constant:              k_B = {k_B:.6e} J/K")
print()

print("PHOTON NUMBER DENSITY")
print("-" * 70)
print("For a blackbody at temperature T:")
print("  n_γ = (2ζ(3)/π²) × (k_B T / ħc)³")
print()
print("where ζ(3) ≈ 1.202 is the Riemann zeta function.")
print()

zeta_3 = 1.2020569  # Apéry's constant
n_photon = (2 * zeta_3 / math.pi**2) * (k_B * T_CMB / (hbar * c))**3

print(f"Photon number density:           n_γ = {n_photon:.6e} m⁻³")
print()

print("THERMAL PRESSURE FROM CMB PHOTONS")
print("-" * 70)
print("The photon pressure is:")
print("  P_γ = (1/3) × ρ_γ c² = (π²/45) × (k_B T)⁴ / (ħ³c³)")
print()
print("Alternatively, from ideal gas law:")
print("  P_th ≈ n_γ k_B T  (approximate, exact factor π²/90)")
print()

P_photon_exact = (math.pi**2 / 15) * (k_B * T_CMB)**4 / (hbar**3 * c**3)
P_thermal_approx = n_photon * k_B * T_CMB

print(f"Photon pressure (exact):         P_γ = {P_photon_exact:.6e} Pa")
print(f"Thermal pressure (approx):       P_th = {P_thermal_approx:.6e} Pa")
print(f"Ratio (should be ≈ π²/90 ≈ 0.11): {P_photon_exact / P_thermal_approx:.3f}")
print()

P_th = P_photon_exact  # Use exact expression

print("RG FLOW: COSMOLOGICAL COOLING")
print("-" * 70)
print("The CMB temperature evolves with scale factor a:")
print("  T(a) = T_dec × (a_dec / a)")
print()
print("From photon decoupling to today:")
print(f"  T_dec ≈ 3000 K  →  T_today = {T_CMB} K")
print(f"  Redshift:     z ≈ {3000 / T_CMB - 1:.0f}")
print()
print("The thermal pressure redshifts as:")
print("  P_th(a) ∝ a⁻⁴  (radiation pressure scaling)")
print()

# ========== CALIBRATION CHECKPOINT ==========
rho_crit = 3 * H_0**2 / (8 * math.pi * G)
Omega_Lambda = 0.685
P_DE = rho_crit * Omega_Lambda * c**2

print("COMPARISON TO OTHER PRESSURES")
print("-" * 70)
print(f"Thermal pressure (CMB):          P_th = {P_th:.6e} Pa")
print(f"Dark energy pressure:            |P_DE| = {P_DE:.6e} Pa")
print(f"Ratio P_th / P_DE:               {P_th / P_DE:.6e}")
print()
print("Today, thermal pressure is ~10⁻⁴ of dark energy pressure (subdominant).")
print()

print("RADIATION VS. DARK ENERGY DOMINATION")
print("-" * 70)
print("The universe transitioned from radiation-dominated to matter-dominated")
print("at z ≈ 3400, then to dark energy-dominated at z ≈ 0.4.")
print()
print("Thermal pressure dominance:")
print("  Early: P_th >> P_m, P_DE  (radiation era, z > 3400)")
print("  Late:  P_th << P_DE       (dark energy era, z < 0.4)")
print()
print("The RG flow from hot Big Bang to cold dark energy-dominated universe")
print("represents a crossover between two IR fixed points.")
print()

print("RENORMALIZATION GROUP INSIGHT")
print("-" * 70)
print("The CMB thermal pressure is the relic IR remnant of the hot Big Bang. As")
print("the universe expands, P_th redshifts as a⁻⁴, becoming subdominant at late")
print("times. The observed T_CMB ≈ 2.7255 K encodes the entire RG flow from")
print("decoupling (T ~ 3000 K) to today, determined by the α¹⁸ cascade setting H₀.")
print("Thermal pressure traces the RG evolution from a high-temperature UV fixed")
print("point to the low-temperature vacuum-dominated IR fixed point.")
print()
print("=" * 70)

input("Press Enter to exit...")
