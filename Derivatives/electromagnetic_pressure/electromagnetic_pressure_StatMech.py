"""
TriPhase V16 — Electromagnetic Pressure (Statistical Mechanics Framework)
==========================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
Electromagnetic pressure arises from radiation and electromagnetic field energy
density. In statistical mechanics, a photon gas in thermal equilibrium follows
Bose-Einstein statistics with zero chemical potential (photons are not conserved).
The partition function Z = Tr[exp(-βH)] for a photon gas yields the famous
Stefan-Boltzmann law: energy density u = aT⁴ where a = π²k_B⁴/(15ℏ³c³), and
radiation pressure P = u/3 = (aT⁴)/3. This factor of 1/3 comes from averaging
over three spatial directions for massless particles traveling at speed c.

The electromagnetic field itself stores energy in E and B fields with density
u_EM = (ε_0E²/2 + B²/(2μ_0)). This energy exerts pressure on boundaries, known
as radiation pressure or Maxwell stress. In cosmology, radiation pressure dominated
the early universe (T > 3000 K), with equation of state w = P/ρ = 1/3. The partition
function for blackbody radiation connects microscopic photon occupation numbers
n(ω) = 1/(exp(βℏω) - 1) to macroscopic thermodynamic quantities. TriPhase connects
electromagnetic pressure to fundamental constants through the impedance of free
space Z_0 = √(μ_0/ε_0) and speed of light c = 1/√(ε_0μ_0).

TAG: (D) — Direct TriPhase derivation from pure wave mechanics
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

# ========== STATISTICAL MECHANICS DERIVATION ==========
print("=" * 70)
print("TriPhase V16: Electromagnetic Pressure (Statistical Mechanics)")
print("=" * 70)
print()

print("STATISTICAL MECHANICS FRAMEWORK")
print("--------------------------------")
print("Ensemble: Canonical (photon gas in thermal equilibrium)")
print("Bose-Einstein distribution: n(ω) = 1/(exp(ℏω/k_BT) - 1)")
print("Stefan-Boltzmann law: u = aT⁴, P = u/3")
print("Equation of state: w = P/ρ = 1/3 (ultrarelativistic)")
print()

print("RADIATION PRESSURE FROM PHOTON GAS")
print("-----------------------------------")
k_B = 1.380649e-23  # Boltzmann constant, J/K

# Stefan-Boltzmann constant
a_rad = math.pi**2 * k_B**4 / (15.0 * hbar**3 * c**3)
print(f"Boltzmann constant k_B = {k_B:.6e} J/K")
print(f"Stefan-Boltzmann constant a = π²k_B⁴/(15ℏ³c³)")
print(f"  a = {a_rad:.6e} J/(m³·K⁴)")
print()

# Example: CMB radiation pressure
T_CMB = 2.7255  # K, CMB temperature today
u_CMB = a_rad * T_CMB**4
P_CMB = u_CMB / 3.0

print(f"CMB temperature T_CMB = {T_CMB:.4f} K")
print(f"CMB energy density u = aT⁴ = {u_CMB:.6e} J/m³")
print(f"CMB radiation pressure P = u/3 = {P_CMB:.6e} Pa")
print()

# Radiation density parameter today
rho_crit = 3.0 * H_0**2 / (8.0 * math.pi * G)
rho_rad_CMB = u_CMB / c**2
Omega_rad_CMB = rho_rad_CMB / rho_crit

print(f"Critical density ρ_c = {rho_crit:.6e} kg/m³")
print(f"CMB mass density ρ_rad = u/c² = {rho_rad_CMB:.6e} kg/m³")
print(f"Radiation fraction Ω_rad = {Omega_rad_CMB:.6e}")
print()

# Electromagnetic field energy density and pressure
# For a monochromatic wave: u = ε_0 E_0² (average over cycle)
# Maxwell stress tensor gives pressure P = u/3 for isotropic radiation
print("ELECTROMAGNETIC FIELD PRESSURE")
print("------------------------------")
print("Field energy density: u_EM = ε_0E²/2 + B²/(2μ_0)")
print("For plane wave: E_0 = cB_0, so u = ε_0E_0² (time-averaged)")
print("Radiation pressure: P_EM = u_EM/3 (isotropic averaging)")
print()

# Example: Solar radiation pressure at Earth
# Solar constant S ~ 1361 W/m²
S_solar = 1361.0  # W/m²
P_solar = S_solar / c
print(f"Solar constant at Earth: S = {S_solar:.1f} W/m²")
print(f"Solar radiation pressure: P = S/c = {P_solar:.6e} Pa")
print(f"  (This is what pushes comet tails and powers solar sails)")
print()

# Vacuum impedance and electromagnetic constants
print("VACUUM ELECTROMAGNETIC PROPERTIES")
print("----------------------------------")
print(f"Permittivity ε_0 = {epsilon_0:.6e} F/m")
print(f"Permeability μ_0 = {mu_0:.6e} H/m")
print(f"Vacuum impedance Z_0 = √(μ_0/ε_0) = {Z_0:.6f} Ω")
print(f"Speed of light c = 1/√(ε_0μ_0) = {c:.6e} m/s")
print()

# ========== CALIBRATION CHECKPOINT ==========
print("CALIBRATION CHECKPOINT")
print("----------------------")
# CODATA 2018: a = 7.565723...e-16 J/(m³·K⁴)
a_rad_codata = 7.5657e-16
deviation_a = (a_rad - a_rad_codata) / a_rad_codata * 1e6
print(f"CODATA 2018 Stefan-Boltzmann a: {a_rad_codata:.6e} J/(m³·K⁴)")
print(f"TriPhase calculation:            {a_rad:.6e} J/(m³·K⁴)")
print(f"Deviation:                       {deviation_a:.0f} ppm")
print()

# CMB temperature and density
T_CMB_COBE = 2.7255  # K (COBE/FIRAS measurement)
print(f"COBE/FIRAS CMB temperature: {T_CMB_COBE:.4f} K")
print(f"TriPhase uses this as input for radiation calculations.")
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT")
print("-----------------------------")
print("Electromagnetic pressure is a quintessential statistical mechanics phenomenon.")
print("The photon gas partition function in the canonical ensemble is:")
print()
print("  Z = Π_k (1 - exp(-βℏω_k))^(-2)")
print()
print("The factor 2 accounts for two polarization states. Taking ln Z and computing")
print("the free energy F = -k_BT ln Z, then using P = -∂F/∂V, we obtain:")
print()
print("  P = (π²/45) (k_BT/ℏc)³ k_BT = u/3")
print()
print("The factor 1/3 arises because photons contribute equally to pressure in all")
print("three spatial directions. For matter (v << c), particles bounce back from")
print("walls, giving P ~ ρv², but photons (v = c) are absorbed and re-emitted,")
print("transferring momentum 2ℏk, yielding P = u/3.")
print()
print("In cosmology, the radiation-dominated era (z > 3400) had equation of state")
print("w = P/ρ = 1/3, causing the scale factor to evolve as a(t) ~ t^(1/2). This is")
print("faster expansion than matter domination (a ~ t^(2/3)), affecting BBN and the")
print("CMB power spectrum. The partition function for the early universe includes:")
print("  • Photons (bosons, massless)")
print("  • Neutrinos (fermions, effectively massless at high T)")
print("  • e⁺e⁻ pairs (below T ~ 0.5 MeV, they annihilate)")
print()
print("The energy density evolves as ρ_rad ∝ (1+z)⁴ due to redshift (frequency drops")
print("as 1/(1+z)) plus volume expansion ((1+z)³). This scaling is encoded in the")
print("partition function's temperature dependence: Z ∝ (k_BT)³ for massless bosons,")
print("giving S ∝ T³ and u ∝ T⁴ (Stefan-Boltzmann).")
print()
print("TriPhase connects EM pressure to vacuum structure through ε_0 and μ_0, which")
print("determine c and Z_0. This suggests radiation pressure may be understood as")
print("vacuum response to energy density — a theme consistent with emergent spacetime")
print("and holographic principles.")
print()
print("=" * 70)

input("Press Enter to exit...")
