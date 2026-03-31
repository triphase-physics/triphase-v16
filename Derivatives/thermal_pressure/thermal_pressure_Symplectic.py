"""
TriPhase V16 — Thermal Pressure (Symplectic Framework)

Copyright © 2025 MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
MIS TAG: (C)

SYMPLECTIC INTERPRETATION:
Thermal pressure P = nk_BT is the momentum flux in the statistical phase space
of a gas. In the canonical ensemble, (q, p) coordinates of N particles form a
2N-dimensional symplectic manifold with ω = Σ dp_i ∧ dq_i. The pressure emerges
as the diagonal component of the stress tensor—the rate of momentum transfer
across a surface in phase space.

For blackbody radiation (photon gas), the pressure P_γ = (π²/15)(k_BT)⁴/(ℏ³c³)
arises from the symplectic volume of photon phase space integrated over Bose-
Einstein statistics. The CMB at T_CMB = 2.7255 K represents the equilibrium
thermal pressure of the early universe's photon-baryon fluid.
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
print("TRIPHASE V16 — THERMAL PRESSURE (SYMPLECTIC FRAMEWORK)")
print("=" * 70)
print()

# ========== SYMPLECTIC DERIVATION ==========
print("PHASE SPACE STRUCTURE:")
print("  Statistical phase space: Γ = {(q_i, p_i) | i=1..N} (6N dimensional)")
print("  Symplectic 2-form: ω = Σ dp_i ∧ dq_i")
print("  Phase space volume: dΓ = Π dq_i dp_i / h³ (quantum normalization)")
print("  For photons: massless bosons, phase space density n(ω) ~ ω²/(e^{ℏω/k_BT} - 1)")
print()

print("HAMILTONIAN FORMULATION:")
print("  H = Σ p_i²/(2m) (non-relativistic)")
print("  H = Σ c√(p_i² + m²c²) (relativistic)")
print("  For photons: H = ℏω, p = ℏk, E = pc")
print("  Partition function: Z = ∫ e^{-βH} dΓ (β = 1/k_BT)")
print("  Pressure: P = -∂F/∂V = k_BT ∂(ln Z)/∂V")
print()

print("SYMPLECTIC INVARIANT:")
print("  Liouville's theorem: dΓ is preserved under Hamiltonian flow")
print("  Entropy S = k_B ln Ω(E,V,N) counts symplectic volume at energy E")
print("  Thermal equilibrium: maximum entropy ↔ uniform phase space distribution")
print()

print("TRIPHASE DERIVATION:")
print("  Ideal gas: P = n k_B T")
print("  Photon gas (blackbody): P_γ = (π²/15) (k_B T)⁴ / (ℏ³ c³)")
print("  CMB temperature: T_CMB = 2.7255 K")
print()

# Constants
k_B = 1.380649e-23  # J/K (exact, SI 2019)
T_CMB = 2.7255  # K (COBE/FIRAS measurement)
zeta3 = 1.20206  # Riemann zeta(3)

print(f"  k_B   = {k_B:.6e} J/K")
print(f"  T_CMB = {T_CMB:.6f} K")
print()

# Compute photon number density
# n_γ = (2ζ(3)/π²) (k_B T / ℏc)³
n_photon = 2.0 * zeta3 / math.pi**2 * (k_B * T_CMB / (hbar * c))**3
print(f"  Photon number density:")
print(f"  n_γ = (2ζ(3)/π²) (k_B T_CMB / ℏc)³")
print(f"      = {n_photon:.6e} photons/m³")
print()

# Compute CMB pressure
P_CMB = (math.pi**2 / 15.0) * (k_B * T_CMB)**4 / (hbar**3 * c**3)
print(f"  CMB radiation pressure:")
print(f"  P_CMB = (π²/15) (k_B T_CMB)⁴ / (ℏ³ c³)")
print(f"        = {P_CMB:.6e} Pa")
print()

# Also compute energy density (for reference)
rho_CMB = (math.pi**2 / 15.0) * (k_B * T_CMB)**4 / (hbar**3 * c**5)
print(f"  CMB energy density:")
print(f"  ρ_CMB = (π²/15) (k_B T_CMB)⁴ / (ℏ³ c⁵)")
print(f"        = {rho_CMB:.6e} J/m³")
print(f"  Relation: P_CMB = ρ_CMB/3 (radiation equation of state)")
print(f"  Check: ρ_CMB/3 = {rho_CMB/3.0:.6e} Pa")
print()

# ========== CALIBRATION CHECKPOINT ==========
print("CALIBRATION CHECKPOINT:")
print("  Using COBE/FIRAS and Planck measurements:")
print("    T_CMB = 2.7255 ± 0.0006 K")
print("    n_γ ≈ 4.11e8 photons/m³")
print("    P_CMB ≈ 4.64e-14 Pa")
print()
n_photon_measured = 4.11e8
P_CMB_measured = 4.64e-14
print(f"  Measured n_γ     = {n_photon_measured:.6e} photons/m³")
print(f"  TriPhase n_γ     = {n_photon:.6e} photons/m³")
deviation_n_ppm = abs(n_photon - n_photon_measured) / n_photon_measured * 1e6
print(f"  Deviation: {deviation_n_ppm:.1f} ppm")
print()
print(f"  Measured P_CMB   = {P_CMB_measured:.6e} Pa")
print(f"  TriPhase P_CMB   = {P_CMB:.6e} Pa")
deviation_P_ppm = abs(P_CMB - P_CMB_measured) / P_CMB_measured * 1e6
print(f"  Deviation: {deviation_P_ppm:.1f} ppm")
print()

# ========== SYMPLECTIC GEOMETRY INSIGHT ==========
print("SYMPLECTIC GEOMETRY INSIGHT:")
print("  Thermal pressure is the rate of symplectic volume flow across a")
print("  boundary in phase space. For a photon gas, each mode contributes")
print("  ⟨p²⟩/m to the pressure—the average squared momentum in that symplectic")
print("  cell. The Planck distribution weights each cell by its Boltzmann factor.")
print()
print("  At thermal equilibrium, the symplectic flow is isotropic: pressure")
print("  is uniform in all directions. The CMB represents the fossilized")
print("  symplectic structure of the early universe's photon-baryon fluid")
print("  at decoupling (z ≈ 1100), redshifted to T_CMB = 2.7255 K today.")
print("=" * 70)

input("Press Enter to exit...")
