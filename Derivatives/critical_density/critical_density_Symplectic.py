"""
TriPhase V16 — Critical Density (Symplectic Framework)

Copyright © 2025 MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
MIS TAG: (D)

SYMPLECTIC INTERPRETATION:
Critical density ρ_crit = 3H²/(8πG) is the energy density at which the universe's
cosmological phase space (a, ȧ) has exactly zero total energy—the boundary between
open (hyperbolic), flat (Euclidean), and closed (spherical) geometries. In the
ADM Hamiltonian formulation, ρ_crit emerges from the Hamiltonian constraint
H = 0, which relates the kinetic energy (ȧ)² to the potential energy (curvature k).

The critical density defines the flat symplectic manifold where spatial curvature
k = 0. The Friedmann equation (ȧ/a)² = (8πG/3)ρ - k/a² reduces to H² = (8πG/3)ρ_crit,
making ρ_crit the canonical scale for cosmological energy density. The density
parameter Ω = ρ/ρ_crit determines the phase space trajectory: Ω < 1 (open),
Ω = 1 (flat), Ω > 1 (closed).
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
print("TRIPHASE V16 — CRITICAL DENSITY (SYMPLECTIC FRAMEWORK)")
print("=" * 70)
print()

# ========== SYMPLECTIC DERIVATION ==========
print("PHASE SPACE STRUCTURE:")
print("  Cosmological phase space: (a, ȧ) where a(t) = scale factor")
print("  Hubble parameter: H = ȧ/a (expansion rate)")
print("  Conjugate momentum: p_a = -3ȧa²/(8πG) (from ADM decomposition)")
print("  Symplectic 2-form: ω = dp_a ∧ da")
print()

print("HAMILTONIAN FORMULATION:")
print("  Friedmann equation (Hamiltonian constraint):")
print("    H² = (ȧ/a)² = (8πG/3)ρ - k/a²")
print("  For flat universe (k=0): H² = (8πG/3)ρ_crit")
print("  Critical density: ρ_crit = 3H²/(8πG)")
print("  Density parameter: Ω_i = ρ_i/ρ_crit")
print("  Total: Ω_total = Ω_m + Ω_r + Ω_Λ + Ω_k")
print("    Ω_k = -k/(a²H²) (curvature contribution)")
print()

print("SYMPLECTIC INVARIANT:")
print("  The Friedmann equation is a Hamiltonian constraint H = 0")
print("  ρ_crit is the energy density where kinetic = potential (critical point)")
print("  Phase space trajectory type determined by Ω:")
print("    Ω < 1: hyperbolic (open), a(t) → ∞ asymptotically")
print("    Ω = 1: flat (critical), a(t) ~ t^{2/3} (matter-dominated)")
print("    Ω > 1: elliptic (closed), a(t) reaches maximum and recollapses")
print()

print("TRIPHASE DERIVATION:")
print("  Formula: ρ_crit = 3 H_0² / (8π G)")
print("  Using H_0 from TriPhase anchor chain:")
print()

# Compute critical density
rho_crit = 3.0 * H_0**2 / (8.0 * math.pi * G)
print(f"  H_0      = {H_0:.6e} Hz")
print(f"  G        = {G:.6e} m³/(kg·s²)")
print(f"  ρ_crit   = 3H_0²/(8πG)")
print(f"           = {rho_crit:.6e} kg/m³")
print()

# Convert to other units
rho_crit_GeV = rho_crit * c**2 / (1.602176634e-10)  # GeV/cm³
print(f"  ρ_crit   = {rho_crit_GeV:.6e} GeV/cm³")
print(f"  ρ_crit   = {rho_crit * c**2:.6e} J/m³")
print()

# Compute Hubble time and distance
t_H = 1.0 / H_0
r_H = c / H_0
print(f"  Hubble time: t_H = 1/H_0 = {t_H:.6e} s")
print(f"                            = {t_H / (365.25 * 24 * 3600):.3e} years")
print(f"  Hubble distance: r_H = c/H_0 = {r_H:.6e} m")
print(f"                                = {r_H / 9.461e15:.3e} light-years")
print()

# Compute age of universe (flat, matter-dominated approximation)
# t_0 = (2/3) t_H for Ω_m = 1, w = 0
t_universe_flat = (2.0 / 3.0) * t_H
print(f"  Age of flat, matter-dominated universe:")
print(f"  t_0 = (2/3) t_H = {t_universe_flat:.6e} s")
print(f"                  = {t_universe_flat / (365.25 * 24 * 3600 * 1e9):.2f} Gyr")
print()

# ========== CALIBRATION CHECKPOINT ==========
print("CALIBRATION CHECKPOINT:")
print("  Using Planck 2018 + H_0 = 67.4 km/s/Mpc:")
print("    H_0 = 2.184e-18 Hz")
print("    ρ_crit = 8.6e-27 kg/m³")
print("    Age of universe = 13.8 Gyr")
print()
H_0_measured = 2.184e-18  # Hz
rho_crit_measured = 8.6e-27  # kg/m³
t_universe_measured = 13.8e9 * 365.25 * 24 * 3600  # seconds
print(f"  Measured H_0        = {H_0_measured:.6e} Hz")
print(f"  TriPhase H_0        = {H_0:.6e} Hz")
deviation_H_ppm = abs(H_0 - H_0_measured) / H_0_measured * 1e6
print(f"  Deviation: {deviation_H_ppm:.1f} ppm")
print()
print(f"  Measured ρ_crit     = {rho_crit_measured:.6e} kg/m³")
print(f"  TriPhase ρ_crit     = {rho_crit:.6e} kg/m³")
deviation_rho_ppm = abs(rho_crit - rho_crit_measured) / rho_crit_measured * 1e6
print(f"  Deviation: {deviation_rho_ppm:.1f} ppm")
print()

# Compare to other fundamental densities
rho_nuclear = 2.3e17  # kg/m³ (nuclear density)
rho_water = 1000.0  # kg/m³
print(f"  Comparison:")
print(f"  ρ_crit / ρ_water   = {rho_crit / rho_water:.6e}")
print(f"  ρ_crit / ρ_nuclear = {rho_crit / rho_nuclear:.6e}")
print()

# ========== SYMPLECTIC GEOMETRY INSIGHT ==========
print("SYMPLECTIC GEOMETRY INSIGHT:")
print("  Critical density is the bifurcation point in cosmological phase space.")
print("  At ρ = ρ_crit, the Hamiltonian constraint H = 0 defines a flat")
print("  symplectic manifold with k = 0 spatial curvature. The phase space")
print("  trajectory (a, ȧ) lies exactly on the separatrix between bounded")
print("  (closed) and unbounded (open) orbits.")
print()
print("  In the (ρ, H) phase diagram, ρ_crit is the critical point where the")
print("  system transitions from one topological class to another. Planck 2018")
print("  measurements show Ω_total = 1.000 ± 0.004, confirming the universe")
print("  sits precisely on this critical symplectic surface—a remarkable")
print("  fine-tuning (the flatness problem).")
print()
print("  ρ_crit sets the natural energy scale for cosmology: it's the energy")
print("  density where gravitational dynamics transitions from local (Newtonian)")
print("  to global (cosmological). All cosmological observables (Ω_m, Ω_Λ, etc.)")
print("  are referenced to this critical symplectic manifold.")
print("=" * 70)

input("Press Enter to exit...")
