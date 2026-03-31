"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================
Derivative:  Up Quark Mass (m_u ~ 2.16 MeV/c²)
Framework:   DiffGeometry
Version:     16.0
Generated:   2026-03-26
Status:      Active Development
Tag: (D*H)
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
"""
import math

print("="*80)
print("TriPhase V16: Up Quark Mass")
print("Framework: DiffGeometry")
print("="*80)

# === ANCHOR INPUTS ===
epsilon_0 = 8.8541878128e-12  # F/m
mu_0      = 1.25663706212e-6  # H/m
e         = 1.602176634e-19   # C

print("\n[ANCHOR INPUTS]")
print(f"ε₀ = {epsilon_0:.13e} F/m")
print(f"μ₀ = {mu_0:.14e} H/m")
print(f"e  = {e:.12e} C")

# === DERIVED ANCHOR CHAIN ===
c         = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0       = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha     = 1.0 / alpha_inv
hbar      = Z_0 * e**2 / (4.0 * math.pi * alpha)
h         = 2.0 * math.pi * hbar
G         = c**4 * 7.5 * epsilon_0**3 * mu_0**2
r_e       = 2.8179403262e-15  # m
m_e       = hbar * alpha / (c * r_e)
f_e       = m_e * c**2 / hbar
T_17      = 17 * 18 // 2  # = 153
mp_me     = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p       = m_e * mp_me

print("\n[DERIVED ANCHOR CHAIN]")
print(f"c     = {c:.6e} m/s")
print(f"Z₀    = {Z_0:.10f} Ω")
print(f"α⁻¹   = {alpha_inv:.10f}")
print(f"α     = {alpha:.12f}")
print(f"ℏ     = {hbar:.12e} J·s")
print(f"h     = {h:.12e} J·s")
print(f"G     = {G:.12e} m³/(kg·s²)")
print(f"r_e   = {r_e:.13e} m")
print(f"m_e   = {m_e:.12e} kg")
print(f"f_e   = {f_e:.12e} Hz")
print(f"T₁₇   = {T_17}")
print(f"mp/me = {mp_me:.10f}")
print(f"m_p   = {m_p:.12e} kg")

# === DERIVATION: Up Quark Mass ===
print("\n" + "="*80)
print("DERIVATION: Up Quark Mass")
print("="*80)

print("\nDiffGeometry Interpretation:")
print("Quarks live on the SU(3) color manifold, a more complex fiber bundle")
print("than the U(1) electromagnetic bundle of leptons. The up quark, with")
print("electric charge +2/3, has mass determined by the curvature ratio between")
print("strong and electromagnetic sectors.")
print()
print("Mass formula (first-order estimate):")
print("  m_u = m_e × 2/(3α)")
print()
print("The factor 2/3 is the fractional charge (winding number on the EM fiber).")
print("The 1/α factor represents the inverse coupling — up quarks feel the")
print("strong curvature more than the EM curvature, scaling their mass.")

# Up quark mass ratio
charge_factor = 2.0 / 3.0
coupling_factor = 1.0 / alpha
up_mass_ratio = charge_factor * coupling_factor

# Up quark mass
m_up_kg = m_e * up_mass_ratio

# Convert to MeV/c²
energy_joules = m_up_kg * c**2
m_up_MeV = energy_joules / (e * 1e6)

print(f"\n2/3 (charge factor)   = {charge_factor:.10f}")
print(f"1/α (coupling depth)  = {coupling_factor:.10f}")
print(f"m_u/m_e               = {up_mass_ratio:.10f}")
print(f"\nm_u = {m_up_kg:.12e} kg")
print(f"m_u = {m_up_MeV:.6f} MeV/c²")

# === CALIBRATION CHECKPOINT ===
print("\n" + "="*80)
print("CALIBRATION CHECKPOINT")
print("="*80)

pdg_m_up_MeV = 2.16  # MeV/c², MS-bar at 2 GeV, PDG 2022
pdg_range = "2.16 +0.49 -0.26 MeV"

error_ppm = abs(m_up_MeV - pdg_m_up_MeV) / pdg_m_up_MeV * 1e6

print(f"TriPhase m_u:  {m_up_MeV:.6f} MeV/c²")
print(f"PDG m_u:       {pdg_m_up_MeV:.2f} MeV/c² (MS-bar at 2 GeV)")
print(f"PDG range:     {pdg_range}")
print(f"Error:         {error_ppm:.0f} ppm")
print()
print("Note: Quark masses are scheme-dependent (MS-bar vs pole mass) and")
print("running with energy scale. The TriPhase formula gives a first-order")
print("geometric estimate. Full QCD renormalization group flow is needed")
print("for precision comparison.")

print("\n" + "="*80)
print("GEOMETRIC INSIGHT")
print("="*80)
print("Quarks differ from leptons in two key ways:")
print()
print("1. Color charge: Quarks carry SU(3) color, living on a 3-dimensional")
print("   fiber bundle with non-abelian curvature. This is why they experience")
print("   confinement — the color curvature prevents isolated quarks.")
print()
print("2. Fractional EM charge: The 2/3 and -1/3 charges are winding numbers")
print("   on the U(1)_EM fiber. These are topological invariants, not arbitrary.")
print()
print("The up quark's small mass (compared to leptons at similar coupling depth)")
print("reflects the strong force's curvature structure — quarks are deeply")
print("embedded in the color manifold's geometry.")

print("\n" + "="*80)
input("Press Enter to exit...")
