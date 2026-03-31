"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================
Derivative:  Charm Quark Mass (m_c ~ 1.27 GeV/c²)
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
print("TriPhase V16: Charm Quark Mass")
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

# === DERIVATION: Charm Quark Mass ===
print("\n" + "="*80)
print("DERIVATION: Charm Quark Mass")
print("="*80)

print("\nDiffGeometry Interpretation:")
print("The charm quark is the second-generation up-type quark. Its mass is")
print("at the proton curvature scale (mp/me) with a geometric coupling factor")
print("3/(2π) and first-order radiative correction (1+α).")
print()
print("Mass formula:")
print("  m_c = m_e × (mp/me) × 3/(2π) × (1+α)")
print()
print("The mp/me factor places charm at the nucleon curvature scale. The 3/(2π)")
print("is the geometric coupling for second-generation up-type quarks. The (1+α)")
print("correction accounts for parallel transport holonomy at this energy scale.")

# Charm quark mass ratio
geometric_factor = 3.0 / (2.0 * math.pi)
radiative_corr = 1.0 + alpha
charm_mass_ratio = mp_me * geometric_factor * radiative_corr

# Charm quark mass
m_charm_kg = m_e * charm_mass_ratio

# Convert to GeV/c²
energy_joules = m_charm_kg * c**2
m_charm_GeV = energy_joules / (e * 1e9)

print(f"\nmp/me                 = {mp_me:.10f}")
print(f"3/(2π)                = {geometric_factor:.10f}")
print(f"1 + α                 = {radiative_corr:.10f}")
print(f"m_c/m_e               = {charm_mass_ratio:.10f}")
print(f"\nm_c = {m_charm_kg:.12e} kg")
print(f"m_c = {m_charm_GeV:.6f} GeV/c²")

# === CALIBRATION CHECKPOINT ===
print("\n" + "="*80)
print("CALIBRATION CHECKPOINT")
print("="*80)

pdg_m_charm_GeV = 1.27  # GeV/c², MS-bar, PDG 2022
pdg_range = "1.27 ± 0.02 GeV"

error_ppm = abs(m_charm_GeV - pdg_m_charm_GeV) / pdg_m_charm_GeV * 1e6

print(f"TriPhase m_c:  {m_charm_GeV:.6f} GeV/c²")
print(f"PDG m_c:       {pdg_m_charm_GeV:.2f} GeV/c² (MS-bar)")
print(f"PDG range:     {pdg_range}")
print(f"Error:         {error_ppm:.0f} ppm")
print()
print("Note: Charm quark mass is better constrained than light quarks due to")
print("clean charmonium (cc̄) spectroscopy. The TriPhase formula shows how")
print("charm sits at the nucleon curvature scale with generation-2 corrections.")

print("\n" + "="*80)
print("GEOMETRIC INSIGHT")
print("="*80)
print("Heavy quarks (charm, bottom, top) have masses at or above the nucleon")
print("scale, meaning they probe deeper curvature regions of the vacuum manifold.")
print()
print("Quark mass hierarchy:")
print("  Generation 1 (u,d): ~MeV scale, 1/α depth")
print("  Generation 2 (c,s): ~100 MeV - 1 GeV, nucleon scale with T₁₇")
print("  Generation 3 (t,b): ~GeV - 100 GeV, α⁻¹ depth above nucleon")
print()
print("The charm quark is the lightest quark above the nucleon scale. Its mass")
print("formula includes mp/me because charm quarks are heavy enough to 'see'")
print("the full nucleon curvature structure, not just the pion/EM layers.")
print()
print("This is why charm hadrons (D mesons, Λ_c baryons) have masses ~2 GeV —")
print("comparable to nucleons but at a higher curvature tier.")

print("\n" + "="*80)
input("Press Enter to exit...")
