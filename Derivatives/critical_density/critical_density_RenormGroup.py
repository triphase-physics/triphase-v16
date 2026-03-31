"""
TriPhase V16 — Critical Density (Renormalization Group Framework)
========================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

RENORMALIZATION GROUP INTERPRETATION
--------------------------------------
The critical density ρ_crit = 3H₀²/(8πG) represents the total energy density
required for spatial flatness in an expanding universe. In the RG framework,
ρ_crit is the IR limit of the total mode density—the sum of all quantum field
excitations (matter, radiation, vacuum) at the cosmic scale after integrating
out all modes from UV (Planck) to IR (horizon).

In quantum cosmology, the Friedmann equation H² = (8πG/3)ρ - k/a² + Λ/3 relates
the expansion rate to the total energy density. For flat spatial geometry (k=0),
the critical density ρ_crit sets the boundary between open (ρ < ρ_crit) and
closed (ρ > ρ_crit) universes. Observations show Ω_tot = ρ_tot/ρ_crit ≈ 1.00
to within 1%, indicating that our universe is nearly spatially flat.

The TriPhase formula H₀ = π√3 × f_e × α¹⁸ connects ρ_crit to particle physics
via the electron mass m_e. The critical density is NOT a free parameter but
emerges from the α¹⁸ cascade: ρ_crit ∝ H₀² ∝ α³⁶ × (electron scale)². This
demonstrates that cosmological density scales are determined by the same RG
flow that generates particle masses, linking the largest and smallest scales
through a unified geometric structure.

TAG: (D) — Pure derivation; critical density from α¹⁸ cascade
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
print("TriPhase V16: Critical Density (RG Framework)")
print("=" * 70)
print()

print("HUBBLE PARAMETER FROM α¹⁸ CASCADE")
print("-" * 70)
print(f"Electron Compton frequency:      f_e = {f_e:.6e} Hz")
print(f"Fine structure constant:         α = {alpha:.10f}")
print(f"Geometric factor:                π√3 = {math.pi * math.sqrt(3):.10f}")
print()
print("TriPhase Hubble parameter:")
print("  H₀ = π√3 × f_e × α¹⁸")
print()
print(f"α¹⁸ suppression:                 α¹⁸ = {alpha**18:.6e}")
print()

H_0_calc = math.pi * math.sqrt(3) * f_e * alpha**18

print(f"Hubble parameter (TriPhase):     H₀ = {H_0_calc:.6e} s⁻¹")
print()

print("CRITICAL DENSITY FROM FRIEDMANN EQUATION")
print("-" * 70)
print("For a flat (k=0) universe, the Friedmann equation gives:")
print("  H² = (8πG/3) ρ_crit")
print()
print("Solving for critical density:")
print("  ρ_crit = 3H₀² / (8πG)")
print()
print(f"Newton's constant:               G = {G:.6e} m³/(kg·s²)")
print()

rho_crit = 3 * H_0_calc**2 / (8 * math.pi * G)

print(f"Critical density (TriPhase):     ρ_crit = {rho_crit:.6e} kg/m³")
print()

print("RG INTERPRETATION: TOTAL MODE DENSITY AT IR")
print("-" * 70)
print("The critical density represents the total quantum field mode density")
print("at the cosmic scale (IR endpoint of RG flow):")
print()
print("  ρ_crit = Σ_i ρ_i  (sum over all field species)")
print()
print("where ρ_i includes matter, radiation, and vacuum energy densities.")
print()
print("In the RG framework:")
print("  ρ_crit(μ) = ρ_crit(M_Pl) + ∫[M_Pl → H₀] β(ρ) d(ln μ)")
print()
print("The α¹⁸ cascade determines the IR endpoint ρ_crit(H₀).")
print()

print("SCALING WITH α³⁶")
print("-" * 70)
print("Since H₀ ∝ α¹⁸, the critical density scales as:")
print("  ρ_crit ∝ H₀² ∝ (α¹⁸)² = α³⁶")
print()
print(f"α³⁶ suppression factor:          α³⁶ = {alpha**36:.6e}")
print()
print("This extreme suppression (~10⁻⁶⁴) explains why cosmic densities are")
print("so small compared to particle physics scales (e.g., nuclear density")
print("ρ_nuclear ~ 10¹⁸ kg/m³).")
print()

rho_nuclear = 2.3e17  # kg/m³ (approximate nuclear density)
print(f"Nuclear density:                 ρ_nuclear ~ {rho_nuclear:.2e} kg/m³")
print(f"Critical density:                ρ_crit = {rho_crit:.2e} kg/m³")
print(f"Ratio ρ_nuclear / ρ_crit:        {rho_nuclear / rho_crit:.2e}")
print()

# ========== CALIBRATION CHECKPOINT ==========
H_0_CODATA = 2.197e-18  # s⁻¹ (67.4 km/s/Mpc, Planck 2018)
rho_crit_CODATA = 3 * H_0_CODATA**2 / (8 * math.pi * G)
deviation_pct = abs(rho_crit - rho_crit_CODATA) / rho_crit_CODATA * 100

print("CALIBRATION vs. CODATA/PLANCK")
print("-" * 70)
print(f"CODATA Hubble parameter:         H₀ = {H_0_CODATA:.6e} s⁻¹")
print(f"CODATA critical density:         ρ_crit = {rho_crit_CODATA:.6e} kg/m³")
print()
print(f"TriPhase Hubble parameter:       H₀ = {H_0_calc:.6e} s⁻¹")
print(f"TriPhase critical density:       ρ_crit = {rho_crit:.6e} kg/m³")
print(f"Deviation:                       {deviation_pct:.2f}%")
print()

print("DENSITY FRACTIONS (PLANCK 2018)")
print("-" * 70)
Omega_Lambda = 0.685  # Dark energy
Omega_m = 0.315       # Matter (baryonic + dark)
Omega_r = 9.2e-5      # Radiation (negligible today)

print(f"Dark energy:                     Ω_Λ = {Omega_Lambda:.3f}")
print(f"Matter:                          Ω_m = {Omega_m:.3f}")
print(f"Radiation:                       Ω_r = {Omega_r:.5f}")
print(f"Total:                           Ω_tot = {Omega_Lambda + Omega_m + Omega_r:.4f}")
print()
print("The observed Ω_tot ≈ 1.00 confirms that our universe is spatially flat,")
print("consistent with the critical density predicted by TriPhase.")
print()

print("RENORMALIZATION GROUP INSIGHT")
print("-" * 70)
print("Critical density is the IR fixed point of the total quantum field mode density.")
print("The α¹⁸ cascade from electron (UV) to cosmic (IR) scale determines H₀, which")
print("in turn fixes ρ_crit. The scaling ρ_crit ∝ α³⁶ explains why cosmic densities")
print("are so small—they're the ultimate IR remnant after 18 RG steps (squared) from")
print("particle physics scales. This connects the largest and smallest scales through")
print("a single unified RG structure, demonstrating that cosmology emerges from particle")
print("physics via deterministic flow, not arbitrary boundary conditions.")
print()
print("=" * 70)

input("Press Enter to exit...")
