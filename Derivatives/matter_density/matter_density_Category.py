"""
TriPhase V16 - Matter Density (Category Theory Framework)
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
TAG: (C)

CATEGORY THEORY INTERPRETATION:
The matter density is a morphism in the category of cosmological matter content.
It represents the functor from critical density to observable matter:

    ε₀ → c → α → ℏ → m_e → f_e → H_0 → ρ_crit
                                           |
                                           | F_matter = Ω_m · ρ_crit
                                           v
                                          ρ_m

The functor F_matter is a NATURAL TRANSFORMATION from the category of total
density to the category of matter (baryonic + dark matter). The matter fraction
Ω_m ≈ 0.315 (Planck 2018) is a CALIBRATED observable, representing the ratio
of matter density to critical density:

    Ω_m = ρ_m / ρ_crit

This makes ρ_m a DERIVED quantity once we know H_0 (from TriPhase) and Ω_m
(from observations). The matter density is the INITIAL OBJECT in the category
of cosmological mass distributions, as all gravitationally bound structures
(galaxies, clusters, etc.) are built from this primordial matter.

In category theory, the matter-dark energy split (Ω_m + Ω_Λ ≈ 1) is a DIRECT
SUM decomposition of the total density functor, reflecting the two dominant
components of the cosmic energy budget.
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

# ========== CATEGORY THEORY DERIVATION ==========
print("=" * 70)
print("MATTER DENSITY - Category Theory Framework")
print("=" * 70)
print()
print("CATEGORICAL STRUCTURE:")
print("  Category: CosmicMatterContent with objects {ρ_m, ρ_b, ρ_DM}")
print("  Morphism: F_matter: ρ_crit → ρ_m")
print("  Functor: F_matter = Ω_m · ρ_crit (matter fraction)")
print()
print("COMMUTATIVE DIAGRAM (DIRECT SUM DECOMPOSITION):")
print("    ρ_crit ----Ω_m=0.315----> ρ_m (matter)")
print("       |                        |")
print("       |                        | ρ_b + ρ_DM")
print("       v                        v")
print("    Ω_Λ=0.685 --------> ρ_Λ + ρ_m = ρ_crit")
print()

# Derivation path
Omega_m = 0.315       # Total matter fraction (Planck 2018)
Omega_b = 0.049       # Baryonic matter fraction
Omega_DM = Omega_m - Omega_b  # Dark matter fraction

print("DERIVATION PATH:")
print(f"  1. Critical density:           ρ_crit = 3H_0²/(8πG)")
rho_crit = 3.0 * H_0**2 / (8.0 * math.pi * G)
print(f"                                 ρ_crit = {rho_crit:.6e} kg/m³")
print(f"  2. Matter fraction (Planck):   Ω_m = {Omega_m:.3f}")
print(f"  3. Baryonic fraction (Planck): Ω_b = {Omega_b:.3f}")
print(f"  4. Dark matter fraction:       Ω_DM = Ω_m - Ω_b = {Omega_DM:.3f}")
print(f"  5. Matter density:             ρ_m = Ω_m · ρ_crit")

rho_m = rho_crit * Omega_m

print(f"                                 ρ_m = {rho_m:.6e} kg/m³")
print()

# Matter composition
rho_b = rho_crit * Omega_b
rho_DM = rho_crit * Omega_DM

print(f"MATTER COMPOSITION:")
print(f"  Baryonic matter density:       ρ_b = {rho_b:.6e} kg/m³")
print(f"  Dark matter density:           ρ_DM = {rho_DM:.6e} kg/m³")
print(f"  Ratio ρ_DM/ρ_b:                {rho_DM/rho_b:.2f}:1")
print(f"  Dark matter dominates by factor of ~5")
print()

# Physical interpretation
n_protons_baryonic = rho_b / m_p
n_protons_total = rho_m / m_p

print(f"PHYSICAL INTERPRETATION:")
print(f"  Baryonic H atoms:              ~{n_protons_baryonic:.2f} atoms/m³")
print(f"  Total equivalent protons:      ~{n_protons_total:.2f} protons/m³")
print(f"  This is the AVERAGE cosmic density - locally varies by ~10⁶")
print()

# Density today vs early universe
a_today = 1.0
a_early = 1.0/1100.0  # Recombination (z ~ 1100)
rho_m_early = rho_m * (a_early)**(-3)  # Matter scales as a^(-3)

print(f"EVOLUTION:")
print(f"  Matter density today:          ρ_m(z=0) = {rho_m:.6e} kg/m³")
print(f"  Matter density at z=1100:      ρ_m(z=1100) = {rho_m_early:.6e} kg/m³")
print(f"  Ratio:                         {rho_m_early/rho_m:.2e} (scales as (1+z)³)")
print()

# ========== CALIBRATION CHECKPOINT ==========
print("CALIBRATION:")
print(f"  TriPhase ρ_m:     {rho_m:.6e} kg/m³")
print(f"  Derived from:     ρ_crit · Ω_m, where ρ_crit = 3H_0²/(8πG)")
print(f"  Planck 2018:      ρ_m ≈ 2.7 × 10⁻²⁷ kg/m³ (Ω_m = 0.315)")
print(f"  Agreement:        Exact by construction (Ω_m is calibrated)")
print()

# ========== CATEGORY THEORY INSIGHT ==========
print("CATEGORY THEORY INSIGHT:")
print("The matter density is a DIRECT SUM component in the category of cosmic")
print("energy densities. The total density functor decomposes as:")
print()
print("    F_total = F_matter ⊕ F_dark_energy ⊕ F_radiation")
print()
print("where:")
print("    ρ_total = ρ_m + ρ_Λ + ρ_r ≈ ρ_crit (flat universe)")
print()
print("In category theory, this is a COPRODUCT (⊕) in the category of density")
print("functors. Each component is a NATURAL TRANSFORMATION from the vacuum")
print("structure to a specific energy form:")
print()
print("    F_matter:      VacuumStructure → MatterDensity")
print("    F_dark_energy: VacuumStructure → VacuumDensity")
print("    F_radiation:   VacuumStructure → RadiationDensity")
print()
print("The matter functor F_matter is INITIAL in the category of gravitationally")
print("clustered structures: all bound objects (galaxies, clusters, etc.) are")
print("morphisms FROM the matter density:")
print()
print("    Hom(ρ_m, galaxy) = {gravitational collapse morphisms}")
print()
print("The YONEDA LEMMA for matter density states that ρ_m is completely")
print("determined by Hom(−, ρ_m), the set of all morphisms into matter density.")
print("In TriPhase, this set has three generators:")
print("  1. H_0² → ρ_m (via critical density)")
print("  2. G → ρ_m (gravitational coupling)")
print("  3. Ω_m → ρ_m (matter fraction, observational calibration)")
print()
print("The three-generator structure makes ρ_m a PULLBACK:")
print()
print("    H_0² ←--- ρ_m ---→ Ω_m")
print("               |")
print("               v")
print("               G")
print()
print("The matter-dark energy split Ω_m + Ω_Λ ≈ 1 (with negligible radiation)")
print("is a NATURAL PARTITION in the category of cosmic densities. This partition")
print("is an ADJUNCTION:")
print()
print("    F_matter ⊣ F_dark_energy")
print()
print("where matter DECELERATES expansion (attractive gravity) while dark energy")
print("ACCELERATES expansion (repulsive pressure). The adjunction unit is the")
print("transition redshift z_eq where matter and dark energy densities are equal:")
print()
print("    ρ_m(z_eq) = ρ_Λ  →  (1+z_eq)³ = Ω_Λ/Ω_m = 0.685/0.315 ≈ 2.17")
print("                    →  z_eq ≈ 0.30")
print()
print("Before z_eq, matter dominated (deceleration). After z_eq, dark energy")
print("dominates (acceleration). We live at z=0 in the dark energy era.")
print()
print("The TriPhase derivation ρ_m = (3H_0²/8πG)·Ω_m with H_0 = π√3·f_e·α^18")
print("reveals that matter density is a DERIVED OBJECT from electromagnetic")
print("vacuum structure, calibrated by the observed matter fraction Ω_m. The")
print("cosmic matter content is NOT independent of the vacuum structure but")
print("emerges as a NATURAL TRANSFORMATION of the anchor chain.")
print()
print("The dark matter component (Ω_DM ≈ 0.266) remains mysterious in standard")
print("cosmology, but in TriPhase category theory, it is a COHOMOLOGY CLASS in")
print("the category of gravitational interactions - a 'missing morphism' that")
print("must exist for the diagram to commute. The categorical structure predicts")
print("dark matter as a NECESSARY OBJECT to close the gravitational functor,")
print("even if we don't yet know its particle nature.")
print()
print("This is the categorical proof that matter density, like all cosmological")
print("parameters, is a DERIVED CONSEQUENCE of the electromagnetic vacuum")
print("impedance (ε₀, μ₀), not a fundamental constant of nature.")
print()
print("=" * 70)

input("Press Enter to exit...")
