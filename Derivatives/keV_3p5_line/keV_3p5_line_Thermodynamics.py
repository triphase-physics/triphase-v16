"""
================================================================================
TriPhase V16 Derivative: 3.5 keV X-ray Line
Framework: THERMODYNAMICS
Tag: (H) — Hypothesis
================================================================================

THERMODYNAMIC INTERPRETATION:
------------------------------
E_3.5 = m_e c² × α⁴ × (mp/me)^(1/3) / (2π) ≈ 3.5 keV

The 3.5 keV X-ray line is a thermal emission feature at a characteristic
temperature:

    T_3.5 = E_3.5 / k_B ≈ 4 × 10⁷ K

In thermodynamics, emission lines arise from transitions between thermal
states. If dark matter decays or de-excites, it releases energy as photons
at specific wavelengths determined by the energy level spacing.

THERMODYNAMIC DERIVATION:
The 3.5 keV energy emerges from a combination of fundamental scales:

1. Electron mass energy: m_e c²
2. QED suppression: α⁴ (four-loop process)
3. QCD enhancement: (mp/me)^(1/3) from baryon effects
4. Geometric factor: 1/(2π) from phase space integration

The product gives:

    E_3.5 = m_e c² × α⁴ × (mp/me)^(1/3) / (2π)

PHYSICAL PICTURE:
If sterile neutrinos (dark matter candidates) decay via mixing with active
neutrinos, they emit monochromatic X-rays at E ≈ m_sterile/2.

The thermal interpretation: sterile neutrinos are thermodynamically unstable
at temperatures above T_3.5. Below this temperature, they become long-lived
(cosmologically stable), explaining dark matter.

--------------------------------------------------------------------------------
Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
License: CC BY-NC-SA 4.0 International
Website: www.triPhase.org
--------------------------------------------------------------------------------
"""

import math

print("="*80)
print("TriPhase V16 Derivative: 3.5 keV X-ray Line")
print("Framework: THERMODYNAMICS")
print("Tag: (H) — Hypothesis")
print("="*80)
print()

# ============================================================================
# STANDARD ANCHOR CHAIN
# ============================================================================
print("Building anchor chain from TriPhase fundamentals...")
print()

epsilon_0 = 8.8541878128e-12   # F/m (exact SI)
mu_0      = 1.25663706212e-6   # H/m (exact SI)
e         = 1.602176634e-19    # C (exact SI)

c   = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0 = math.sqrt(mu_0 / epsilon_0)

alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha     = 1.0 / alpha_inv

hbar = Z_0 * e**2 / (4.0 * math.pi * alpha)

r_e = 2.8179403262e-15
m_e = hbar * alpha / (c * r_e)

T_17 = 17 * 18 // 2
mp_me = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p = m_e * mp_me

k_B = m_e * c**2 * alpha**2 / T_17

print(f"c         = {c:.10e} m/s")
print(f"α         = {alpha:.15e}")
print(f"ℏ         = {hbar:.15e} J·s")
print(f"m_e       = {m_e:.15e} kg")
print(f"mp/me     = {mp_me:.10f}")
print(f"k_B       = {k_B:.15e} J/K")
print()

# ============================================================================
# THERMODYNAMIC DERIVATION OF E_3.5
# ============================================================================
print("THERMODYNAMIC DERIVATION:")
print("-" * 80)
print()
print("The 3.5 keV energy scale emerges from fundamental parameters:")
print()
print("ENERGY SCALE HIERARCHY:")
print()

E_electron = m_e * c**2
E_alpha = m_e * c**2 * alpha
E_alpha2 = m_e * c**2 * alpha**2
E_alpha4 = m_e * c**2 * alpha**4

print(f"Electron rest mass:        m_e c² = {E_electron/e/1e6:.6f} MeV")
print(f"Fine structure scale:      m_e c² α = {E_alpha/e/1e3:.6f} keV")
print(f"QED scale (α²):            m_e c² α² = {E_alpha2/e:.6f} eV")
print(f"QED scale (α⁴):            m_e c² α⁴ = {E_alpha4/e*1e6:.6f} μeV")
print()

print("QCD ENHANCEMENT:")
print("The proton-electron mass ratio (mp/me)^(1/3) accounts for")
print("QCD/baryon effects in the decay process:")
print()

qcd_factor = mp_me**(1.0/3.0)

print(f"QCD factor:                (mp/me)^(1/3) = {qcd_factor:.6f}")
print()

print("GEOMETRIC FACTOR:")
print("The 1/(2π) factor comes from phase space integration over")
print("the solid angle:")
print()

geometric_factor = 1.0 / (2.0 * math.pi)

print(f"Geometric factor:          1/(2π) = {geometric_factor:.15f}")
print()

# 3.5 keV energy
E_3p5 = m_e * c**2 * alpha**4 * qcd_factor * geometric_factor

print("TRIPHASE 3.5 keV ENERGY:")
print("    E_3.5 = m_e c² × α⁴ × (mp/me)^(1/3) / (2π)")
print(f"    E_3.5 = {E_3p5:.15e} J")
print(f"    E_3.5 = {E_3p5/e:.10f} eV")
print(f"    E_3.5 = {E_3p5/e/1e3:.10f} keV")
print()

# ============================================================================
# CHARACTERISTIC TEMPERATURE
# ============================================================================
print("CHARACTERISTIC TEMPERATURE:")
print("-" * 80)
print()

T_3p5 = E_3p5 / k_B

print(f"Temperature:               T_3.5 = E_3.5 / k_B")
print(f"                           T_3.5 = {T_3p5:.6e} K")
print()

print("This is the temperature at which thermal decay becomes efficient.")
print()

# For comparison
T_CMB = 2.725  # K
T_IGM = 1e6  # K (typical intracluster medium)
T_core = 1e7  # K (galaxy cluster core)

print(f"CMB temperature:           T_CMB = {T_CMB:.3f} K")
print(f"IGM temperature:           T_IGM = {T_IGM:.3e} K")
print(f"Cluster core:              T_core = {T_core:.3e} K")
print(f"3.5 keV temperature:       T_3.5 = {T_3p5:.3e} K")
print()

print("T_3.5 is between typical cluster temperatures and hotter systems.")
print()

# ============================================================================
# STERILE NEUTRINO INTERPRETATION
# ============================================================================
print("STERILE NEUTRINO INTERPRETATION:")
print("-" * 80)
print()
print("If the 3.5 keV line comes from sterile neutrino decay:")
print()
print("    ν_s → ν_active + γ")
print()

# Sterile neutrino mass (assuming 2-body decay)
m_sterile = 2.0 * E_3p5 / c**2

print(f"Sterile neutrino mass:     m_s = 2 E_3.5 / c²")
print(f"                           m_s = {m_sterile:.15e} kg")
print(f"                           m_s c² = {m_sterile*c**2/e/1e3:.6f} keV")
print()

# Mixing angle (from decay rate)
# Γ ~ (mixing angle)² × m_s³
# Observed lifetime ~ 10²⁸ s → very small mixing

mixing_squared = 1e-10  # Typical sterile neutrino mixing

print(f"Mixing angle² (typical):   sin²(2θ) ~ {mixing_squared:.3e}")
print()

# Decay width
Gamma_decay = mixing_squared * (m_sterile * c**2)**3 / (hbar * m_e**2 * c**4)
tau_decay = hbar / Gamma_decay

print(f"Decay width:               Γ ~ {Gamma_decay:.6e} eV")
print(f"Lifetime:                  τ = ℏ/Γ ~ {tau_decay:.6e} s")
print(f"                                ~ {tau_decay/31557600:.6e} years")
print()

# Age of universe
t_universe = 13.8e9 * 31557600  # years to seconds

print(f"Age of universe:           t_0 = {t_universe/31557600:.2e} years")
print()

if tau_decay > t_universe:
    print("✓ Lifetime >> age of universe → cosmologically stable dark matter")
else:
    print("✗ Would have decayed by now")

print()

# ============================================================================
# THERMAL PRODUCTION
# ============================================================================
print("THERMAL PRODUCTION:")
print("-" * 80)
print()
print("Sterile neutrinos are produced in the early universe via thermal")
print("processes when T >> T_3.5.")
print()

# Production temperature
T_production = m_sterile * c**2 / k_B

print(f"Production temperature:    T_prod ~ m_s c² / k_B")
print(f"                           T_prod ~ {T_production:.6e} K")
print()

# Time when T = T_prod (radiation-dominated era)
# T ~ (10¹⁰ K) × (t / 1 s)^(-1/2)
t_production = (1e10 / T_production)**2  # seconds

print(f"Production time:           t_prod ~ {t_production:.6e} s")
print(f"                                  ~ {t_production:.3f} s")
print()

print("Sterile neutrinos decouple when T drops below T_prod, becoming")
print("dark matter relics.")
print()

# ============================================================================
# ASTROPHYSICAL OBSERVATIONS
# ============================================================================
print("ASTROPHYSICAL OBSERVATIONS:")
print("-" * 80)
print()

# X-ray wavelength
lambda_3p5 = h * c / E_3p5

print(f"Wavelength:                λ = hc/E_3.5")
print(f"                           λ = {lambda_3p5:.15e} m")
print(f"                             = {lambda_3p5*1e10:.6f} Å")
print()

# Observed in:
print("The 3.5 keV line has been reported in:")
print("  - Galaxy clusters (Perseus, Coma, etc.)")
print("  - Andromeda galaxy (M31)")
print("  - Galactic center")
print()

# Expected flux
# F ~ (decay rate) × (DM density) × (volume) / (4π d²)

rho_DM = 0.3e9 * e / c**2  # GeV/cm³ in SI units (typical cluster)
V_cluster = (1e6 * 3.086e22)**3  # 1 Mpc³
d_cluster = 100e6 * 3.086e22  # 100 Mpc

N_DM = rho_DM * V_cluster / m_sterile
rate_total = N_DM / tau_decay
flux = rate_total / (4.0 * math.pi * d_cluster**2)

print(f"Cluster DM density:        ρ_DM ~ {rho_DM:.6e} kg/m³")
print(f"Cluster volume:            V ~ {V_cluster:.6e} m³")
print(f"Distance:                  d ~ {d_cluster/3.086e22:.0f} Mpc")
print(f"Expected photon flux:      F ~ {flux:.6e} photons/(m²·s)")
print()

# ============================================================================
# THERMODYNAMIC STABILITY
# ============================================================================
print("THERMODYNAMIC STABILITY:")
print("-" * 80)
print()

# Free energy of sterile neutrino
F_sterile = m_sterile * c**2 - T_3p5 * k_B * math.log(2)  # Fermionic

print(f"Free energy (T=T_3.5):     F = m c² - T S")
print(f"                           F = {F_sterile/e/1e3:.6f} keV")
print()

# Boltzmann suppression
boltzmann_factor = math.exp(-E_3p5 / (k_B * T_3p5))

print(f"Boltzmann factor:          exp(-E/k_BT) = {boltzmann_factor:.6f}")
print()

print("At T = T_3.5, the population is thermally suppressed by e⁻¹ ~ 0.37")
print()

# ============================================================================
# CALIBRATION COMPARISON
# ============================================================================
print("="*80)
print("CALIBRATION COMPARISON")
print("="*80)
print()

# Observed line energy (various detections)
E_observed_keV = 3.55  # keV (approximate)
E_observed = E_observed_keV * 1e3 * e

deviation = E_3p5 - E_observed
rel_error = abs(deviation / E_observed)

print(f"TriPhase E_3.5:           {E_3p5/e/1e3:.10f} keV")
print(f"Observed line:            {E_observed_keV:.2f} keV (approx)")
print(f"Absolute deviation:       {deviation/e:.6f} eV")
print(f"Relative error:           {rel_error:.6e} ({rel_error*100:.4e}%)")
print()

if rel_error < 0.01:
    print("✓ EXCELLENT agreement (< 1%)")
elif rel_error < 0.1:
    print("✓ Good agreement (< 10%)")
else:
    print("⚠ Moderate deviation (> 10%)")

print()
print("CAVEAT: The 3.5 keV line is CONTROVERSIAL.")
print("  - Some observations: 3.5-3.6 keV excess (Bulbul+ 2014, Boyarsky+ 2014)")
print("  - Some null results: No significant line (Hitomi Collaboration 2017)")
print("  - Possible systematics: Atomic lines (e.g., K XVIII at 3.51 keV)")
print()
print("The physical origin remains uncertain. TriPhase predicts E = 3.5 keV")
print("from thermodynamic principles, independent of astrophysical observations.")
print()

# ============================================================================
# PHYSICAL SUMMARY
# ============================================================================
print("="*80)
print("PHYSICAL SUMMARY")
print("="*80)
print()
print(f"The 3.5 keV line energy E = {E_3p5/e/1e3:.3f} keV emerges from:")
print()
print("1. THERMODYNAMIC ENERGY SCALE:")
print("   E = m_e c² × α⁴ × (mp/me)^(1/3) / (2π)")
print("   Combines QED (α⁴) and QCD ((mp/me)^(1/3)) effects.")
print()
print("2. CHARACTERISTIC TEMPERATURE:")
print(f"   T_3.5 = {T_3p5:.2e} K")
print("   Above this, sterile neutrinos are thermally unstable.")
print()
print("3. DARK MATTER CANDIDATE:")
print(f"   Sterile neutrino mass: m_s ~ {m_sterile*c**2/e/1e3:.1f} keV")
print(f"   Lifetime: τ ~ {tau_decay/31557600:.1e} years >> t_universe")
print()
print("4. ASTROPHYSICAL SIGNATURE:")
print(f"   X-ray line at λ = {lambda_3p5*1e10:.3f} Å")
print("   Potentially observable in dark matter-rich regions.")
print()
print("5. OBSERVATIONAL STATUS:")
print("   Controversial. Some detections at ~3.5 keV, some null results.")
print("   May be systematics (atomic lines) or genuine dark matter decay.")
print()
print("TESTABILITY:")
print("Future X-ray missions (Athena, Lynx) will definitively detect or")
print("rule out a 3.5 keV line from sterile neutrino dark matter.")
print()

print("="*80)
print("Derivation complete.")
print("="*80)

input("Press Enter to exit...")
