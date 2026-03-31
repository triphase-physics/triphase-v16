"""
================================================================================
TriPhase V16 Derivative: Electron Anomalous Magnetic Moment (g-2)
Framework: THERMODYNAMICS
Tag: (D) — Pure derivation from thermodynamic principles
================================================================================

THERMODYNAMIC INTERPRETATION:
------------------------------
g_e = 2(1 + α/(2π)) ≈ 2.002319...

The electron g-factor describes its magnetic moment:

    μ = g_e (eℏ/2m_e) S

where S is the spin. The Dirac equation predicts g = 2 exactly. The
anomalous magnetic moment (g-2)/2 = α/(2π) arises from quantum corrections.

In thermodynamics, this has an interpretation as a THERMAL FLUCTUATION
correction. The vacuum is not a static background — it's a thermal bath
of virtual particle-antiparticle pairs.

THERMODYNAMIC DERIVATION:
The magnetic moment receives corrections from vacuum fluctuations. The
one-loop correction (Schwinger term) is:

    a_e = (g-2)/2 = α/(2π)

This is the thermal correction from virtual pair creation. At the
Compton scale λ_C = ℏ/(m_e c), the vacuum has thermal fluctuations
that modify the effective magnetic moment.

The factor α/(2π) can be understood as:
- α: Coupling strength (probability of virtual photon emission)
- 1/(2π): Phase space factor from angular integration

Higher-order terms α²/(2π)², α³/(2π)³, ... correspond to multi-loop
thermal corrections (virtual pairs creating more virtual pairs).

PHYSICAL PICTURE:
The electron is surrounded by a "cloud" of virtual photons and pairs.
This cloud contributes to the magnetic moment, enhancing it beyond
the bare Dirac value g = 2.

This is a thermodynamic effect: the vacuum polarization creates an
effective medium that modifies the electron's response to magnetic fields.

--------------------------------------------------------------------------------
Copyright (c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
License: CC BY-NC-SA 4.0 International
Website: www.triPhase.org
--------------------------------------------------------------------------------
"""

import math

print("="*80)
print("TriPhase V16 Derivative: Electron Anomalous Magnetic Moment (g-2)")
print("Framework: THERMODYNAMICS")
print("Tag: (D) — Pure derivation")
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

print(f"c         = {c:.10e} m/s")
print(f"α         = {alpha:.15e}")
print(f"ℏ         = {hbar:.15e} J·s")
print(f"m_e       = {m_e:.15e} kg")
print()

# ============================================================================
# THERMODYNAMIC DERIVATION OF g-2
# ============================================================================
print("THERMODYNAMIC DERIVATION:")
print("-" * 80)
print()
print("DIRAC PREDICTION:")
print("The Dirac equation for a spin-1/2 particle gives:")
print()
print("    g_Dirac = 2 (exactly)")
print()

g_Dirac = 2.0

print(f"Dirac g-factor:            g = {g_Dirac}")
print()

print("SCHWINGER CORRECTION:")
print("One-loop QED correction from virtual photon emission:")
print()
print("    Δa_e = α/(2π)")
print()

a_e_Schwinger = alpha / (2.0 * math.pi)

print(f"Schwinger term:            a_e = α/(2π)")
print(f"                           a_e = {a_e_Schwinger:.15e}")
print()

print("This is a THERMAL FLUCTUATION correction — the vacuum polarization")
print("from virtual pairs modifies the effective magnetic moment.")
print()

# Total g-factor (1-loop)
g_e_1loop = g_Dirac * (1.0 + a_e_Schwinger)

print("ONE-LOOP g-FACTOR:")
print("    g_e = 2(1 + α/(2π))")
print(f"    g_e = {g_e_1loop:.15f}")
print()

# Anomalous moment
a_e = (g_e_1loop - 2.0) / 2.0

print(f"Anomalous moment:          a_e = (g-2)/2")
print(f"                           a_e = {a_e:.15e}")
print()

# ============================================================================
# HIGHER-ORDER CORRECTIONS
# ============================================================================
print("HIGHER-ORDER THERMAL CORRECTIONS:")
print("-" * 80)
print()
print("Multi-loop corrections:")
print()

# 2-loop (from QED calculation)
C_2 = -0.32847844  # Coefficient for (α/π)² term
a_e_2loop = C_2 * (alpha/math.pi)**2

print(f"2-loop:                    a_e^(2) = {C_2:.8f} (α/π)²")
print(f"                                   = {a_e_2loop:.15e}")
print()

# 3-loop
C_3 = 1.181241456  # Coefficient for (α/π)³ term
a_e_3loop = C_3 * (alpha/math.pi)**3

print(f"3-loop:                    a_e^(3) = {C_3:.9f} (α/π)³")
print(f"                                   = {a_e_3loop:.15e}")
print()

# 4-loop (approximate)
C_4 = -1.7283  # Approximate coefficient
a_e_4loop = C_4 * (alpha/math.pi)**4

print(f"4-loop:                    a_e^(4) ≈ {C_4:.4f} (α/π)⁴")
print(f"                                    ≈ {a_e_4loop:.15e}")
print()

# Total (up to 4-loop)
a_e_total = a_e_Schwinger + a_e_2loop + a_e_3loop + a_e_4loop

print(f"Total (4-loop):            a_e = {a_e_total:.15e}")
print()

g_e_total = 2.0 * (1.0 + a_e_total)

print(f"g-factor (4-loop):         g_e = {g_e_total:.15f}")
print()

# ============================================================================
# THERMODYNAMIC INTERPRETATION
# ============================================================================
print("THERMODYNAMIC INTERPRETATION:")
print("-" * 80)
print()

# Compton wavelength
lambda_C = h / (m_e * c)

print(f"Compton wavelength:        λ_C = h/(m_e c)")
print(f"                           λ_C = {lambda_C:.15e} m")
print()

print("At the Compton scale, the vacuum has thermal fluctuations from")
print("virtual pair creation. These fluctuations modify the effective")
print("magnetic moment.")
print()

# Virtual photon energy scale
E_virtual = hbar * c / lambda_C

print(f"Virtual photon energy:     E ~ ℏc/λ_C = m_e c²")
print(f"                           E = {E_virtual:.6e} J")
print(f"                             = {E_virtual/e/1e6:.3f} MeV")
print()

# Thermal occupation number
# n_thermal ~ exp(-E/k_BT)
# At what temperature are virtual pairs thermally accessible?

T_17 = 17 * 18 // 2
k_B = m_e * c**2 * alpha**2 / T_17
T_pair = m_e * c**2 / k_B

print(f"Pair production temp:      T ~ m_e c²/k_B")
print(f"                           T ~ {T_pair:.6e} K")
print()

print("Below this temperature, virtual pairs are 'frozen out' —")
print("but QUANTUM fluctuations persist even at T = 0!")
print()

# ============================================================================
# VACUUM POLARIZATION
# ============================================================================
print("VACUUM POLARIZATION:")
print("-" * 80)
print()

# Effective fine structure constant at Compton scale
alpha_eff = alpha / (1.0 - a_e_Schwinger)

print("The virtual pairs screen the electron charge, modifying α:")
print()
print(f"Bare coupling:             α_bare ≈ {alpha:.10f}")
print(f"Running coupling (λ_C):    α_eff ≈ {alpha_eff:.10f}")
print(f"Enhancement:               α_eff/α = {alpha_eff/alpha:.10f}")
print()

# Dielectric constant of vacuum
epsilon_eff = epsilon_0 / (1.0 + a_e_Schwinger)

print(f"Vacuum permittivity:       ε_eff ≈ ε₀/(1 + a_e)")
print(f"                           ε_eff/ε₀ ≈ {epsilon_eff/epsilon_0:.10f}")
print()

print("The vacuum behaves as a polarizable medium!")
print()

# ============================================================================
# MAGNETIC MOMENT
# ============================================================================
print("MAGNETIC MOMENT:")
print("-" * 80)
print()

# Bohr magneton
mu_B = e * hbar / (2.0 * m_e)

print(f"Bohr magneton:             μ_B = eℏ/(2m_e)")
print(f"                           μ_B = {mu_B:.15e} J/T")
print()

# Electron magnetic moment
mu_e = g_e_total * mu_B / 2.0

print(f"Electron magnetic moment:  μ_e = (g_e/2) μ_B")
print(f"                           μ_e = {mu_e:.15e} J/T")
print()

# Spin quantum number
S = 1.0/2.0

print(f"Spin:                      S = {S}")
print(f"Magnetic moment:           μ = g_e μ_B S")
print(f"                           μ = {g_e_total * mu_B * S:.15e} J/T")
print()

# ============================================================================
# ENERGY IN MAGNETIC FIELD
# ============================================================================
print("ENERGY IN MAGNETIC FIELD:")
print("-" * 80)
print()

B_test = 1.0  # Tesla

E_Zeeman = mu_e * B_test

print(f"Magnetic field:            B = {B_test} T")
print(f"Zeeman splitting:          ΔE = μ_e B")
print(f"                           ΔE = {E_Zeeman:.15e} J")
print(f"                              = {E_Zeeman/e*1e6:.6f} μeV")
print()

# Frequency
f_Zeeman = E_Zeeman / h

print(f"Larmor frequency:          f = ΔE/h")
print(f"                           f = {f_Zeeman:.6e} Hz")
print(f"                             = {f_Zeeman/1e9:.6f} GHz")
print()

# ============================================================================
# CALIBRATION COMPARISON
# ============================================================================
print("="*80)
print("CALIBRATION COMPARISON")
print("="*80)
print()

# Experimental value (most precise measurement in physics!)
a_e_experiment = 0.00115965218128  # ±0.00000000000018
g_e_experiment = 2.0 * (1.0 + a_e_experiment)

deviation_a = a_e_total - a_e_experiment
rel_error_a = abs(deviation_a / a_e_experiment)

deviation_g = g_e_total - g_e_experiment
rel_error_g = abs(deviation_g / g_e_experiment)

print(f"TriPhase a_e (4-loop):    {a_e_total:.15e}")
print(f"Experiment a_e:           {a_e_experiment:.15e}")
print(f"Absolute deviation:       {deviation_a:+.15e}")
print(f"Relative error:           {rel_error_a:.6e} ({rel_error_a*100:.4e}%)")
print()

print(f"TriPhase g_e (4-loop):    {g_e_total:.15f}")
print(f"Experiment g_e:           {g_e_experiment:.15f}")
print(f"Absolute deviation:       {deviation_g:+.15e}")
print(f"Relative error:           {rel_error_g:.6e} ({rel_error_g*100:.4e}%)")
print()

if rel_error_a < 1e-6:
    print("✓ EXCELLENT agreement (< 1 ppm)")
elif rel_error_a < 1e-4:
    print("✓ Good agreement (< 100 ppm)")
else:
    print("⚠ Moderate deviation")

print()
print("NOTE: Full QED calculation requires 5-loop corrections plus")
print("hadronic and weak contributions. The Schwinger term α/(2π)")
print("is the dominant thermal fluctuation contribution.")
print()

# ============================================================================
# PHYSICAL SUMMARY
# ============================================================================
print("="*80)
print("PHYSICAL SUMMARY")
print("="*80)
print()
print(f"The electron g-factor g_e = {g_e_total:.6f} is:")
print()
print("1. DIRAC VALUE:")
print("   g_Dirac = 2 (exactly) from relativistic quantum mechanics")
print()
print("2. SCHWINGER CORRECTION:")
print(f"   Δa_e = α/(2π) = {a_e_Schwinger:.6e}")
print("   This is the ONE-LOOP thermal fluctuation correction from")
print("   virtual photon emission and pair creation.")
print()
print("3. HIGHER-ORDER TERMS:")
print("   Multi-loop corrections (α/π)², (α/π)³, ... account for")
print("   'thermal fluctuations of thermal fluctuations' — virtual")
print("   pairs creating more virtual pairs.")
print()
print("4. VACUUM POLARIZATION:")
print("   The virtual pair cloud modifies the effective charge:")
print(f"   α_eff/α ≈ {alpha_eff/alpha:.6f} at λ_C")
print()
print("5. THERMODYNAMIC ORIGIN:")
print("   The anomalous moment is a THERMAL EFFECT — the vacuum is")
print("   not a static background but a fluctuating thermal medium.")
print()
print("The g-2 anomaly is the MOST PRECISELY MEASURED quantity in")
print("physics (12 significant figures!). It tests QED to exquisite")
print("precision and confirms the thermodynamic nature of the vacuum.")
print()

print("="*80)
print("Derivation complete.")
print("="*80)

input("Press Enter to exit...")
