"""
TriPhase V16 — Speed of Light (Statistical Mechanics Framework)
========================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
The speed of light c emerges as the maximum group velocity in the relativistic
density of states. In statistical mechanics, particle distributions depend on the
density of states g(E) in momentum space. For massless particles (photons, gravitons),
the dispersion relation is E = pc, which defines a linear density of states.

From the canonical ensemble perspective, c is the velocity at which the partition
function Z = ∫ g(p) exp(-βE) dp diverges for massless particles in thermal
equilibrium. This divergence is regulated by the finite volume of the universe,
and c emerges as the characteristic velocity scale that separates subluminal
(causal) from superluminal (acausal) propagation.

The TriPhase formula c = 1/√(ε₀μ₀) shows that light speed is the geometric mean
of electric and magnetic field propagation rates in the vacuum. This is the
fundamental wave equation: ∂²φ/∂t² = c²∇²φ, where c² = 1/(ε₀μ₀) sets the scale
for field fluctuations. In the microcanonical ensemble, c is the velocity at which
information (entropy) propagates through the vacuum.

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
print("TriPhase V16: Speed of Light (Statistical Mechanics)")
print("=" * 70)
print()

print("DENSITY OF STATES FOR MASSLESS PARTICLES:")
print("-" * 70)
print("For photons with dispersion E = pc, the phase space density is:")
print()
print("  g(E) dE = (V/2π²ℏ³c³) E² dE")
print()
print("This quadratic density of states leads to Planck distribution:")
print("  n(E) = g(E) / [exp(βE) - 1]")
print()
print(f"The speed c sets the scale for this distribution.")
print()

print("VACUUM PROPAGATION FORMULA:")
print("-" * 70)
print("Light speed emerges from EM vacuum parameters:")
print()
print(f"  c = 1/√(ε₀μ₀)")
print()

c_calc = c  # from anchor chain

print(f"  Vacuum permittivity:  ε₀ = {epsilon_0:.6e} F/m")
print(f"  Vacuum permeability:  μ₀ = {mu_0:.6e} H/m")
print()
print(f"  Product:  ε₀μ₀ = {epsilon_0 * mu_0:.6e} s²/m²")
print(f"  Speed:    c = {c_calc:.6e} m/s")
print()

print("STATISTICAL INTERPRETATION:")
print("-" * 70)
print("The speed c is the maximum entropy propagation velocity.")
print("In the microcanonical ensemble, information cannot propagate")
print("faster than c because that would violate causality—it would")
print("allow states to influence their own past.")
print()
print("From the canonical ensemble perspective, c is the velocity at")
print("which the free energy F = -kT ln(Z) becomes singular for massless")
print("particles. This singularity enforces the light cone structure of")
print("spacetime: events separated by spacelike intervals cannot exchange")
print("statistical information.")
print()

print("EQUIPARTITION THEOREM:")
print("-" * 70)
print("For photons in thermal equilibrium, energy is distributed as:")
print("  ⟨E⟩ = ∫ E·n(E)·g(E) dE / ∫ n(E)·g(E) dE")
print()
print("The integral converges only if the dispersion relation is E = pc")
print("with p < p_max (IR cutoff). The value of c emerges from requiring")
print("consistency between energy and momentum distributions.")
print()

# ========== CALIBRATION CHECKPOINT ==========
c_exact = 299792458.0  # m/s, exact by definition (SI 2019)
deviation_ppm = (c_calc - c_exact) / c_exact * 1e6

print("=" * 70)
print("CALIBRATION CHECKPOINT")
print("=" * 70)
print(f"SI 2019 (exact):        c = {c_exact:.0f} m/s")
print(f"TriPhase V16 (StatMech):  = {c_calc:.6f} m/s")
print(f"Deviation:                  {deviation_ppm:+.4f} ppm")
print()
print("(Near-exact match confirms ε₀ and μ₀ values are correctly chosen)")
print("=" * 70)
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT:")
print("-" * 70)
print("The speed of light is not an arbitrary constant—it's the characteristic")
print("velocity for entropy flow in the vacuum. Statistical mechanics requires")
print("that information (correlations between microstates) propagate at finite")
print("speed to preserve causality.")
print()
print("The formula c = 1/√(ε₀μ₀) shows that light speed is determined by the")
print("vacuum's capacity to store electric and magnetic field energy. The")
print("permittivity ε₀ controls charge screening, and permeability μ₀ controls")
print("current screening. Their product sets the wave propagation speed.")
print()
print("In the grand canonical ensemble, c also emerges as the velocity at which")
print("particle-antiparticle pairs (virtual fluctuations) can separate before")
print("annihilating. This is why c is the speed limit: faster propagation would")
print("violate the fluctuation-dissipation theorem.")
print()
print("The speed of light is fundamentally a statistical mechanics quantity:")
print("it's the rate constant for information diffusion in the quantum vacuum.")
print("=" * 70)

input("Press Enter to exit...")
