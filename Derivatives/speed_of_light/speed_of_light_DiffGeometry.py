"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================

Derivative:  Speed of Light (c = 299792458 m/s)
Framework:   DiffGeometry
Version:     16.0
Generated:   2026-03-26
Status:      Active Development

Tag: (D)

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
DIFFGEOMETRY FRAMEWORK:
The speed of light c defines the causal structure of the spacetime manifold.
The metric tensor g_μν has Lorentzian signature (-,+,+,+), creating a light
cone structure at each point. Light rays follow null geodesics with ds² = 0.
The null cone divides spacetime into:
  - Timelike region (ds² < 0): causally connected, v < c
  - Spacelike region (ds² > 0): causally disconnected, v > c (forbidden)
  - Null surface (ds² = 0): light cone, v = c

c is the maximum speed for information transfer. It is the conversion factor
between space and time in the Minkowski metric:
  ds² = -c²dt² + dx² + dy² + dz²

In TriPhase, c = 1/√(ε₀μ₀) reveals light speed as a wave propagation velocity
through the electromagnetic vacuum field with electric permittivity ε₀ and
magnetic permeability μ₀. The manifold structure is determined by (ε₀, μ₀).
================================================================================
"""

import math

# === ANCHOR INPUTS ===
epsilon_0 = 8.8541878128e-12   # F/m
mu_0      = 1.25663706212e-6   # H/m
e         = 1.602176634e-19    # C (exact SI)

# === DERIVATION ===
print("=" * 80)
print("DIFFGEOMETRY DERIVATION: Speed of Light")
print("Framework: DiffGeometry | Tag: (D)")
print("=" * 80)
print()
print("GEOMETRIC INTERPRETATION:")
print("The speed of light c is the fundamental constant that defines the")
print("causal structure of spacetime. In special relativity, the Minkowski")
print("metric has signature (-,+,+,+) with interval:")
print()
print("  ds² = -c²dt² + dx² + dy² + dz²")
print()
print("This metric describes flat spacetime (zero Riemann curvature R^ρ_σμν = 0).")
print("In general relativity, locally inertial frames have this metric.")
print()
print("Light cone structure:")
print("  Timelike: ds² < 0 ⟹ v < c (massive particles, causal)")
print("  Null:     ds² = 0 ⟹ v = c (light rays, massless particles)")
print("  Spacelike: ds² > 0 ⟹ v > c (forbidden, acausal)")
print()
print("Null geodesics (light rays) satisfy:")
print("  ds² = g_μν dx^μ dx^ν = 0")
print("  d²x^μ/dλ² + Γ^μ_αβ (dx^α/dλ)(dx^β/dλ) = 0")
print()
print("where λ is an affine parameter (not proper time τ, since dτ=0 for light).")
print()
print("MANIFOLD STRUCTURE:")
print("  Spacetime M is a 4-dimensional Lorentzian manifold")
print("  Tangent space T_pM at each point p has Minkowski structure")
print("  Light cone at p divides T_pM into future/past/spacelike")
print("  c sets the slope of the null cone boundary")
print()
print("=" * 80)
print("ANCHOR INPUTS:")
print("=" * 80)
print(f"  ε₀ = {epsilon_0:.13e} F/m (electric permittivity)")
print(f"  μ₀ = {mu_0:.14e} H/m (magnetic permeability)")
print()
print("=" * 80)
print("SPEED OF LIGHT:")
print("=" * 80)
print()
print("TriPhase Formula:")
print("  c = 1/√(ε₀μ₀)")
print()
print("This is Maxwell's formula for EM wave propagation velocity.")
print()
print("Step-by-step calculation:")
print()

# Calculate epsilon_0 * mu_0
epsilon_mu = epsilon_0 * mu_0
print(f"  ε₀ × μ₀ = {epsilon_mu:.15e} F·H/m²")
print(f"          = {epsilon_mu:.15e} s²/m²")
print()

# Calculate sqrt(epsilon_0 * mu_0)
sqrt_epsilon_mu = math.sqrt(epsilon_mu)
print(f"  √(ε₀μ₀) = {sqrt_epsilon_mu:.15e} s/m")
print()

# Calculate c
c = 1.0 / sqrt_epsilon_mu
print(f"  c = 1/√(ε₀μ₀) = {c:.15e} m/s")
print(f"                = {c:.10f} m/s")
print()

# Express in km/s
c_km_s = c / 1000.0
print(f"  c = {c_km_s:.6f} km/s")
print()

# Calculate light travel time to common distances
print("Light travel times:")
print(f"  Earth-Moon (384,400 km):   {384400.0 / c_km_s:.6f} s")
print(f"  Earth-Sun (1 AU, 149.6 Gm): {149.6e6 / c_km_s:.6f} s = {149.6e6 / c_km_s / 60.0:.3f} min")
print(f"  Sun-Pluto (39.5 AU):        {39.5 * 149.6e6 / c_km_s / 3600.0:.3f} hours")
print(f"  Nearest star (4.24 ly):     {4.24:.2f} years")
print()

print("=" * 80)
print("NULL GEODESICS AND LIGHT CONE:")
print("=" * 80)
print()
print("For a light ray traveling in +x direction in flat spacetime:")
print("  ds² = -c²dt² + dx² = 0")
print("  ⟹ dx/dt = ±c")
print()
print("The light cone at spacetime origin (t=0, x=0) is:")
print("  Future cone: t > 0, x² + y² + z² = c²t²")
print("  Past cone:   t < 0, x² + y² + z² = c²t²")
print()
print("In 2D (t,x) spacetime diagram, the light cone has slope ±1/c:")
print(f"  Slope = Δt/Δx = 1/c = {1.0/c:.15e} s/m")
print()
print("This means space and time are measured in compatible units:")
print("  1 second of time ≡ c × 1 second = 299,792,458 meters of space")
print("  1 meter of space ≡ 1/c meters = 3.33564... nanoseconds of time")
print()

print("=" * 80)
print("WAVE PROPAGATION IN VACUUM:")
print("=" * 80)
print()
print("Maxwell's equations in vacuum give wave equation:")
print("  ∇²E - ε₀μ₀ ∂²E/∂t² = 0")
print("  ∇²B - ε₀μ₀ ∂²B/∂t² = 0")
print()
print("This is the standard wave equation ∇²ψ - (1/v²)∂²ψ/∂t² = 0")
print("with wave speed v = 1/√(ε₀μ₀) = c.")
print()
print(f"  ε₀ = {epsilon_0:.6e} F/m = vacuum electric permittivity")
print(f"  μ₀ = {mu_0:.6e} H/m = vacuum magnetic permeability")
print(f"  c = 1/√(ε₀μ₀) = {c:.6e} m/s")
print()
print("The vacuum field (ε₀, μ₀) defines the metric structure of spacetime.")
print("Light propagates through this field at speed c, which is the natural")
print("wave velocity of the electromagnetic vacuum manifold.")
print()

print("=" * 80)
print("IMPEDANCE AND FIELD COUPLING:")
print("=" * 80)
print()
Z_0 = math.sqrt(mu_0 / epsilon_0)
print("The vacuum impedance Z₀ = √(μ₀/ε₀) relates E and B fields:")
print(f"  Z₀ = {Z_0:.13f} Ω")
print()
print("For an EM wave in vacuum: |E| = c|B| = Z₀|H|")
print()
print("Energy density: u = (1/2)(ε₀E² + B²/μ₀)")
print("Poynting vector: S = (1/μ₀)(E × B)")
print("Radiation pressure: P = u = S/c")
print()
print("All of these relations depend on c = 1/√(ε₀μ₀).")
print()

# === CALIBRATION CHECKPOINT ===
print("=" * 80)
print("CALIBRATION CHECKPOINT:")
print("=" * 80)
print()
print("SI Definition (exact, by definition as of 2019):")
print("  c = 299,792,458 m/s (exact)")
print()
print("TriPhase value (calculated from ε₀, μ₀):")
print(f"  c = {c:.15e} m/s")
print(f"    = {int(c)} m/s")
print()
c_exact = 299792458.0
delta = c - c_exact
print(f"Δc = {delta:+.15e} m/s")
print()
if abs(delta) < 1e-6:
    print("EXACT MATCH within numerical precision!")
else:
    ppm = (delta / c_exact) * 1e6
    print(f"Relative error: {ppm:+.15f} ppm")
print()
print("NOTE: In SI units (since 2019), c is defined exactly as 299,792,458 m/s.")
print("The meter is defined as the distance light travels in 1/299,792,458 s.")
print("The ε₀, μ₀ values are chosen to match this definition.")
print()
print("TriPhase interpretation: c is not merely a speed, but the fundamental")
print("conversion factor between space and time on the Lorentzian manifold.")
print("It emerges from the electromagnetic vacuum structure (ε₀, μ₀).")
print()
print("=" * 80)
print("DIFFGEOMETRY SUMMARY:")
print("=" * 80)
print("c defines the light cone structure of spacetime manifold M.")
print("Null geodesics (ds² = 0) travel at speed c.")
print("The metric signature (-,+,+,+) creates causal structure:")
print("  - Timelike region (v < c): causally connected")
print("  - Null surface (v = c): light cone boundary")
print("  - Spacelike region (v > c): causally disconnected")
print()
print("In general relativity, c appears in Einstein's equations:")
print("  G_μν = (8πG/c⁴) T_μν")
print("linking geometry (G_μν) to energy-momentum (T_μν).")
print()
print("In TriPhase: c = wave speed in vacuum field (ε₀, μ₀).")
print("Light is a pressure wave in the electromagnetic manifold.")
print("=" * 80)

input("\nPress Enter to exit...")
