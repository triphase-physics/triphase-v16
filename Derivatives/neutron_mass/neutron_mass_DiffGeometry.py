"""
================================================================================
TriPhase V16 Python Derivative Script
================================================================================
Derivative:  Neutron Mass (m_n = 1.67493e-27 kg = 939.565 MeV/c²)
Framework:   DiffGeometry
Version:     16.0
Generated:   2026-03-26
Status:      Active Development
Tag: (D*)
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
================================================================================
"""
import math

print("=" * 80)
print("TRIPHASE V16 - NEUTRON MASS DERIVATION")
print("Framework: DiffGeometry")
print("=" * 80)

# === ANCHOR INPUTS ===
print("\n[ANCHOR INPUTS]")
epsilon_0 = 8.8541878128e-12
mu_0      = 1.25663706212e-6
e         = 1.602176634e-19
r_e       = 2.8179403262e-15

print(f"ε₀ = {epsilon_0:.13e} F/m")
print(f"μ₀ = {mu_0:.14e} H/m")
print(f"e  = {e:.12e} C")
print(f"r_e = {r_e:.13e} m")

# === DERIVED ANCHOR CHAIN ===
print("\n[DERIVED ANCHOR CHAIN]")
c         = 1.0 / math.sqrt(epsilon_0 * mu_0)
Z_0       = math.sqrt(mu_0 / epsilon_0)
alpha_inv = 137.0 + math.log(137.0) / 137.0
alpha     = 1.0 / alpha_inv
hbar      = Z_0 * e**2 / (4.0 * math.pi * alpha)
h         = 2.0 * math.pi * hbar
G         = c**4 * 7.5 * epsilon_0**3 * mu_0**2
m_e       = hbar * alpha / (c * r_e)
f_e       = m_e * c**2 / hbar
T_17      = 17 * 18 // 2  # = 153
mp_me     = 4.0 * 27.0 * 17.0 * (1.0 + 5.0 * alpha**2 / math.pi)
m_p       = m_e * mp_me

print(f"c     = {c:.10e} m/s")
print(f"Z₀    = {Z_0:.10e} Ω")
print(f"α⁻¹   = {alpha_inv:.10f}")
print(f"α     = {alpha:.12e}")
print(f"ħ     = {hbar:.10e} J·s")
print(f"h     = {h:.10e} J·s")
print(f"G     = {G:.10e} m³/(kg·s²)")
print(f"m_e   = {m_e:.10e} kg")
print(f"f_e   = {f_e:.10e} Hz")
print(f"T₁₇   = {T_17}")
print(f"mp/me = {mp_me:.10f}")
print(f"m_p   = {m_p:.11e} kg")

# === DERIVATION ===
print("\n[NEUTRON MASS DERIVATION]")
print("\nStep 1: Neutron-proton mass splitting")
print("  m_n = m_p × (1 + α/(2π×T₁₇))")
print("  Isospin symmetry breaking correction")

correction = alpha / (2.0 * math.pi * T_17)
print(f"  α/(2π×T₁₇) = {correction:.12e}")
print(f"  Correction factor: 1 + {correction:.12e} = {1.0 + correction:.12f}")

print("\nStep 2: Neutron mass")
m_n = m_p * (1.0 + correction)

print(f"  m_n = {m_n:.11e} kg")

# Convert to MeV/c²
m_n_MeV = m_n * c**2 / (e * 1e6)
m_p_MeV = m_p * c**2 / (e * 1e6)

print(f"  m_n = {m_n_MeV:.6f} MeV/c²")

print("\nStep 3: Mass difference")
delta_m = m_n - m_p
delta_m_MeV = m_n_MeV - m_p_MeV

print(f"  Δm = m_n - m_p = {delta_m:.11e} kg")
print(f"  Δm = {delta_m_MeV:.6f} MeV/c²")

# === DIFFERENTIAL GEOMETRY INTERPRETATION ===
print("\n[DIFFERENTIAL GEOMETRY INTERPRETATION]")
print("\nNeutron-proton mass splitting as curvature perturbation:")
print("  - The isospin manifold has gauge group SU(2)_isospin")
print("  - Proton (uud) and neutron (udd) are geodesics on this manifold")
print("  - In perfect isospin symmetry, m_p = m_n (zero curvature difference)")
print("  - Electromagnetic and quark mass effects break the symmetry")
print("  - This creates a metric perturbation g_μν → g_μν + δg_μν")
print("  - The perturbation magnitude δg/g ~ α/(2π×T₁₇)")
print("  - Factor α → electromagnetic coupling (photon exchange)")
print("  - Factor 1/(2π×T₁₇) → TriPhase suppression through 153 modes")
print("\nGeometric picture:")
print("  - Proton and neutron are nearby points on the isospin manifold")
print("  - Geodesic distance Δs corresponds to mass difference Δm")
print("  - Riemann curvature tensor R^μ_νρσ encodes u-d quark distinction")
print("  - Mass splitting Δm ~ ∫ R ds along the isospin geodesic")
print("\nPhysical mechanism:")
print("  - u quark (charge +2e/3) vs d quark (charge -e/3)")
print("  - Proton (net +e) has lower electromagnetic self-energy")
print("  - Neutron (net 0) has higher quark mass contribution")
print("  - Result: m_n > m_p by ~1.29 MeV (enables beta decay)")

# === CALIBRATION CHECKPOINT ===
print("\n[CALIBRATION CHECKPOINT]")
m_n_CODATA = 1.67492749804e-27  # kg
m_n_CODATA_MeV = 939.56542052  # MeV/c²
delta_m_CODATA_MeV = 1.29333236  # MeV/c²

print(f"TriPhase m_n:  {m_n:.11e} kg = {m_n_MeV:.6f} MeV/c²")
print(f"CODATA m_n:    {m_n_CODATA:.11e} kg = {m_n_CODATA_MeV:.6f} MeV/c²")

error_kg = abs(m_n - m_n_CODATA) / m_n_CODATA * 100
error_MeV = abs(m_n_MeV - m_n_CODATA_MeV) / m_n_CODATA_MeV * 100

print(f"Error:         {error_kg:.6f}% (kg)")
print(f"               {error_MeV:.6f}% (MeV)")

print(f"\nTriPhase Δm:   {delta_m_MeV:.6f} MeV/c²")
print(f"CODATA Δm:     {delta_m_CODATA_MeV:.6f} MeV/c²")
error_delta = abs(delta_m_MeV - delta_m_CODATA_MeV) / delta_m_CODATA_MeV * 100
print(f"Error (Δm):    {error_delta:.6f}%")

if error_kg < 0.1:
    print("\nSTATUS: EXCELLENT AGREEMENT (<0.1% error)")
elif error_kg < 1.0:
    print("\nSTATUS: GOOD AGREEMENT (<1% error)")
else:
    print("\nSTATUS: CALIBRATION NEEDED (>1% error)")

print("\n" + "=" * 80)
print("DERIVATION COMPLETE")
print("=" * 80)

input("\nPress Enter to exit...")
