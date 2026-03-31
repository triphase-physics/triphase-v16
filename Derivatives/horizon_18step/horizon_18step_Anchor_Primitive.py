"""
TriPhase V16 Derivative - 18-Step Horizon Distance (ANCHOR_PRIMITIVE)
Framework: Anchor_Primitive
Tag: (D)
DOI: 10.5281/zenodo.17855383

ANCHOR PRIMITIVE DERIVATION:
- ONLY inputs: epsilon_0 = 8.8541878128e-12 F/m, mu_0 = 1.25663706212e-6 H/m
- e = 1.602176634e-19 C (exact SI definition)
- c = 1/sqrt(epsilon_0 * mu_0)
- Z_0 = sqrt(mu_0/epsilon_0)
- alpha_inv = 137 + ln(137)/137, alpha = 1/alpha_inv
- hbar = Z_0 * e^2 / (4*pi*alpha)
- m_e from hbar, alpha, c chain
- f_e = m_e * c^2 / h (electron frequency)
- H_0 = pi*sqrt(3) * f_e * alpha^18 (18-step Hubble constant)
- d_horizon = c / H_0 (Hubble radius)

The 18-step cascade: alpha^18 connects electron frequency to cosmic expansion.

PURE ANCHOR CHAIN - No shortcuts, full derivation from electromagnetic constants.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
"""

import math

print("="*70)
print("TriPhase V16 - 18-Step Horizon Distance")
print("Framework: Anchor_Primitive | Tag: (D)")
print("="*70)
print()

# ============================================================================
# ANCHOR PRIMITIVE INPUTS (ONLY)
# ============================================================================
epsilon_0 = 8.8541878128e-12  # F/m (permittivity of free space)
mu_0 = 1.25663706212e-6       # H/m (permeability of free space)
e = 1.602176634e-19           # C (elementary charge, exact SI definition)

print("ANCHOR PRIMITIVE INPUTS:")
print(f"  epsilon_0 = {epsilon_0:.13e} F/m")
print(f"  mu_0      = {mu_0:.14e} H/m")
print(f"  e         = {e:.12e} C (exact SI)")
print()

# ============================================================================
# DERIVED FUNDAMENTAL CONSTANTS
# ============================================================================
print("DERIVED CONSTANTS FROM ANCHOR CHAIN:")
print()

# Speed of light
c = 1.0 / math.sqrt(epsilon_0 * mu_0)
print(f"c = 1/sqrt(epsilon_0 * mu_0)")
print(f"  = {c:.10e} m/s")
print()

# Impedance of free space
Z_0 = math.sqrt(mu_0 / epsilon_0)
print(f"Z_0 = sqrt(mu_0/epsilon_0)")
print(f"    = {Z_0:.10f} Ohms")
print()

# Fine structure constant (TriPhase derivation)
alpha_inv = 137.0 + math.log(137.0)/137.0
alpha = 1.0 / alpha_inv
print(f"alpha_inv = 137 + ln(137)/137 = {alpha_inv:.10f}")
print(f"alpha     = 1/alpha_inv = {alpha:.12e}")
print()

# Reduced Planck constant
hbar = Z_0 * e**2 / (4.0 * math.pi * alpha)
print(f"hbar = Z_0 * e^2 / (4*pi*alpha)")
print(f"     = {hbar:.15e} J·s")
print()

# Planck constant
h = 2.0 * math.pi * hbar
print(f"h = 2*pi*hbar")
print(f"  = {h:.15e} J·s")
print()

# Compton wavelength of electron
lambda_C = 2.0 * math.pi / (alpha**2 * 137.0)
print(f"lambda_C = 2*pi / (alpha^2 * 137)")
print(f"         = {lambda_C:.15e} m")
print()

# Electron mass
m_e = hbar / (c * lambda_C)
print(f"m_e = hbar / (c * lambda_C)")
print(f"    = {m_e:.15e} kg")
print()

# Electron frequency (rest energy frequency)
f_e = m_e * c**2 / h
print(f"f_e = m_e * c^2 / h")
print(f"    = {f_e:.10e} Hz")
print()

# ============================================================================
# 18-STEP HUBBLE CONSTANT DERIVATION
# ============================================================================
print("="*70)
print("18-STEP HUBBLE CONSTANT DERIVATION:")
print("="*70)
print()

print("The 18-step cascade connects quantum to cosmic scales:")
print("  H_0 = pi*sqrt(3) * f_e * alpha^18")
print()
print("Each alpha step represents a scale transition.")
print("18 steps: electron -> nucleus -> atom -> ... -> cosmic horizon")
print()

# 18-step scaling factor
alpha_18 = alpha**18
print(f"alpha^18 = {alpha:.12e}^18")
print(f"         = {alpha_18:.15e}")
print()

# Geometric prefactor
geo_factor = math.pi * math.sqrt(3.0)
print(f"Geometric prefactor: pi*sqrt(3) = {geo_factor:.10f}")
print()

# Hubble constant (18-step)
H_0 = geo_factor * f_e * alpha_18
print(f"H_0 = pi*sqrt(3) * f_e * alpha^18")
print(f"    = {H_0:.15e} Hz")
print(f"    = {H_0:.15e} s^-1")
print()

# Convert to km/s/Mpc (standard cosmology units)
Mpc_to_m = 3.0857e22  # meters per megaparsec
H_0_km_s_Mpc = H_0 * Mpc_to_m / 1000.0
print(f"H_0 = {H_0_km_s_Mpc:.4f} km/s/Mpc")
print()

# ============================================================================
# HORIZON DISTANCE CALCULATION
# ============================================================================
print("="*70)
print("HORIZON DISTANCE (HUBBLE RADIUS):")
print("="*70)
print()

print("Hubble radius: d_horizon = c / H_0")
print("(Comoving distance to particle horizon)")
print()

# Hubble radius
d_horizon = c / H_0
print(f"d_horizon = c / H_0")
print(f"          = {d_horizon:.10e} m")
print()

# Convert to light-years
ly_to_m = 9.4607304725808e15  # meters per light-year
d_horizon_ly = d_horizon / ly_to_m
print(f"d_horizon = {d_horizon_ly:.10e} light-years")
print(f"          = {d_horizon_ly/1e9:.4f} billion light-years (Gly)")
print()

# Convert to Gpc (gigaparsecs)
Gpc_to_m = 3.0857e25  # meters per gigaparsec
d_horizon_Gpc = d_horizon / Gpc_to_m
print(f"d_horizon = {d_horizon_Gpc:.6f} Gpc")
print()

# Universe age estimate (1/H_0 for flat universe)
t_universe = 1.0 / H_0
t_universe_years = t_universe / (365.25 * 24.0 * 3600.0)
t_universe_Gyr = t_universe_years / 1e9

print(f"Universe age estimate: t = 1/H_0")
print(f"                         = {t_universe:.10e} s")
print(f"                         = {t_universe_Gyr:.4f} Gyr")
print()

# ============================================================================
# OBSERVATIONAL COMPARISON (CALIBRATION CHECKPOINT)
# ============================================================================
print("="*70)
print("OBSERVATIONAL COMPARISON (CALIBRATION CHECKPOINT):")
print("="*70)

# Planck 2018 values
H_0_Planck = 67.66  # km/s/Mpc (Planck 2018 CMB)
H_0_SH0ES = 73.04   # km/s/Mpc (SH0ES/HST 2022)
d_horizon_Planck = c / (H_0_Planck * 1000.0 / Mpc_to_m)
d_horizon_SH0ES = c / (H_0_SH0ES * 1000.0 / Mpc_to_m)

print(f"TriPhase H_0:  {H_0_km_s_Mpc:.4f} km/s/Mpc")
print(f"Planck 2018:   {H_0_Planck:.4f} km/s/Mpc (CMB)")
print(f"SH0ES 2022:    {H_0_SH0ES:.4f} km/s/Mpc (Local distance ladder)")
print()
print(f"TriPhase matches: {'Planck' if abs(H_0_km_s_Mpc - H_0_Planck) < abs(H_0_km_s_Mpc - H_0_SH0ES) else 'SH0ES'}")
print()
print(f"Horizon distances:")
print(f"  TriPhase: {d_horizon_Gpc:.6f} Gpc")
print(f"  Planck:   {d_horizon_Planck/Gpc_to_m:.6f} Gpc")
print(f"  SH0ES:    {d_horizon_SH0ES/Gpc_to_m:.6f} Gpc")
print()
print(f"Observable universe diameter ~ 2 * d_horizon:")
print(f"  TriPhase: {2.0*d_horizon_Gpc:.4f} Gpc = {2.0*d_horizon_ly/1e9:.2f} Gly")
print()

# ============================================================================
# 18-STEP CASCADE BREAKDOWN
# ============================================================================
print("="*70)
print("18-STEP CASCADE STRUCTURE:")
print("="*70)
print()
print("Each alpha step represents a characteristic scale transition:")
print()

scales = [
    (0, "Electron Compton frequency", f_e),
    (1, "Fine structure", f_e * alpha),
    (2, "Atomic scale", f_e * alpha**2),
    (3, "Molecular scale", f_e * alpha**3),
    (6, "Chemical/thermal", f_e * alpha**6),
    (9, "Macroscopic", f_e * alpha**9),
    (12, "Planetary", f_e * alpha**12),
    (15, "Stellar", f_e * alpha**15),
    (18, "Cosmic horizon (H_0)", f_e * alpha**18),
]

for step, name, freq in scales:
    wavelength = c / freq
    print(f"  Step {step:2d}: {name:25s} f = {freq:.4e} Hz, λ = {wavelength:.4e} m")

print()
print("="*70)
print("Derivation complete. Quantum-to-cosmic bridge via 18 alpha steps.")
print("Strong foundations: Horizon distance traced to epsilon_0, mu_0.")
print("="*70)

input("Press Enter to exit...")
