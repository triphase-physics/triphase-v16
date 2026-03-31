"""
TriPhase V16 - Lyman Alpha Wavelength (QFT Framework)
======================================================

QFT INTERPRETATION:
The Lyman-alpha transition (n=2→1 in hydrogen) is the first resonance in the
bound state spectrum of QED:
- Arises from solving the Dirac equation with Coulomb potential
- Fine structure splitting from spin-orbit coupling (QED radiative corrections)
- Lamb shift from vacuum polarization modifies energy levels
- Bethe's calculation of the Lamb shift was an early triumph of renormalized QED

TriPhase's formula λ_Lyα = h/(m_e c α) connects the wavelength to the electron
Compton wavelength scaled by α⁻¹, revealing the Lyman-alpha photon energy as
the electromagnetic binding energy at the Bohr radius scale.

(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
TAG: (D) - Direct derivation from Rydberg formula
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

# ========== QFT DERIVATION: LYMAN ALPHA WAVELENGTH ==========
print("=" * 70)
print("TriPhase V16 - Lyman Alpha Wavelength")
print("QFT Framework: Bound State Spectrum in QED")
print("=" * 70)
print()

print("QFT CONTEXT:")
print("The hydrogen spectrum arises from solving QED for a bound electron-proton")
print("system. The Lyman-alpha line (n=2→1 transition) at ~121.6 nm is the most")
print("prominent feature in the UV spectrum. QED corrections include:")
print("- Lamb shift: vacuum polarization shifts 2S₁/₂ level by ~1 GHz")
print("- Anomalous magnetic moment: modifies hyperfine structure")
print("- Recoil corrections: finite proton mass effects")
print()

print("TRIPHASE DERIVATION:")
print("λ_Lyα = h / (m_e × c × α)")
print()
print(f"Planck constant:      h = {h:.10e} J·s")
print(f"Electron mass:        m_e = {m_e:.10e} kg")
print(f"Speed of light:       c = {c:.10e} m/s")
print(f"Fine structure:       α = {alpha:.12f}")
print(f"m_e × c × α =         {m_e * c * alpha:.10e}")
print(f"λ_Lyα (TriPhase):     {h / (m_e * c * alpha):.10e} m")
print(f"                      {h / (m_e * c * alpha) * 1e9:.6f} nm")
print()

# Relation to Compton wavelength
lambda_C = h / (m_e * c)
lambda_Ly = h / (m_e * c * alpha)
print(f"Compton wavelength:   λ_C = h/(m_e c) = {lambda_C:.10e} m")
print(f"Lyman-alpha:          λ_Lyα = λ_C / α = {lambda_Ly:.10e} m")
print(f"Ratio λ_Lyα / λ_C:    {lambda_Ly / lambda_C:.6f} ≈ α⁻¹")
print()

# ========== CALIBRATION CHECKPOINT ==========
codata_lambda_Ly = 1.21567e-7  # m (approximate)
deviation_ppm = (lambda_Ly - codata_lambda_Ly) / codata_lambda_Ly * 1e6

print("CALIBRATION CHECKPOINT:")
print(f"Observed:             {codata_lambda_Ly:.10e} m ({codata_lambda_Ly*1e9:.4f} nm)")
print(f"TriPhase:             {lambda_Ly:.10e} m ({lambda_Ly*1e9:.4f} nm)")
print(f"Deviation:            {deviation_ppm:+.2f} ppm")
print()

# ========== QFT INSIGHT ==========
print("QFT INSIGHT:")
print("The factor α⁻¹ ≈ 137 amplification from Compton to Lyman-alpha wavelength")
print("reflects the electromagnetic binding energy structure. In QFT terms:")
print("- The Compton wavelength sets the scale for virtual photon exchange")
print("- The Bohr radius a₀ ~ λ_C/α is where the Coulomb potential balances kinetic energy")
print("- The Rydberg energy Ry ~ m_e c² α² / 2 sets the binding energy scale")
print()
print("The Lyman-alpha photon energy E_Lyα ~ Ry represents the difference between")
print("n=2 and n=1 bound states, a pure QED effect from non-perturbative solving")
print("of the Dirac equation in an external Coulomb field.")
print()
print("=" * 70)

input("Press Enter to exit...")
