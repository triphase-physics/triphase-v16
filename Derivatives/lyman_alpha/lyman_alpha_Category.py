"""
TriPhase V16: Lyman Alpha Wavelength - Category Theory Framework
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

Category Theory Interpretation:
The Lyman alpha wavelength λ_Ly = h/(m_e·c·α) is a morphism in the category of
atomic spectral lines. It represents the colimit of the composition: electron
mass → momentum → wavelength, scaled by the fine structure constant α. The
functor F: ParticleMass → SpectralLine maps from the category of particle
properties to observable photon wavelengths. The factor α appears as a natural
transformation accounting for the first excited state of hydrogen. This reveals
atomic spectra as emergent from the adjunction between particle mass and
electromagnetic coupling, with λ_Ly as the initial object in the hydrogen
spectral series.

TAG: (D)
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
print("CATEGORY THEORY: Lyman Alpha Wavelength")
print("=" * 70)
print()

print("CATEGORICAL STRUCTURE:")
print("  Object P: Particle properties (m_e)")
print("  Object S: Spectral lines (λ_Ly)")
print("  Morphism λ: P → S (mass → wavelength)")
print("  Functor F: ParticleMass → AtomicSpectrum")
print("  Natural transformation α: ground state → excited state")
print()

print("COMMUTATIVE DIAGRAM:")
print("       m_e ──────×c──────→ p = m_e·c (momentum)")
print("        │                     │")
print("        │ ×α                  │ ×α (excitation)")
print("        ↓                     ↓")
print("   m_e·c·α ────h/x────→ λ_Ly = h/(m_e·c·α)")
print()
print("  de Broglie relation: λ = h/p, scaled by α for n=2 → n=1")
print()

print("DERIVATION:")
print(f"  Planck constant:      h   = {h:.12e} J·s")
print(f"  Electron mass:        m_e = {m_e:.12e} kg")
print(f"  Speed of light:       c   = {c:.6e} m/s")
print(f"  Fine structure:       α   = {alpha:.10f}")
print()
print(f"  Characteristic momentum scale:")
print(f"    m_e × c               = {m_e * c:.12e} kg·m/s")
print(f"    m_e × c × α           = {m_e * c * alpha:.12e} kg·m/s")
print()

lambda_Ly = h / (m_e * c * alpha)

print(f"  λ_Ly = h/(m_e·c·α)    = {lambda_Ly:.12e} m")
print(f"  λ_Ly                  = {lambda_Ly * 1e9:.6f} nm")
print()

# ========== CALIBRATION CHECKPOINT ==========
lambda_Ly_NIST = 1.21567e-7  # m (121.567 nm)
error_ppm = abs(lambda_Ly - lambda_Ly_NIST) / lambda_Ly_NIST * 1e6

print("CALIBRATION:")
print(f"  NIST (Lyman α)        = {lambda_Ly_NIST:.12e} m")
print(f"  NIST (Lyman α)        = {lambda_Ly_NIST * 1e9:.6f} nm")
print(f"  Error                 = {error_ppm:.3f} ppm")
print()

# Show energy of this transition
E_Ly = h * c / lambda_Ly
E_Ly_eV = E_Ly / 1.602176634e-19

print(f"  Transition energy:")
print(f"    E = hc/λ              = {E_Ly:.12e} J")
print(f"    E                     = {E_Ly_eV:.6f} eV")
print("    (Hydrogen n=2 → n=1 transition: 10.2 eV)")
print()

# ========== CATEGORY THEORY INSIGHT ==========
print("CATEGORY THEORY INSIGHT:")
print("Lyman alpha is the initial object in the category of hydrogen spectral")
print("lines, corresponding to the n=2 → n=1 transition. The morphism")
print("λ: m_e → λ_Ly is a functor from particle mass to photon wavelength,")
print("mediated by the natural transformation α (fine structure constant).")
print("This reveals atomic spectra as emergent from the adjunction between")
print("matter (m_e) and radiation (λ). The factor α appears because the")
print("transition involves a change in electromagnetic potential - it's the")
print("coupling constant between atomic states. The de Broglie relation")
print("λ = h/p is a special case of this functor for free particles (α → 1).")
print("The categorical perspective: all hydrogen transitions factor through")
print("λ_Ly via morphisms scaled by n² differences (Rydberg formula). This")
print("demonstrates spectroscopy as a branch of category theory, with each")
print("spectral line a morphism in the observable photon category.")
print("=" * 70)

input("Press Enter to exit...")
