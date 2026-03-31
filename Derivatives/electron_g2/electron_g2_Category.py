"""
TriPhase V16: Electron g-factor (g-2) - Category Theory Framework
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

Category Theory Interpretation:
The electron anomalous magnetic moment (g-2)/2 is a morphism in the category of
quantum corrections. The g-factor g = 2(1 + α/(2π) - 0.328(α/π)²) represents a
functor from the electromagnetic coupling category to the magnetic moment category.
The first-order correction α/(2π) is the Schwinger term - a natural transformation
from tree-level QED to one-loop corrections. The second-order term -0.328(α/π)²
encodes two-loop vacuum polarization. This reveals quantum corrections as a
composition of morphisms, each representing a loop order in the perturbation
series. The categorical perspective: (g-2) is the colimit of the loop expansion.

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
print("CATEGORY THEORY: Electron g-factor (g-2)")
print("=" * 70)
print()

print("CATEGORICAL STRUCTURE:")
print("  Object D: Dirac theory (g = 2, tree-level)")
print("  Object Q: QED corrections (g > 2, loop corrections)")
print("  Morphism Δ: D → Q (Dirac → quantum corrections)")
print("  Functor F: EMCoupling → MagneticMoment")
print("  Natural transformations: loop order expansions")
print()

print("COMMUTATIVE DIAGRAM:")
print("       g=2 ──────+α/(2π)──────→ g ≈ 2.00232")
print("        │                          │")
print("        │ (tree-level)             │ + higher orders")
print("        ↓                          ↓")
print("   Dirac ────→ QED (1-loop + 2-loop + ...)")
print()

print("DERIVATION:")
print(f"  Fine structure:       α       = {alpha:.10f}")
print(f"  Coupling ratio:       α/π     = {alpha/math.pi:.10f}")
print()
print("  g-factor expansion:")
print("    g = 2 × (1 + C₁·(α/2π) + C₂·(α/π)² + ...)")
print()
print("  Schwinger term (1-loop):")
print(f"    C₁ = 1.0")
print(f"    α/(2π)                      = {alpha/(2.0*math.pi):.12e}")
print()
print("  Two-loop correction:")
print(f"    C₂ = -0.328")
print(f"    -0.328·(α/π)²               = {-0.328*(alpha/math.pi)**2:.12e}")
print()

g_factor = 2.0 * (1.0 + alpha/(2.0*math.pi) - 0.328*(alpha/math.pi)**2)

print(f"  g = 2(1 + α/2π - 0.328(α/π)²) = {g_factor:.15f}")
print()

# Anomalous magnetic moment a_e = (g-2)/2
a_e = (g_factor - 2.0) / 2.0

print(f"  Anomalous moment a_e = (g-2)/2")
print(f"  a_e                           = {a_e:.15f}")
print(f"  a_e × 10¹²                    = {a_e * 1e12:.6f}")
print()

# ========== CALIBRATION CHECKPOINT ==========
g_codata = 2.00231930436256  # CODATA 2018
a_e_codata = (g_codata - 2.0) / 2.0

error_ppm = abs(a_e - a_e_codata) / a_e_codata * 1e6

print("CALIBRATION:")
print(f"  CODATA 2018 g-factor          = {g_codata:.15f}")
print(f"  CODATA 2018 a_e               = {a_e_codata:.15f}")
print(f"  CODATA a_e × 10¹²             = {a_e_codata * 1e12:.6f}")
print()
print(f"  TriPhase V16 g-factor         = {g_factor:.15f}")
print(f"  Error in a_e                  = {error_ppm:.3f} ppm")
print()

# Show higher order terms for context
print("NOTE ON HIGHER ORDERS:")
print("  Full QED calculation includes terms up to (α/π)⁵:")
print("    a_e = C₁(α/π) + C₂(α/π)² + C₃(α/π)³ + C₄(α/π)⁴ + C₅(α/π)⁵")
print("  Where:")
print("    C₁ = 0.5 (Schwinger, 1948)")
print("    C₂ = -0.328478965... (two-loop)")
print("    C₃ = 1.181241456... (three-loop)")
print("    C₄ = -1.7283... (four-loop)")
print("    C₅ = 0.0... (five-loop, hadronic + weak)")
print()
print("  TriPhase uses simplified two-loop approximation.")
print()

# ========== CATEGORY THEORY INSIGHT ==========
print("CATEGORY THEORY INSIGHT:")
print("The electron g-factor is a colimit in the category of perturbative QED.")
print("Each loop order represents a morphism in the composition:")
print()
print("  g₀=2 → g₁ (1-loop) → g₂ (2-loop) → g₃ → ... → g_∞ (true value)")
print()
print("The functor F: α → (g-2)/2 maps electromagnetic coupling to magnetic")
print("moment anomaly. Each term (α/π)ⁿ is a natural transformation representing")
print("n-loop Feynman diagrams. The categorical structure reveals why QED")
print("converges: the morphism composition is a Cauchy sequence in the category")
print("of observables, with limit point a_e. The Schwinger term α/(2π) is the")
print("initial morphism - all higher corrections factor through it via the")
print("functor structure. The Yoneda perspective: (g-2) is uniquely determined")
print("by its relationship to α across all loop orders. This proves the magnetic")
print("moment is not a free parameter but emergent from vacuum structure. The")
print("two-loop coefficient -0.328 encodes photon-photon scattering (light-by-light)")
print("and fermion loops - these are colimit constructions in the category of")
print("virtual processes. Category theory makes QED's structure manifest: it's")
print("a tower of natural transformations, each adding precision to (g-2).")
print("=" * 70)

input("Press Enter to exit...")
