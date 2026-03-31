"""
TriPhase V16 — Down Quark Mass (Symplectic Framework)
======================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

SYMPLECTIC INTERPRETATION
--------------------------
The down quark is the second-lightest quark, with electric charge -1/3 e.
In symplectic geometry, quark masses emerge from the triangular harmonic
structure of QCD phase space, scaled by charge fraction, α, and T₁₇.

Phase Space: (x, p, color) = (position, momentum, SU(3) color charge)
Hamiltonian: H = √(p²c² + m_d²c⁴) + V_QCD(color)

Like all quarks, the down quark is confined and has two mass scales:
current mass (bare) and constituent mass (in hadrons).

QUARK PROPERTIES
----------------
Down quark:
  Charge: -1/3 e
  Spin: 1/2
  Color: r, g, or b
  Generation: 1st (lightest down-type quark)

ISOSPIN SYMMETRY
----------------
Up and down quarks form an isospin doublet:
  u: I₃ = +1/2 (isospin up)
  d: I₃ = -1/2 (isospin down)

In the limit of zero quark masses, QCD has approximate SU(2) isospin symmetry.
The mass difference m_d - m_u breaks this symmetry.

COLOR PHASE SPACE
-----------------
Like the up quark, the down quark lives in extended phase space (x, p, color).
Symplectic form: ω = dp ∧ dx + ω_color

The SU(3) color structure is identical for all quarks.

TRIANGULAR HARMONIC
-------------------
The triangular number T₁₇ = 153 appears in quark mass formulas:
m_d ~ m_e × (charge magnitude) × α × T₁₇

For down quark with charge -1/3:
m_d ~ m_e × (4/3) × α × T₁₇

The factor 4/3 (not 1/3) accounts for the symplectic structure of the
charge-conjugate state. In TriPhase, both up and down use similar
scaling but with different charge-related factors.

POISSON BRACKET
---------------
{x, p} = 1 (position-momentum)
{color_i, color_j} ~ f_ijk color_k (SU(3) algebra)

LIOUVILLE'S THEOREM
-------------------
Phase space volume preserved in spacetime and color space.

TRIPHASE FORMULA
----------------
m_d ~ m_e × (4/3) × α × T₁₇

TAG: (D*H) — Direct TriPhase derivation with heuristic elements
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

# ========== SYMPLECTIC DERIVATION ==========
print("=" * 70)
print("TriPhase V16: Down Quark Mass (Symplectic)")
print("=" * 70)
print()

print("PHASE SPACE STRUCTURE")
print("-" * 70)
print("Quark phase space: (x, p, color)")
print("  x = position")
print("  p = momentum")
print("  color = SU(3) color charge (r, g, b)")
print("Symplectic form: ω = dp ∧ dx + ω_color")
print()

print("QUARK PROPERTIES")
print("-" * 70)
print("Down quark:")
print("  Charge: -1/3 e")
print("  Spin: 1/2")
print("  Color: r, g, or b (superposition)")
print("  Generation: 1st (lightest down-type)")
print()

print("ISOSPIN DOUBLET")
print("-" * 70)
print("Up and down form an isospin doublet:")
print("  u: charge +2/3, I₃ = +1/2")
print("  d: charge -1/3, I₃ = -1/2")
print()
print("Approximate isospin symmetry: m_u ≈ m_d (up to EM corrections)")
print("Actual: m_d > m_u (isospin breaking)")
print()

print("QUARK CONFINEMENT")
print("-" * 70)
print("Quarks are confined by QCD — never observed free")
print("Two mass scales:")
print("  1. Current mass (bare): m_d ~ few MeV")
print("  2. Constituent mass (in hadrons): M_d ~ 300 MeV")
print()
print("TriPhase predicts current mass (chiral symmetry breaking scale)")
print()

print("COLOR PHASE SPACE")
print("-" * 70)
print("SU(3) color group: r, g, b")
print("Gluons mediate strong force (8 gluon types)")
print("Poisson bracket: {color_i, color_j} ~ f_ijk color_k")
print()

print("TRIANGULAR HARMONIC")
print("-" * 70)
print(f"T₁₇ = 17 × 18 / 2 = {T_17}")
print("T₁₇ appears in quark mass formulas")
print()

print("CHARGE FRACTION SCALING")
print("-" * 70)
charge_fraction = 4.0 / 3.0
print(f"Down quark charge: -1/3 e")
print(f"Scaling factor: 4/3 (includes symplectic charge structure)")
print(f"Charge factor: {charge_fraction:.6f}")
print()

print("TRIPHASE DERIVATION")
print("-" * 70)
print(f"m_d ~ m_e × (4/3) × α × T₁₇")
print(f"")
print(f"m_e   = {m_e:.12e} kg")
print(f"α     = {alpha:.12e}")
print(f"T₁₇   = {T_17}")
print(f"4/3   = {charge_fraction:.6f}")
print(f"")

m_d = m_e * charge_fraction * alpha * T_17

print(f"m_d   = {m_d:.12e} kg")
print()

# Rest energy
E_d_J = m_d * c**2
E_d_eV = E_d_J / e
E_d_MeV = E_d_eV / 1e6

print(f"Rest energy (current mass):")
print(f"  E = m_d c² = {E_d_J:.6e} J")
print(f"  E = {E_d_eV:.3e} eV")
print(f"  E = {E_d_MeV:.6f} MeV")
print()

# Mass ratio
mass_ratio = m_d / m_e
print(f"Mass ratio: m_d / m_e = {mass_ratio:.6f}")
print()

# Up-down mass ratio
m_u = m_e * (2.0/3.0) * alpha * T_17
ratio_d_u = m_d / m_u
print(f"m_d / m_u ratio = {ratio_d_u:.3f} = (4/3) / (2/3) = 2.0")
print()

# Compton wavelength
lambda_C_d = h / (m_d * c)
lambda_C_d_fm = lambda_C_d * 1e15  # femtometers

print(f"Compton wavelength:")
print(f"  λ_C = {lambda_C_d:.6e} m")
print(f"  λ_C = {lambda_C_d_fm:.6f} fm")
print()

# ========== CALIBRATION CHECKPOINT ==========
# MS-bar mass at 2 GeV scale
m_d_current_MeV = 4.67  # MeV (PDG MS-bar mass)
m_d_current_kg = m_d_current_MeV * 1e6 * e / c**2

deviation_percent = (m_d - m_d_current_kg) / m_d_current_kg * 100.0

print("CALIBRATION CHECKPOINT")
print("-" * 70)
print(f"TriPhase m_d:    {E_d_MeV:.6f} MeV")
print(f"PDG m_d (MS̄):    {m_d_current_MeV:.2f} MeV (current mass, 2 GeV scale)")
print(f"Deviation:       {deviation_percent:+.1f} %")
print()
print("Note: Quark masses are scale-dependent (running masses in QCD).")
print("PDG quotes MS-bar mass at μ = 2 GeV renormalization scale.")
print()

print("CONSTITUENT VS CURRENT MASS")
print("-" * 70)
M_d_constituent_MeV = 300  # Typical constituent mass
print(f"Current mass m_d:      ~{m_d_current_MeV:.1f} MeV (TriPhase: {E_d_MeV:.1f} MeV)")
print(f"Constituent mass M_d:  ~{M_d_constituent_MeV} MeV (effective in hadrons)")
print()
print("The constituent mass includes QCD binding energy.")
print()

print("NEUTRON COMPOSITION")
print("-" * 70)
print("Neutron: udd (one up quark, two down quarks)")
m_n_measured_MeV = 939.565  # Neutron mass
print(f"Neutron mass: m_n ~ {m_n_measured_MeV:.1f} MeV")
print(f"Quark masses: m_u + 2m_d ~ {m_u*c**2/e/1e6 + 2*m_d*c**2/e/1e6:.1f} MeV (current)")
print()
print("Like the proton, most neutron mass comes from QCD binding energy.")
print()

print("NEUTRON-PROTON MASS DIFFERENCE")
print("-" * 70)
m_n_kg = m_n_measured_MeV * 1e6 * e / c**2
Delta_m_MeV = (m_n_kg - m_p) * c**2 / (e * 1e6)
print(f"Δm = m_n - m_p = {Delta_m_MeV:.3f} MeV")
print()
print("This mass difference is crucial for nuclear stability:")
print("  - Free neutron decays: n → p + e⁻ + ν̄_e")
print("  - Lifetime: τ_n ≈ 880 seconds")
print()
print("The mass difference arises from:")
print("  1. Quark mass difference: m_d - m_u")
print("  2. Electromagnetic self-energy (isospin breaking)")
print()

print("SYMPLECTIC GEOMETRY INSIGHT")
print("-" * 70)
print("The down quark mass emerges from the same triangular symplectic")
print("structure as the up quark, with a different charge-related scaling:")
print()
print("  m_u ~ m_e × (2/3) × α × T₁₇")
print("  m_d ~ m_e × (4/3) × α × T₁₇")
print()
print("The ratio m_d/m_u = 2 is built into the TriPhase charge structure.")
print("This suggests that the 'charge fractions' (2/3 for u, 4/3 for d)")
print("encode the symplectic geometry of quark phase space.")
print()
print("In symplectic geometry, both quarks live in extended phase space")
print("(x, p, color) with SU(3) color symmetry. The charge-dependent")
print("factors (2/3, 4/3) arise from the way electric charge couples to")
print("the triangular phase space lattice T₁₇.")
print()
print("The isospin doublet (u, d) represents a symplectic pair in flavor")
print("phase space, related by an SU(2) isospin transformation. The mass")
print("splitting m_d > m_u breaks the isospin symmetry, revealing the")
print("underlying electromagnetic structure of the TriPhase lattice.")
print()
print("This unifies up and down quarks within a single symplectic framework,")
print("where both emerge as excitations in the triangular phase space,")
print("scaled by their respective charge couplings to the EM field.")
print()
print("=" * 70)

input("Press Enter to exit...")
