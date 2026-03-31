#!/usr/bin/env python3
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
UP QUARK MASS FROM GROUP THEORY (D*H)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DERIVATION TAG: (D*H) = Pure derivation with hypothetical group structure

GROUP THEORY INTERPRETATION:
The up quark sits in multiple representations simultaneously:
  - SU(3)_color: Fundamental 3 representation
  - SU(2)_L: Weak isospin doublet (u, d)_L
  - U(1)_Y: Hypercharge Y = 1/3

Its mass comes from the PRODUCT of Casimir eigenvalues across all groups.

CASIMIR OPERATORS:
  - SU(3)_c: C₂(3) = (N²-1)/(2N) = 4/3 for fundamental
  - SU(2)_L: C₂(1/2) = 3/4 for doublet
  - U(1)_Y: C₁ = Y² = (1/3)² = 1/9

MASS FORMULA:
    m_u ~ m_e × α × √[C₂(SU(3)) × C₂(SU(2)) × C₁(U(1))]
    m_u ~ m_e × α × √[(4/3) × (3/4) × (1/9)]
    m_u ~ m_e × α / 3

The up quark is the LIGHTEST quark because it's in the smallest non-trivial
representations of all three gauge groups.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

import math

def derive_up_quark_mass_group_theory():
    """
    Derive up quark mass from SU(3)×SU(2)×U(1) representations.

    GROUP THEORY FRAMEWORK:
    - SU(3)_color: 3 (fundamental), Casimir = 4/3
    - SU(2)_weak: 2 (doublet), Casimir = 3/4
    - U(1)_hypercharge: Y = 1/3, Casimir = 1/9

    DERIVATION STRATEGY:
    1. Start from electron mass
    2. Apply Casimir factors from each gauge group
    3. Include QCD binding corrections
    4. Account for chiral symmetry breaking
    """

    print("=" * 80)
    print("UP QUARK MASS FROM GROUP THEORY (SU(3)×SU(2)×U(1))")
    print("=" * 80)
    print()
    print("GAUGE GROUPS:")
    print("  - SU(3)_color: fundamental 3, C₂ = 4/3")
    print("  - SU(2)_L: doublet 2, C₂ = 3/4")
    print("  - U(1)_Y: Y = 1/3, C₁ = 1/9")
    print()
    print("MASS FORMULA: m_u ~ m_e × α × √(C₂(3) × C₂(2) × C₁(U(1)))")
    print()
    print("=" * 80)
    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 1: ANCHOR CONSTANTS
    # ═══════════════════════════════════════════════════════════════════════════

    print("STEP 1: FUNDAMENTAL CONSTANTS")
    print("-" * 80)

    epsilon_0 = 8.8541878128e-12
    mu_0      = 1.25663706212e-6
    e         = 1.602176634e-19
    r_e       = 2.8179403262e-15

    print(f"  ε₀ = {epsilon_0:.13e} F/m")
    print(f"  μ₀ = {mu_0:.11e} H/m")
    print(f"  e  = {e:.12e} C")
    print(f"  r_e = {r_e:.10e} m")
    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 2: DERIVE STANDARD ANCHOR CHAIN
    # ═══════════════════════════════════════════════════════════════════════════

    print("STEP 2: STANDARD ANCHOR CHAIN")
    print("-" * 80)

    c = 1.0 / math.sqrt(epsilon_0 * mu_0)
    Z_0 = math.sqrt(mu_0 / epsilon_0)
    alpha_inv = 137.0 + math.log(137.0) / 137.0
    alpha = 1.0 / alpha_inv
    hbar = Z_0 * e**2 / (4.0 * math.pi * alpha)
    h = 2.0 * math.pi * hbar
    m_e = hbar * alpha / (c * r_e)
    f_e = m_e * c**2 / hbar

    print(f"  c     = {c:.10e} m/s")
    print(f"  α     = {alpha:.15f}")
    print(f"  ℏ     = {hbar:.15e} J·s")
    print(f"  m_e   = {m_e:.15e} kg")
    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 3: SU(3) COLOR CASIMIR
    # ═══════════════════════════════════════════════════════════════════════════

    print("STEP 3: SU(3)_color CASIMIR (FUNDAMENTAL REPRESENTATION)")
    print("-" * 80)

    # SU(N) fundamental: C₂(N) = (N²-1)/(2N)
    N_color = 3
    C2_SU3 = (N_color**2 - 1) / (2.0 * N_color)

    print(f"  SU(3) fundamental representation")
    print(f"  N = {N_color}")
    print(f"  C₂(3) = (N²-1)/(2N) = ({N_color}²-1)/(2×{N_color})")
    print(f"        = {C2_SU3:.10f}")
    print()
    print("INTERPRETATION:")
    print("  - Up quark carries color charge (red, green, or blue)")
    print("  - Fundamental 3 representation of SU(3)")
    print("  - Casimir = 4/3 is the color charge squared")
    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 4: SU(2) WEAK ISOSPIN CASIMIR
    # ═══════════════════════════════════════════════════════════════════════════

    print("STEP 4: SU(2)_L WEAK ISOSPIN CASIMIR (DOUBLET)")
    print("-" * 80)

    # SU(2) doublet: C₂(j) = j(j+1) with j=1/2
    j_weak = 0.5
    C2_SU2 = j_weak * (j_weak + 1.0)

    print(f"  SU(2)_L doublet: (u, d)_L")
    print(f"  Isospin j = {j_weak}")
    print(f"  C₂(1/2) = j(j+1) = {j_weak} × {j_weak + 1}")
    print(f"          = {C2_SU2:.10f}")
    print()
    print("INTERPRETATION:")
    print("  - Up quark is in weak isospin doublet with down quark")
    print("  - T₃(u) = +1/2, T₃(d) = -1/2")
    print("  - Casimir = 3/4 from isospin algebra")
    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 5: U(1) HYPERCHARGE CASIMIR
    # ═══════════════════════════════════════════════════════════════════════════

    print("STEP 5: U(1)_Y HYPERCHARGE CASIMIR")
    print("-" * 80)

    # Hypercharge for left-handed quark doublet: Y = 1/3
    Y_quark = 1.0 / 3.0
    C1_U1 = Y_quark**2

    print(f"  Left-handed quark doublet (u, d)_L")
    print(f"  Hypercharge Y = {Y_quark:.10f}")
    print(f"  C₁(U(1)) = Y² = ({Y_quark:.4f})²")
    print(f"           = {C1_U1:.10f}")
    print()
    print("INTERPRETATION:")
    print("  - Hypercharge Y relates to electric charge:")
    print("    Q = T₃ + Y/2")
    print("  - For up quark: Q = +1/2 + 1/6 = +2/3")
    print("  - For down quark: Q = -1/2 + 1/6 = -1/3")
    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 6: COMBINED CASIMIR FACTOR
    # ═══════════════════════════════════════════════════════════════════════════

    print("STEP 6: COMBINED CASIMIR FACTOR")
    print("-" * 80)

    # Product of Casimirs from all groups
    casimir_product = C2_SU3 * C2_SU2 * C1_U1
    casimir_factor = math.sqrt(casimir_product)

    print(f"  Product: C₂(SU(3)) × C₂(SU(2)) × C₁(U(1))")
    print(f"         = {C2_SU3:.6f} × {C2_SU2:.6f} × {C1_U1:.6f}")
    print(f"         = {casimir_product:.10f}")
    print()
    print(f"  Casimir factor = √(product)")
    print(f"                 = {casimir_factor:.10f}")
    print()
    print("GROUP THEORY INTERPRETATION:")
    print("  - Quark sits in tensor product: 3 ⊗ 2 ⊗ Y")
    print("  - Mass scale from combined Casimir eigenvalues")
    print("  - Factor √(product) from quantum fluctuations")
    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 7: QCD COUPLING CORRECTION
    # ═══════════════════════════════════════════════════════════════════════════

    print("STEP 7: QCD COUPLING CORRECTION")
    print("-" * 80)

    # QCD coupling α_s at low energy ~ 0.3
    # Provides additional mass scale
    alpha_s = 0.3  # Approximate strong coupling at 1 GeV

    qcd_correction = 1.0 + alpha_s / math.pi

    print(f"  α_s(1 GeV) ≈ {alpha_s}")
    print(f"  QCD correction = 1 + α_s/π")
    print(f"                 = {qcd_correction:.10f}")
    print()
    print("INTERPRETATION:")
    print("  - Strong coupling α_s >> α (electromagnetic)")
    print("  - Gluon loops provide mass enhancement")
    print("  - Running coupling depends on energy scale")
    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 8: CHIRAL SYMMETRY BREAKING
    # ═══════════════════════════════════════════════════════════════════════════

    print("STEP 8: CHIRAL SYMMETRY BREAKING CORRECTION")
    print("-" * 80)

    # QCD spontaneously breaks chiral symmetry
    # Generates constituent quark mass ~ Λ_QCD
    # Current quark mass (from Yukawa) is much smaller
    # Correction factor to get current mass from constituent

    chiral_factor = 0.25  # m_current / m_constituent

    print(f"  Chiral symmetry broken by QCD vacuum")
    print(f"  Constituent mass ~ Λ_QCD ~ 300 MeV")
    print(f"  Current mass ~ Yukawa coupling × v_Higgs")
    print(f"  Chiral factor = m_current / m_constituent")
    print(f"                = {chiral_factor:.6f}")
    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 9: COMPUTE UP QUARK MASS
    # ═══════════════════════════════════════════════════════════════════════════

    print("STEP 9: COMPUTE UP QUARK MASS")
    print("-" * 80)

    # Base formula: m_u ~ m_e × α × Casimir × corrections
    m_u = m_e * alpha * casimir_factor * qcd_correction * chiral_factor

    print(f"  m_u = m_e × α × √(Casimirs) × QCD × chiral")
    print(f"      = {m_e:.6e} × {alpha:.6f} × {casimir_factor:.6f}")
    print(f"        × {qcd_correction:.6f} × {chiral_factor:.6f}")
    print(f"      = {m_u:.15e} kg")
    print()

    # Convert to MeV/c²
    m_u_MeV = m_u * c**2 / (e * 1e6)

    print(f"  m_u = {m_u_MeV:.10f} MeV/c²")
    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 10: CALIBRATION CHECKPOINT
    # ═══════════════════════════════════════════════════════════════════════════

    print("STEP 10: CALIBRATION CHECKPOINT")
    print("-" * 80)

    m_u_pdg = 2.16  # MeV/c² (PDG 2020, MS-bar at 2 GeV)
    diff = m_u_MeV - m_u_pdg
    rel_err = (m_u_MeV - m_u_pdg) / m_u_pdg * 100

    print(f"  Derived:  m_u = {m_u_MeV:.10f} MeV/c²")
    print(f"  PDG:      m_u = {m_u_pdg:.10f} MeV/c² (MS-bar, 2 GeV)")
    print(f"  Diff:          {diff:.6f} MeV/c²")
    print(f"  Rel. err:      {rel_err:.3f}%")
    print()

    if abs(rel_err) < 10.0:
        print("  ✓ REASONABLE AGREEMENT (< 10%)")
    elif abs(rel_err) < 30.0:
        print("  ✓ ORDER OF MAGNITUDE CORRECT")
    else:
        print("  ⚠ SIGNIFICANT DEVIATION")
    print()
    print("NOTE: Quark masses are scheme-dependent (MS-bar vs pole)")
    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # SUMMARY
    # ═══════════════════════════════════════════════════════════════════════════

    print("=" * 80)
    print("SUMMARY: UP QUARK MASS FROM GROUP THEORY")
    print("=" * 80)
    print()
    print("GAUGE GROUP QUANTUM NUMBERS:")
    print("  - SU(3)_color: Fundamental 3 (red, green, blue)")
    print("  - SU(2)_L: Doublet with down quark, T₃ = +1/2")
    print("  - U(1)_Y: Hypercharge Y = +1/3")
    print("  - U(1)_EM: Electric charge Q = +2/3")
    print()
    print("CASIMIR EIGENVALUES:")
    print(f"  - C₂(SU(3)) = {C2_SU3:.4f} (color charge²)")
    print(f"  - C₂(SU(2)) = {C2_SU2:.4f} (isospin²)")
    print(f"  - C₁(U(1))  = {C1_U1:.4f} (hypercharge²)")
    print()
    print("MASS GENERATION:")
    print("  1. Yukawa coupling to Higgs: y_u × v_Higgs")
    print("  2. Group theory: mass ~ √(product of Casimirs)")
    print("  3. QCD corrections from gluon loops")
    print("  4. Chiral symmetry breaking: m_current << m_constituent")
    print()
    print("WHY UP QUARK IS LIGHTEST:")
    print("  - Smallest non-trivial representations of all groups")
    print("  - SU(3): fundamental 3 (not adjoint 8)")
    print("  - SU(2): doublet 2 (not triplet 3)")
    print("  - U(1): Y = 1/3 (minimal hypercharge for quarks)")
    print()
    print("COMPARISON TO LEPTONS:")
    print(f"  - Electron: no color, m_e ~ {m_e * c**2 / (e * 1e6):.4f} MeV")
    print(f"  - Up quark: has color, m_u ~ {m_u_MeV:.2f} MeV")
    print(f"  - Ratio m_u/m_e ~ {m_u_MeV / (m_e * c**2 / (e * 1e6)):.1f}")
    print()
    print("DERIVED VALUE:")
    print(f"  m_u = {m_u_MeV:.10f} MeV/c²")
    print(f"  PDG: {m_u_pdg:.10f} MeV/c² (MS-bar, 2 GeV)")
    print(f"  Agreement: {100 - abs(rel_err):.1f}%")
    print()
    print("=" * 80)
    print()

    return m_u

if __name__ == "__main__":
    m_u = derive_up_quark_mass_group_theory()
    input("Press Enter to exit...")
