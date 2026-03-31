"""
TriPhase V16 — Proton Mass (Statistical Mechanics Framework)
=============================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
The proton mass (~938.27 MeV) is one of the most profound emergent phenomena in
particle physics. Unlike the Higgs mechanism which gives mass to elementary quarks
and leptons, 99% of the proton's mass arises from QCD binding energy — the gluon
and quark kinetic energies plus the confining potential energy. In statistical
mechanics terms, the proton mass is the free energy of a three-quark system (uud)
in the confined phase of QCD, computed at zero temperature and zero chemical potential.
The partition function sums over all possible gluon field configurations binding the
valence quarks, including contributions from virtual quark-antiquark pairs (sea quarks)
and gluon self-interactions.

The confinement phase transition at T_c ~ 170 MeV provides the key statistical
mechanics context. Below T_c, the QCD partition function exhibits color confinement
with a linearly rising potential V(r) ~ σr, where σ ~ 1 GeV/fm is the string tension.
The proton mass emerges as the ground state energy of this system, computed via
lattice QCD Monte Carlo simulations that sample the path integral. In TriPhase, the
proton mass is derived from the electron mass via the ratio m_p/m_e ≈ 1836.15, which
encodes the electromagnetic and strong coupling hierarchies. The factor 4×27×17 captures
color (3²), flavor, and geometric factors in the three-body bound state problem.

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
print("TriPhase V16: Proton Mass (Statistical Mechanics)")
print("=" * 70)
print()

print("STATISTICAL MECHANICS FRAMEWORK")
print("--------------------------------")
print("Ensemble: Microcanonical (isolated three-quark system)")
print("Microstates: Gluon + sea quark configurations in confined phase")
print("Partition function: Z = ∫Dg exp(-S_QCD[g, u, d])")
print("Observable: Ground state energy E₀ = m_p c²")
print("Phase: Below T_c ~ 170 MeV (confined, color singlet only)")
print()

print("QCD CONFINEMENT AND EMERGENT MASS")
print("----------------------------------")
print(f"Electron mass m_e = {m_e:.6e} kg")
print(f"Fine structure α  = {alpha:.10f}")
print()

# TriPhase formula for proton mass
print("TriPhase proton-to-electron mass ratio:")
print(f"  m_p/m_e = 4 × 27 × 17 × (1 + 5α²/π)")
print()

factor_base = 4.0 * 27.0 * 17.0
factor_qed = 1.0 + 5.0 * alpha**2 / math.pi
mp_me_calc = factor_base * factor_qed

print(f"  Base factor (4 × 27 × 17)    = {factor_base:.4f}")
print(f"  QED correction (1 + 5α²/π)   = {factor_qed:.10f}")
print(f"  Total ratio                  = {mp_me_calc:.6f}")
print()

m_p_calc = m_e * mp_me_calc
print(f"m_p (TriPhase) = {m_p_calc:.15e} kg")
print()

# Convert to MeV/c²
m_p_MeV = m_p_calc * c**2 / (1.602176634e-13)
print(f"m_p (TriPhase) = {m_p_MeV:.6f} MeV/c²")
print()

# ========== CALIBRATION CHECKPOINT ==========
print("CALIBRATION CHECKPOINT")
print("----------------------")
# CODATA 2018: m_p = 1.67262192369(51)e-27 kg
m_p_codata = 1.67262192369e-27  # kg
m_p_codata_MeV = 938.27208816  # MeV/c²

deviation_kg = (m_p_calc - m_p_codata) / m_p_codata * 1e9
deviation_MeV = (m_p_MeV - m_p_codata_MeV) / m_p_codata_MeV * 1e9

print(f"CODATA 2018:        {m_p_codata:.15e} kg")
print(f"TriPhase:           {m_p_calc:.15e} kg")
print(f"Deviation:          {deviation_kg:.3f} ppb")
print()
print(f"CODATA 2018:        {m_p_codata_MeV:.8f} MeV/c²")
print(f"TriPhase:           {m_p_MeV:.8f} MeV/c²")
print(f"Deviation:          {deviation_MeV:.3f} ppb")
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT")
print("-----------------------------")
print("The proton mass is the canonical example of emergent mass from pure energy.")
print("In the QCD partition function at T = 0, the ground state consists of three")
print("valence quarks (uud) whose bare masses (~5 MeV total) contribute only 1% of")
print("the total. The remaining 99% comes from:")
print("  • Gluon field energy (kinetic + self-interaction)")
print("  • Quark kinetic energy (confined in ~1 fm volume)")
print("  • Sea quark vacuum polarization")
print("  • Confining potential energy (string tension)")
print()
print("The TriPhase formula m_p/m_e = 4×27×17×(1 + 5α²/π) encodes this:")
print("  • Factor 4: Two quarks (u) + antiquark (d̄) pairing structure")
print("  • Factor 27: 3³ color permutations in SU(3) gauge theory")
print("  • Factor 17: Geometric factor for three-body bound state")
print("  • QED correction: Electromagnetic contribution to binding energy")
print()
print("This ratio emerges from the density of states in the confined QCD vacuum,")
print("where the characteristic scale is set by Λ_QCD ~ 200 MeV. The partition")
print("function below T_c is dominated by colorless hadron states, with the proton")
print("being the lightest baryon. Statistical mechanics thus reveals the proton mass")
print("as a collective phenomenon — not a property of quarks, but of the QCD vacuum.")
print()
print("=" * 70)

input("Press Enter to exit...")
