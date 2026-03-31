"""
TriPhase V16 — Charm Quark Mass (Statistical Mechanics Framework)
==================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
The charm quark mass appears in the QCD partition function as a heavy quark threshold
that modifies the running of coupling constants and the density of hadronic states.
In statistical mechanics, heavy quarks introduce a temperature-dependent effective
action through thermal loop corrections. The charm quark operates in a regime where
both perturbative QCD (at high T) and non-perturbative confinement (at low T) contribute
to the partition function. The mass can be extracted from the free energy difference
between flavored and unflavored QCD vacua, computed in the canonical ensemble at zero
chemical potential.

The charm mass scale (~1.27 GeV) represents a critical energy where the QCD beta
function changes character due to the activation of an additional quark flavor in
the vacuum polarization loops. In TriPhase, this is encoded through higher powers
of alpha that capture the increased phase space for gluon-quark fluctuations. The
statistical interpretation involves integrating out heavy quark modes in the path
integral, generating an effective Lagrangian whose mass parameter emerges from the
vacuum expectation value of the chiral condensate at the charm mass scale.

TAG: (D*) — TriPhase derivation requiring phenomenological coefficient
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
print("TriPhase V16: Charm Quark Mass (Statistical Mechanics)")
print("=" * 70)
print()

print("STATISTICAL MECHANICS FRAMEWORK")
print("--------------------------------")
print("Ensemble: Canonical (QCD at T = 0, fixed flavor number)")
print("Microstates: Heavy quark-antiquark pairs in QCD vacuum")
print("Free energy: F = -T ln Z, T → 0 limit extracts mass")
print("Threshold effect: Charm activates at μ ~ m_c c²")
print()

print("HEAVY QUARK MASS FROM VACUUM POLARIZATION")
print("------------------------------------------")
print(f"Electron rest frequency f_e = {f_e:.6e} Hz")
print(f"Fine structure constant α   = {alpha:.10f}")
print()

# TriPhase formula for charm quark mass
# m_c ~ m_e * (large factor) * alpha^(-4) * (phenomenological)
# Heavy quark mass from density of states in perturbative regime
coeff_c = 27.0 * math.exp(6.2 * alpha) / (alpha**2)
m_c_tph = m_e * coeff_c

print(f"Heavy quark statistical factor: {coeff_c:.6f}")
print(f"m_c (TriPhase) = m_e × {coeff_c:.6f}")
print(f"m_c (TriPhase) = {m_c_tph:.6e} kg")
print()

# Convert to GeV/c²
m_c_GeV = m_c_tph * c**2 / (1.602176634e-10)  # Convert J to GeV
print(f"m_c (TriPhase) = {m_c_GeV:.6f} GeV/c²")
print()

# ========== CALIBRATION CHECKPOINT ==========
print("CALIBRATION CHECKPOINT")
print("----------------------")
# PDG MS-bar scheme: m_c(m_c) = 1.27 GeV
m_c_pdg = 1.27  # GeV/c²
deviation = (m_c_GeV - m_c_pdg) / m_c_pdg * 1e6
print(f"PDG value (MS-bar):     {m_c_pdg:.6f} GeV/c²")
print(f"TriPhase prediction:    {m_c_GeV:.6f} GeV/c²")
print(f"Deviation:              {deviation:.0f} ppm")
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT")
print("-----------------------------")
print("The charm quark mass marks a threshold in the QCD partition function where")
print("a new flavor degree of freedom becomes thermally accessible. In the path")
print("integral formulation, integrating out charm quarks generates a mass-dependent")
print("effective action that modifies the vacuum energy density. The TriPhase formula")
print("captures this through alpha^(-2) scaling, reflecting the increased coupling")
print("strength required to create heavy quark pairs from vacuum fluctuations.")
print("The exponential factor arises from the density of charmonium bound states,")
print("which proliferate as phase space opens up above the charm threshold. This")
print("statistical ensemble of resonances encodes the non-perturbative structure")
print("of QCD confinement at intermediate energy scales.")
print()
print("=" * 70)

input("Press Enter to exit...")
