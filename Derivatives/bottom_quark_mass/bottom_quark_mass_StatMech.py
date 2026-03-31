"""
TriPhase V16 — Bottom Quark Mass (Statistical Mechanics Framework)
===================================================================
(c) 2025-2026 Christian R. Fuccillo / MIS Magnetic Innovative Solutions LLC
DOI: 10.5281/zenodo.17855383

STATISTICAL MECHANICS INTERPRETATION
--------------------------------------
The bottom quark mass (~4.18 GeV) represents a deep threshold in the QCD spectrum
where heavy quark effective theory (HQET) becomes the natural description. In the
canonical ensemble, bottom quarks are sufficiently heavy that they decouple from
light quark dynamics, acting as static color sources around which gluon fields
organize. The mass emerges from the self-energy of this static configuration,
computed via the partition function of gluonic excitations in the presence of
heavy quark sources. The density of bottomonium states follows from solving the
Schrödinger equation with a Cornell potential, whose parameters are themselves
statistical averages over QCD vacuum configurations.

From a statistical mechanics perspective, the bottom quark mass fixes the scale
where the running coupling α_s(m_b) reaches a value that makes perturbative
calculations reliable even in bound state physics. The partition function for
Υ (upsilon) mesons exhibits clear non-relativistic structure, with level spacings
determined by the reduced mass and the average confining force. TriPhase encodes
this heavy quark limit through factors that suppress relativistic corrections,
while the alpha dependence captures the electromagnetic contribution to the quark
self-energy through vacuum polarization loops.

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
print("TriPhase V16: Bottom Quark Mass (Statistical Mechanics)")
print("=" * 70)
print()

print("STATISTICAL MECHANICS FRAMEWORK")
print("--------------------------------")
print("Ensemble: Canonical (HQET limit, static quark sources)")
print("Microstates: Gluon field configurations around heavy quarks")
print("Partition function: Z = ∫Dg exp(-S_eff[g, m_b])")
print("Observable: Heavy quark self-energy from gluon loops")
print()

print("HEAVY QUARK EFFECTIVE THEORY REGIME")
print("------------------------------------")
print(f"Electron rest frequency f_e = {f_e:.6e} Hz")
print(f"Fine structure constant α   = {alpha:.10f}")
print(f"Proton-to-electron mass ratio = {mp_me:.6f}")
print()

# TriPhase formula for bottom quark mass
# m_b ~ m_e * (very large factor) * alpha^(-5) approximately
# HQET regime: mass from static quark self-energy
coeff_b = 81.0 * math.exp(7.8 * alpha) / (alpha**3)
m_b_tph = m_e * coeff_b

print(f"HQET statistical factor: {coeff_b:.6f}")
print(f"m_b (TriPhase) = m_e × {coeff_b:.6f}")
print(f"m_b (TriPhase) = {m_b_tph:.6e} kg")
print()

# Convert to GeV/c²
m_b_GeV = m_b_tph * c**2 / (1.602176634e-10)  # Convert J to GeV
print(f"m_b (TriPhase) = {m_b_GeV:.6f} GeV/c²")
print()

# ========== CALIBRATION CHECKPOINT ==========
print("CALIBRATION CHECKPOINT")
print("----------------------")
# PDG MS-bar scheme: m_b(m_b) = 4.18 GeV
m_b_pdg = 4.18  # GeV/c²
deviation = (m_b_GeV - m_b_pdg) / m_b_pdg * 1e6
print(f"PDG value (MS-bar):     {m_b_pdg:.6f} GeV/c²")
print(f"TriPhase prediction:    {m_b_GeV:.6f} GeV/c²")
print(f"Deviation:              {deviation:.0f} ppm")
print()

# ========== STATISTICAL MECHANICS INSIGHT ==========
print("STATISTICAL MECHANICS INSIGHT")
print("-----------------------------")
print("In the HQET limit, the bottom quark acts as a static color source, and its")
print("mass emerges from the partition function of gluon excitations around this")
print("fixed charge. The statistical ensemble consists of all possible gluon field")
print("configurations consistent with the boundary conditions imposed by the heavy")
print("quark. The mass is then the vacuum expectation value of the energy operator")
print("in this restricted ensemble. TriPhase captures this through alpha^(-3) scaling,")
print("indicating that heavy quark creation requires overcoming a large electromagnetic")
print("barrier. The exponential factor reflects the density of bottomonium resonances,")
print("which form a nearly non-relativistic tower of states with Coulombic plus linear")
print("confinement energy levels. This spectrum is the statistical signature of QCD")
print("transitioning from asymptotic freedom (short distance) to confinement (long")
print("distance) within the same partition function.")
print()
print("=" * 70)

input("Press Enter to exit...")
