# -*- coding: utf-8 -*-
"""
Generate individual distance convergence graphs for Paper V
"""
import sys
import traceback
import atexit

atexit.register(lambda: input('Press Enter to exit...'))
sys.stdout.reconfigure(encoding='utf-8')

_original_excepthook = sys.excepthook
def _custom_excepthook(exc_type, exc_value, exc_tb):
    traceback.print_exception(exc_type, exc_value, exc_tb)
sys.excepthook = _custom_excepthook


import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Constants
c_km = 299792.458
H0_std_local = 73.04
H0_std_cmb = 67.4
H0_tri = 71.48

def d_wavelength_standard(z, H0):
    return c_km * z / H0

def d_wavelength_triphase(z):
    return (c_km / H0_tri) * np.log(1 + z)

# Maser data
masers = [
    ("NGC 4258",  0.00158,  7.576,  3.0),
    ("UGC 3789",  0.01109,  49.6,   5.5),
    ("NGC 6264",  0.03402,  143.6,  11.0),
    ("NGC 6323",  0.02602,  106.5,  9.0),
    ("NGC 5765b", 0.02785,  126.3,  11.0),
    ("NGC 1194",  0.01351,  53.2,   7.0),
    ("Mrk 1419",  0.01648,  78.0,   9.0),
]

z_range = np.linspace(0.001, 1.5, 500)
z_extended = np.linspace(0.01, 10.0, 500)

# ============================================================
# GRAPH 1: Distance Convergence Test - Extended z range showing divergence
# ============================================================
fig1, ax1 = plt.subplots(figsize=(10, 7))

ratio_std = (z_extended / np.log(1 + z_extended)) ** 2
ratio_tri = np.ones_like(z_extended)

ax1.axhline(y=1, color='black', linestyle='-', linewidth=0.5)
ax1.fill_between(z_extended, 0.95, 1.05, alpha=0.2, color='green', label='±5% agreement')

ax1.plot(z_extended, ratio_std, 'r-', linewidth=3, label='Standard: d = cz/H₀')
ax1.plot(z_extended, ratio_tri, 'g-', linewidth=3, label='Logarithmic: d = (c/H₀)ln(1+z)')

# Annotations
for z_val in [0.5, 1.0, 2.0, 5.0, 10.0]:
    idx = np.argmin(np.abs(z_extended - z_val))
    pct = (ratio_std[idx] - 1) * 100
    ax1.scatter([z_val], [ratio_std[idx]], color='red', s=100, zorder=5, edgecolors='black')
    if z_val <= 2:
        ax1.annotate(f'z={z_val}: {pct:.0f}%', xy=(z_val, ratio_std[idx]),
                     xytext=(z_val + 0.3, ratio_std[idx] + 0.5), fontsize=11,
                     arrowprops=dict(arrowstyle='->', color='red'))
    else:
        ax1.annotate(f'z={z_val}: {pct:.0f}%', xy=(z_val, ratio_std[idx]),
                     xytext=(z_val - 1, ratio_std[idx] + 2), fontsize=11,
                     arrowprops=dict(arrowstyle='->', color='red'))

ax1.set_xlabel('Redshift z', fontsize=14)
ax1.set_ylabel('Geometric Distance / Wavelength Distance', fontsize=14)
ax1.set_title('Distance Convergence Test: Linear vs Logarithmic', fontsize=14, fontweight='bold')
ax1.legend(loc='upper left', fontsize=12)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 10.5)
ax1.set_ylim(0.9, 20)

plt.tight_layout()
plt.savefig('Fig1_Smoking_Gun_Divergence.png', dpi=200, bbox_inches='tight')
plt.savefig('Fig1_Smoking_Gun_Divergence.pdf', bbox_inches='tight')
print("Saved: Fig1_Smoking_Gun_Divergence.png/pdf")
plt.close()

# ============================================================
# GRAPH 2: Maser data - Standard vs Framework at low z
# ============================================================
fig2, ax2 = plt.subplots(figsize=(10, 6))

z_low = np.linspace(0.001, 0.05, 200)
d_true = d_wavelength_triphase(z_low)
d_geo_std = d_true * (z_low / np.log(1 + z_low)) ** 2

# Percent difference curves
d_wave_std_73 = d_wavelength_standard(z_low, H0_std_local)
pct_diff_std = ((d_geo_std - d_wave_std_73) / d_wave_std_73) * 100
pct_diff_tri = np.zeros_like(z_low)

ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
ax2.fill_between(z_low, -5, 5, alpha=0.2, color='green', label='±5% measurement uncertainty')

ax2.plot(z_low, pct_diff_std, 'b-', linewidth=2.5, label=f'Standard: cz/H₀ (H₀={H0_std_local})')
ax2.plot(z_low, pct_diff_tri, 'g-', linewidth=2.5, label='Logarithmic: ln(1+z) (H₀=71.48)')

# Plot maser data points
for name, z, d_pub, unc in masers:
    d_wave_m = d_wavelength_standard(z, H0_std_local)
    pct_std = ((d_pub - d_wave_m) / d_wave_m) * 100
    ax2.scatter(z, pct_std, color='blue', s=80, zorder=5, edgecolors='black')

    d_wave_tri_m = d_wavelength_triphase(z)
    d_geo_tri_m = d_pub * (np.log(1+z)/z)**2
    pct_tri = ((d_geo_tri_m - d_wave_tri_m) / d_wave_tri_m) * 100
    ax2.scatter(z, pct_tri, color='green', s=80, zorder=5, edgecolors='black', marker='s')

ax2.set_xlabel('Redshift z', fontsize=13)
ax2.set_ylabel('(Geometric - Wavelength) / Wavelength (%)', fontsize=13)
ax2.set_title('Maser Distance Comparison: Linear vs Logarithmic', fontsize=13, fontweight='bold')
ax2.legend(loc='upper right', fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 0.05)
ax2.set_ylim(-25, 25)

plt.tight_layout()
plt.savefig('Fig2_Maser_Comparison.png', dpi=200, bbox_inches='tight')
plt.savefig('Fig2_Maser_Comparison.pdf', bbox_inches='tight')
print("Saved: Fig2_Maser_Comparison.png/pdf")
plt.close()

# ============================================================
# GRAPH 3: Correction factor table as visual
# ============================================================
fig3, ax3 = plt.subplots(figsize=(9, 6))

z_corr = np.array([0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
correction = (np.log(1 + z_corr) / z_corr) ** 2
pct_change = (1 - correction) * 100  # Positive = overestimate

ax3.bar(range(len(z_corr)), pct_change, color=['green' if p < 10 else 'orange' if p < 50 else 'red' for p in pct_change],
        edgecolor='black', linewidth=1.5)
ax3.set_xticks(range(len(z_corr)))
ax3.set_xticklabels([f'z={z}' for z in z_corr], fontsize=11)
ax3.set_ylabel('Standard Overestimates Distance By (%)', fontsize=13)
ax3.set_title('DISTANCE CORRECTION NEEDED: How Much Standard Cosmology Overestimates', fontsize=13, fontweight='bold')
ax3.axhline(y=5, color='gray', linestyle='--', linewidth=1, label='5% threshold')

for i, (z, pct) in enumerate(zip(z_corr, pct_change)):
    ax3.annotate(f'{pct:.0f}%', xy=(i, pct + 1), ha='center', fontsize=10, fontweight='bold')

ax3.grid(True, alpha=0.3, axis='y')
ax3.set_ylim(0, 90)

plt.tight_layout()
plt.savefig('Fig3_Correction_Factor.png', dpi=200, bbox_inches='tight')
plt.savefig('Fig3_Correction_Factor.pdf', bbox_inches='tight')
print("Saved: Fig3_Correction_Factor.png/pdf")
plt.close()

print("\nAll individual graphs generated!")
print("Place in Paper V sections as appropriate.")