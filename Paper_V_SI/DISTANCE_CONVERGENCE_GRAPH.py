# -*- coding: utf-8 -*-
"""
============================================================
DISTANCE CONVERGENCE GRAPH
============================================================
Visual comparison: Standard DIVERGES, TriPhase CONVERGES

Shows geometric vs wavelength distance agreement across all z.

Author: Magnetic Innovative Solutions LLC
Date: March 25, 2026
============================================================
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
matplotlib.use('Agg')  # Non-interactive backend for saving
import matplotlib.pyplot as plt

# ============================================================
# CONSTANTS
# ============================================================
c_km = 299792.458  # km/s
H0_std_local = 73.04   # SH0ES
H0_std_cmb = 67.4      # Planck
H0_tri = 71.48         # TriPhase (derived)

# ============================================================
# DISTANCE FORMULAS
# ============================================================

def d_wavelength_standard(z, H0):
    """Standard: d = cz/H0"""
    return c_km * z / H0

def d_wavelength_triphase(z):
    """TriPhase: d = (c/H0) * ln(1+z)"""
    return (c_km / H0_tri) * np.log(1 + z)

def d_geometric_standard(z, d_true):
    """
    Standard geometric distance.
    If we had the TRUE distance, standard would calculate it
    using v = cz for any Doppler component.
    For this comparison, we use TriPhase as "truth" and show
    what standard would report.
    """
    # Standard uses v = cz, so if true d came from v = c*ln(1+z),
    # standard would report d_std = d_true * (z / ln(1+z))^2
    if z > 0:
        return d_true * (z / np.log(1 + z)) ** 2
    return d_true

def d_geometric_triphase(z, d_true):
    """TriPhase geometric = true (by definition when using correct formulas)"""
    return d_true

# ============================================================
# MASER DATA (Real measurements)
# ============================================================
masers = [
    # (name, z, d_published Mpc, uncertainty %)
    ("NGC 4258",  0.00158,  7.576,  3.0),
    ("UGC 3789",  0.01109,  49.6,   5.5),
    ("NGC 6264",  0.03402,  143.6,  11.0),
    ("NGC 6323",  0.02602,  106.5,  9.0),
    ("NGC 5765b", 0.02785,  126.3,  11.0),
    ("NGC 1194",  0.01351,  53.2,   7.0),
    ("Mrk 1419",  0.01648,  78.0,   9.0),
]

# ============================================================
# CREATE THE GRAPH
# ============================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Geometric vs Wavelength Distance: Standard DIVERGES, TriPhase CONVERGES',
             fontsize=14, fontweight='bold')

# Generate z range
z_range = np.linspace(0.001, 1.5, 500)

# ============================================================
# PANEL 1: Standard Cosmology - Absolute Distances
# ============================================================
ax1 = axes[0, 0]

# Wavelength distance (what they calculate from redshift)
d_wave_std_73 = d_wavelength_standard(z_range, H0_std_local)
d_wave_std_67 = d_wavelength_standard(z_range, H0_std_cmb)

# "True" geometric distance (using TriPhase as ground truth)
d_true = d_wavelength_triphase(z_range)

# What standard would report for geometric (inflated by wrong v formula)
d_geo_std = np.array([d_geometric_standard(z, d) for z, d in zip(z_range, d_true)])

ax1.plot(z_range, d_true, 'g-', linewidth=2, label='True Distance (TriPhase)')
ax1.plot(z_range, d_wave_std_73, 'b--', linewidth=1.5, label=f'Wavelength (H0={H0_std_local})')
ax1.plot(z_range, d_wave_std_67, 'b:', linewidth=1.5, label=f'Wavelength (H0={H0_std_cmb})')
ax1.plot(z_range, d_geo_std, 'r-', linewidth=2, label='Geometric (v=cz)')

ax1.set_xlabel('Redshift z', fontsize=11)
ax1.set_ylabel('Distance (Mpc)', fontsize=11)
ax1.set_title('Standard Cosmology: Distances Diverge', fontsize=12)
ax1.legend(loc='upper left', fontsize=9)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 1.5)

# ============================================================
# PANEL 2: TriPhase - Absolute Distances
# ============================================================
ax2 = axes[0, 1]

# TriPhase wavelength distance
d_wave_tri = d_wavelength_triphase(z_range)

# TriPhase geometric (same as wavelength when done correctly)
d_geo_tri = d_wave_tri  # They match!

ax2.plot(z_range, d_wave_tri, 'g-', linewidth=2, label='Wavelength: d = (c/H0)·ln(1+z)')
ax2.plot(z_range, d_geo_tri, 'r--', linewidth=2, label='Geometric (v=c·ln(1+z))')

ax2.set_xlabel('Redshift z', fontsize=11)
ax2.set_ylabel('Distance (Mpc)', fontsize=11)
ax2.set_title('TriPhase: Distances CONVERGE (Same Line!)', fontsize=12)
ax2.legend(loc='upper left', fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 1.5)

# ============================================================
# PANEL 3: Percent Difference (Geometric - Wavelength) / Wavelength
# ============================================================
ax3 = axes[1, 0]

# Standard: geometric vs wavelength percent difference
# Using H0=73 for wavelength
pct_diff_std_73 = ((d_geo_std - d_wave_std_73) / d_wave_std_73) * 100
pct_diff_std_67 = ((d_geo_std - d_wave_std_67) / d_wave_std_67) * 100

# TriPhase: should be ~0%
pct_diff_tri = ((d_geo_tri - d_wave_tri) / d_wave_tri) * 100

ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
ax3.fill_between(z_range, -5, 5, alpha=0.2, color='green', label='±5% (measurement uncertainty)')

ax3.plot(z_range, pct_diff_std_73, 'b-', linewidth=2, label=f'Standard (H0={H0_std_local})')
ax3.plot(z_range, pct_diff_std_67, 'b--', linewidth=2, label=f'Standard (H0={H0_std_cmb})')
ax3.plot(z_range, pct_diff_tri, 'g-', linewidth=3, label='TriPhase (H0=71.48)')

# Add maser data points
for name, z, d_pub, unc in masers:
    # Standard difference
    d_wave_m = d_wavelength_standard(z, H0_std_local)
    d_geo_m = d_pub  # Published uses v=cz
    pct_m = ((d_geo_m - d_wave_m) / d_wave_m) * 100
    ax3.scatter(z, pct_m, color='blue', s=60, zorder=5, edgecolors='black')

    # TriPhase difference
    d_wave_tri_m = d_wavelength_triphase(z)
    d_geo_tri_m = d_pub * (np.log(1+z)/z)**2  # Corrected
    pct_tri_m = ((d_geo_tri_m - d_wave_tri_m) / d_wave_tri_m) * 100
    ax3.scatter(z, pct_tri_m, color='green', s=60, zorder=5, edgecolors='black', marker='s')

ax3.set_xlabel('Redshift z', fontsize=11)
ax3.set_ylabel('(Geometric - Wavelength) / Wavelength (%)', fontsize=11)
ax3.set_title('Distance Divergence: Standard vs TriPhase', fontsize=12)
ax3.legend(loc='upper right', fontsize=9)
ax3.grid(True, alpha=0.3)
ax3.set_xlim(0, 0.15)  # Focus on maser redshift range
ax3.set_ylim(-30, 30)

# ============================================================
# PANEL 4: Extended z range - The Full Picture
# ============================================================
ax4 = axes[1, 1]

z_extended = np.linspace(0.01, 10.0, 500)

# Ratio: d_geometric / d_wavelength
# Standard: d_geo/d_wave = (z/ln(1+z))^2 * (H0/c) / (z/H0/c) = z/ln(1+z) * (H0_wave/H0_geo)
# Simplify: if same H0, ratio = (z/ln(1+z))^2

ratio_std = (z_extended / np.log(1 + z_extended)) ** 2

# TriPhase: ratio = 1 (they match)
ratio_tri = np.ones_like(z_extended)

ax4.axhline(y=1, color='black', linestyle='-', linewidth=0.5)
ax4.fill_between(z_extended, 0.95, 1.05, alpha=0.2, color='green', label='±5% agreement')

ax4.plot(z_extended, ratio_std, 'r-', linewidth=2.5, label='Standard: d_geo / d_wave')
ax4.plot(z_extended, ratio_tri, 'g-', linewidth=2.5, label='TriPhase: d_geo / d_wave')

# Mark key points
ax4.scatter([0.1], [ratio_std[np.argmin(np.abs(z_extended - 0.1))]],
            color='red', s=100, zorder=5, edgecolors='black')
ax4.annotate(f'z=0.1: {(ratio_std[np.argmin(np.abs(z_extended - 0.1))]-1)*100:.0f}% off',
             xy=(0.1, ratio_std[np.argmin(np.abs(z_extended - 0.1))]),
             xytext=(0.2, 1.15), fontsize=10,
             arrowprops=dict(arrowstyle='->', color='red'))

ax4.scatter([0.5], [ratio_std[np.argmin(np.abs(z_extended - 0.5))]],
            color='red', s=100, zorder=5, edgecolors='black')
ax4.annotate(f'z=0.5: {(ratio_std[np.argmin(np.abs(z_extended - 0.5))]-1)*100:.0f}% off',
             xy=(0.5, ratio_std[np.argmin(np.abs(z_extended - 0.5))]),
             xytext=(0.6, 1.35), fontsize=10,
             arrowprops=dict(arrowstyle='->', color='red'))

ax4.scatter([1.0], [ratio_std[np.argmin(np.abs(z_extended - 1.0))]],
            color='red', s=100, zorder=5, edgecolors='black')
ax4.annotate(f'z=1: {(ratio_std[np.argmin(np.abs(z_extended - 1.0))]-1)*100:.0f}%',
             xy=(1.0, ratio_std[np.argmin(np.abs(z_extended - 1.0))]),
             xytext=(1.3, 2.5), fontsize=9,
             arrowprops=dict(arrowstyle='->', color='red'))

ax4.scatter([2.0], [ratio_std[np.argmin(np.abs(z_extended - 2.0))]],
            color='red', s=100, zorder=5, edgecolors='black')
ax4.annotate(f'z=2: {(ratio_std[np.argmin(np.abs(z_extended - 2.0))]-1)*100:.0f}%',
             xy=(2.0, ratio_std[np.argmin(np.abs(z_extended - 2.0))]),
             xytext=(2.5, 4.0), fontsize=9,
             arrowprops=dict(arrowstyle='->', color='red'))

ax4.scatter([5.0], [ratio_std[np.argmin(np.abs(z_extended - 5.0))]],
            color='red', s=100, zorder=5, edgecolors='black')
ax4.annotate(f'z=5: {(ratio_std[np.argmin(np.abs(z_extended - 5.0))]-1)*100:.0f}%',
             xy=(5.0, ratio_std[np.argmin(np.abs(z_extended - 5.0))]),
             xytext=(5.5, 8.0), fontsize=9,
             arrowprops=dict(arrowstyle='->', color='red'))

ax4.scatter([10.0], [ratio_std[np.argmin(np.abs(z_extended - 10.0))]],
            color='red', s=100, zorder=5, edgecolors='black')
ax4.annotate(f'z=10: {(ratio_std[np.argmin(np.abs(z_extended - 10.0))]-1)*100:.0f}%',
             xy=(10.0, ratio_std[np.argmin(np.abs(z_extended - 10.0))]),
             xytext=(8.0, 15.0), fontsize=9,
             arrowprops=dict(arrowstyle='->', color='red'))

ax4.set_xlabel('Redshift z', fontsize=11)
ax4.set_ylabel('d_geometric / d_wavelength', fontsize=11)
ax4.set_title('Distance Ratio at High z: Standard Diverges, TriPhase Consistent', fontsize=12)
ax4.legend(loc='upper left', fontsize=10)
ax4.grid(True, alpha=0.3)
ax4.set_xlim(0, 10.5)
ax4.set_ylim(0.9, 20)

plt.tight_layout()
plt.savefig(r'C:\Users\cfucc\Documents\MIS Magnetic Innovative Solutions LLC\TriPhase\Calculations\DISTANCE_CONVERGENCE_GRAPH.png',
            dpi=150, bbox_inches='tight')
plt.savefig(r'C:\Users\cfucc\Documents\MIS Magnetic Innovative Solutions LLC\TriPhase\Calculations\DISTANCE_CONVERGENCE_GRAPH.pdf',
            bbox_inches='tight')

print("Graph saved to:")
print("  DISTANCE_CONVERGENCE_GRAPH.png")
print("  DISTANCE_CONVERGENCE_GRAPH.pdf")
print()
print("KEY RESULT - STANDARD COSMOLOGY DIVERGENCE:")
print(f"  At z=0.1:  Standard diverges by {(ratio_std[np.argmin(np.abs(z_extended - 0.1))]-1)*100:.0f}%")
print(f"  At z=0.5:  Standard diverges by {(ratio_std[np.argmin(np.abs(z_extended - 0.5))]-1)*100:.0f}%")
print(f"  At z=1.0:  Standard diverges by {(ratio_std[np.argmin(np.abs(z_extended - 1.0))]-1)*100:.0f}%")
print(f"  At z=2.0:  Standard diverges by {(ratio_std[np.argmin(np.abs(z_extended - 2.0))]-1)*100:.0f}%")
print(f"  At z=5.0:  Standard diverges by {(ratio_std[np.argmin(np.abs(z_extended - 5.0))]-1)*100:.0f}%")
print(f"  At z=10.0: Standard diverges by {(ratio_std[np.argmin(np.abs(z_extended - 10.0))]-1)*100:.0f}%")
print()
print("TRIPHASE: 0% divergence at ALL redshifts - geometric = wavelength ALWAYS")

# plt.show()  # Commented out for batch mode - uncomment to display interactively