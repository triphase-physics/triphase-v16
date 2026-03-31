"""
============================================================
GEOMETRIC vs REDSHIFT DISTANCE COMPARISON
============================================================
Compare actual geometric distances to redshift-derived distances
to see if TriPhase light slowdown can explain discrepancies.

Geometric methods (independent of cosmology):
- Megamasers: D from orbital geometry
- Eclipsing binaries: D from stellar radii + angular size
- Lensing time delays: D from light path geometry

Standard cosmology predicts D(z) from redshift.
If light slows down, standard D(z) is WRONG.
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

# =============================================================
# PHYSICAL CONSTANTS
# =============================================================
c = 299792.458  # km/s
H0_planck = 67.4  # km/s/Mpc - CMB-derived
H0_local = 73.0   # km/s/Mpc - Local distance ladder
H0_triphase = 71.4  # km/s/Mpc - TriPhase derived

# =============================================================
# TRIPHASE DRAG FORMULA
# =============================================================
def drag_triphase(z, n=5):
    """
    TriPhase drag formula from Bush-CMB connection.
    Drag_n(z) = (z/1100) * 5/(36 * n(n-1))

    At z=1100, n=5: Drag = 0.00694 ≈ 0.7% (matches Bush-CMB shift)
    """
    return (z / 1100) * 5 / (36 * n * (n - 1))

def cumulative_drag_factor(z, n=5):
    """
    Cumulative correction factor for distance.
    If light slows, it takes longer to arrive.
    During that extra time, more expansion occurs.
    So observed z is HIGHER than "true" expansion z.

    This means: D_true < D_standard for given observed z.
    """
    drag = drag_triphase(z, n)
    # First-order correction: D_true ≈ D_standard * (1 - drag)
    return 1 - drag

# =============================================================
# STANDARD COSMOLOGY D(z) - Simplified for low z
# =============================================================
def D_hubble(z, H0=H0_local):
    """Simple Hubble law distance for low z: D = cz/H0"""
    return c * z / H0  # Mpc

def D_comoving_flat(z, H0=H0_local, Om=0.3):
    """
    Comoving distance in flat ΛCDM.
    D_c = (c/H0) * integral[0 to z] dz'/E(z')
    where E(z) = sqrt(Om(1+z)^3 + (1-Om))

    For simplicity, use numerical integration.
    """
    from scipy.integrate import quad

    def E(zp):
        return np.sqrt(Om * (1 + zp)**3 + (1 - Om))

    integrand = lambda zp: 1 / E(zp)
    result, _ = quad(integrand, 0, z)
    return (c / H0) * result

def D_angular_diameter(z, H0=H0_local, Om=0.3):
    """Angular diameter distance: D_A = D_c / (1+z)"""
    D_c = D_comoving_flat(z, H0, Om)
    return D_c / (1 + z)

def D_luminosity(z, H0=H0_local, Om=0.3):
    """Luminosity distance: D_L = D_c * (1+z)"""
    D_c = D_comoving_flat(z, H0, Om)
    return D_c * (1 + z)

# =============================================================
# GEOMETRIC DISTANCE DATA
# =============================================================
# Format: name, z (redshift), D_geometric (Mpc), D_geo_err, method

GEOMETRIC_DATA = [
    # MEGAMASER DISTANCES (angular diameter distance from geometry)
    # Source: Megamaser Cosmology Project
    {
        'name': 'NGC 4258',
        'z': 470 / c,  # v_sys = 470 km/s
        'D_geo': 7.576,
        'D_geo_err': 0.11,
        'method': 'Megamaser',
        'ref': 'Reid+ 2019'
    },
    {
        'name': 'UGC 3789',
        'z': 3325 / c,  # v ~ 3325 km/s (from Hubble flow at ~50 Mpc with H0~67)
        'D_geo': 49.6,
        'D_geo_err': 5.1,
        'method': 'Megamaser',
        'ref': 'MCP IV'
    },
    {
        'name': 'NGC 6264',
        'z': 10200 / c,  # v ~ 10200 km/s (at ~144 Mpc)
        'D_geo': 144,
        'D_geo_err': 19,
        'method': 'Megamaser',
        'ref': 'Kuo+ 2013'
    },

    # ECLIPSING BINARY DISTANCES
    {
        'name': 'LMC',
        'z': 262 / c,  # v ~ 262 km/s
        'D_geo': 0.04997,  # 49.97 kpc in Mpc
        'D_geo_err': 0.00111,
        'method': 'Eclipsing Binary',
        'ref': 'Pietrzyński+ 2019'
    },
    {
        'name': 'SMC',
        'z': 158 / c,  # v ~ 158 km/s
        'D_geo': 0.0621,  # ~62 kpc
        'D_geo_err': 0.002,
        'method': 'Eclipsing Binary',
        'ref': 'Graczyk+ 2020'
    },

    # GRAVITATIONAL LENSING TIME DELAYS
    # These give time-delay distances, related to D_A
    {
        'name': 'B1608+656',
        'z_lens': 0.6304,
        'z_source': 1.394,
        'D_dt': 5156,  # Time-delay distance Mpc
        'D_dt_err': 296,
        'method': 'Lensing Time Delay',
        'ref': 'H0LiCOW'
    },
    {
        'name': 'HE0435-1223',
        'z_lens': 0.4546,
        'z_source': 1.693,
        'D_dt': 2612,  # Time-delay distance Mpc
        'D_dt_err': 200,
        'method': 'Lensing Time Delay',
        'ref': 'H0LiCOW Wong+ 2017'
    },
]

# =============================================================
# ANALYSIS
# =============================================================
print("=" * 80)
print("GEOMETRIC vs REDSHIFT-DERIVED DISTANCE COMPARISON")
print("=" * 80)
print()

print("TRIPHASE DRAG FORMULA CHECK:")
print(f"  At z=1100 (CMB), n=5: Drag = {drag_triphase(1100, 5):.6f} = {drag_triphase(1100, 5)*100:.3f}%")
print(f"  At z=1 (moderate): Drag = {drag_triphase(1, 5):.6e} = {drag_triphase(1, 5)*100:.6f}%")
print(f"  At z=0.01 (local): Drag = {drag_triphase(0.01, 5):.6e} = {drag_triphase(0.01, 5)*100:.8f}%")
print()

print("-" * 80)
print(f"{'Object':<15} {'z':>8} {'D_geo':>10} {'D_Planck':>10} {'D_Local':>10} {'Geo/Plk':>8} {'Geo/Loc':>8}")
print(f"{'':15} {'':>8} {'(Mpc)':>10} {'(Mpc)':>10} {'(Mpc)':>10} {'ratio':>8} {'ratio':>8}")
print("-" * 80)

results = []

for obj in GEOMETRIC_DATA:
    if 'z' in obj:
        z = obj['z']
        D_geo = obj['D_geo']

        # Standard cosmology predictions
        if z < 0.01:
            # Very low z - use simple Hubble law
            D_planck = D_hubble(z, H0_planck)
            D_local = D_hubble(z, H0_local)
        else:
            # Use full ΛCDM
            try:
                D_planck = D_angular_diameter(z, H0_planck)
                D_local = D_angular_diameter(z, H0_local)
            except:
                D_planck = D_hubble(z, H0_planck)
                D_local = D_hubble(z, H0_local)

        ratio_planck = D_geo / D_planck if D_planck > 0 else np.nan
        ratio_local = D_geo / D_local if D_local > 0 else np.nan

        print(f"{obj['name']:<15} {z:>8.5f} {D_geo:>10.3f} {D_planck:>10.3f} {D_local:>10.3f} {ratio_planck:>8.3f} {ratio_local:>8.3f}")

        results.append({
            'name': obj['name'],
            'z': z,
            'D_geo': D_geo,
            'D_planck': D_planck,
            'D_local': D_local,
            'ratio_planck': ratio_planck,
            'ratio_local': ratio_local,
            'method': obj['method']
        })

    elif 'z_lens' in obj:
        # Lensing - more complex, skip for now
        print(f"{obj['name']:<15} z_s={obj['z_source']:.3f} -- Lensing (time-delay distance) --")

print("-" * 80)
print()

# =============================================================
# STATISTICAL ANALYSIS
# =============================================================
if results:
    ratios_planck = [r['ratio_planck'] for r in results if not np.isnan(r['ratio_planck'])]
    ratios_local = [r['ratio_local'] for r in results if not np.isnan(r['ratio_local'])]

    print("STATISTICAL SUMMARY:")
    print(f"  Using Planck H0 = {H0_planck} km/s/Mpc:")
    print(f"    Mean D_geo/D_Planck = {np.mean(ratios_planck):.4f}")
    print(f"    Std dev = {np.std(ratios_planck):.4f}")
    print()
    print(f"  Using Local H0 = {H0_local} km/s/Mpc:")
    print(f"    Mean D_geo/D_Local = {np.mean(ratios_local):.4f}")
    print(f"    Std dev = {np.std(ratios_local):.4f}")
    print()

# =============================================================
# TRIPHASE INTERPRETATION
# =============================================================
print("=" * 80)
print("TRIPHASE INTERPRETATION")
print("=" * 80)
print()
print("The Hubble tension: Local measurements give H0 ~ 73, CMB gives H0 ~ 67")
print("Ratio: 73/67 = 1.09 -- a 9% discrepancy")
print()
print("If TriPhase is correct:")
print("  - Light slows down as it travels through vacuum")
print("  - Observed redshift includes contribution from energy loss (drag)")
print("  - So for given z, actual expansion-based distance is DIFFERENT")
print()
print("GEOMETRIC distances don't assume H0 - they measure D directly!")
print("If geometric D consistently differs from redshift D, that's evidence.")
print()

# =============================================================
# KEY TEST: Does H0 appear to change with distance?
# =============================================================
print("=" * 80)
print("KEY TEST: Implied H0 vs Distance")
print("=" * 80)
print()
print("If light slows, the 'implied H0' should vary with distance:")
print(f"{'Object':<15} {'z':>8} {'D_geo':>10} {'H0_implied':>12}")
print("-" * 50)

for r in results:
    # H0 = cz/D (for low z)
    H0_implied = c * r['z'] / r['D_geo']
    print(f"{r['name']:<15} {r['z']:>8.5f} {r['D_geo']:>10.3f} {H0_implied:>12.1f}")

print("-" * 50)
print()
print("If H0_implied increases with distance, that's consistent with TriPhase")
print("(more distant objects have accumulated more drag, so z is 'too high'")
print(" relative to D, making H0_implied appear larger)")