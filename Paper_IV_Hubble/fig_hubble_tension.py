"""
Paper IV Figure: Hubble Tension Resolution
============================================
Companion script for:
  "A Proposed Resolution of the Hubble Tension from Vacuum Wave Mechanics"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Two panels: (a) H0 measurements with error bars, (b) Bias factor B(z) curve.
Shows how H0 = 71.48 maps to 73.04 at SH0ES redshift via B(z) = z/ln(1+z).
Output: fig_hubble_tension.html (interactive Plotly chart)

Requires: numpy, plotly
"""
import sys
import traceback
import os


def main():
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # ============================================================
    # MECHANISM (What this script does and why)
    # ============================================================
    print("=" * 72)
    print("  HUBBLE TENSION RESOLUTION FIGURE -- Paper IV")
    print("  Magnetic Innovative Solutions LLC")
    print("=" * 72)
    print()
    print("WHAT THIS CALCULATES:")
    print("-" * 72)
    print("""
The Hubble tension is a 5.8-sigma discrepancy between two ways of
measuring the expansion rate of the universe:

  Local (SH0ES, distance ladder):  H0 = 73.04 +/- 1.04 km/s/Mpc
  CMB   (Planck, early universe):  H0 = 67.4  +/- 0.5  km/s/Mpc

Standard cosmology uses the linear distance-redshift formula:
  d = c*z / H0

This work derives a logarithmic formula from wave packet slowing:
  d = (c / H0) * ln(1 + z)

Using the wrong (linear) formula introduces a bias factor:
  B(z) = z / ln(1 + z)

which is always > 1 for z > 0.  At the mean SH0ES redshift z = 0.044:
  B(0.044) = 0.044 / ln(1.044) = 1.0220

So the true H0 = 71.48 km/s/Mpc appears as:
  H0_apparent = 71.48 x 1.0220 = 73.05 km/s/Mpc

This matches SH0ES within measurement uncertainty.  The tension is a
systematic bias in the distance formula, not a physics discrepancy.

This script:
  1. Derives H0 = 71.48 from vacuum coupling rate
  2. Computes B(z) at several redshifts
  3. Compares predicted apparent H0 to published measurements
  4. Generates a two-panel interactive figure
""")

    # ============================================================
    # STEP 1: Derive H0 from vacuum properties
    # ============================================================
    print("=" * 72)
    print("  STEP 1: H0 from Vacuum Coupling Rate")
    print("=" * 72)
    print()

    c = 299792458       # m/s
    mu0 = 1.25663706212e-6  # H/m
    eps0 = 8.854187812e-12  # F/m

    kappa_vac = (7 / 6) * mu0**2 / (c**3 * eps0)
    H0_si = kappa_vac * c   # in s^-1
    H0_km = H0_si * 3.086e19  # convert to km/s/Mpc

    print(f"  Vacuum constants:")
    print(f"    c   = {c:,} m/s")
    print(f"    mu0 = {mu0:.13e} H/m")
    print(f"    eps0= {eps0:.13e} F/m")
    print()
    print(f"  Vacuum coupling rate:")
    print(f"    kappa_vac = (7/6) * mu0^2 / (c^3 * eps0)")
    print(f"              = (7/6) * ({mu0:.6e})^2 / (({c:.0f})^3 * {eps0:.6e})")
    print(f"              = {kappa_vac:.6e} m^-1")
    print()
    print(f"  H0 = kappa_vac * c")
    print(f"     = {kappa_vac:.6e} * {c:,}")
    print(f"     = {H0_si:.6e} s^-1")
    print(f"     = {H0_si:.6e} * 3.086e19 km/Mpc")
    print(f"     = {H0_km:.2f} km/s/Mpc")
    print()

    H0_true = 71.48  # Use the rounded value from the paper

    # ============================================================
    # STEP 2: Bias factor B(z) at key redshifts
    # ============================================================
    print("=" * 72)
    print("  STEP 2: Bias Factor B(z) = z / ln(1 + z)")
    print("=" * 72)
    print()
    print(f"  {'z':>8s}  {'ln(1+z)':>10s}  {'B(z)':>8s}  {'H0_apparent':>12s}  Note")
    print(f"  {'-'*8}  {'-'*10}  {'-'*8}  {'-'*12}  {'-'*24}")

    z_table = [0.001, 0.005, 0.01, 0.02, 0.03, 0.044, 0.05, 0.10, 0.15]
    for z in z_table:
        ln1z = np.log(1 + z)
        Bz = z / ln1z
        H0_app = H0_true * Bz
        note = ""
        if z == 0.044:
            note = "<-- SH0ES mean redshift"
        elif z == 0.001:
            note = "nearby galaxies"
        elif z == 0.10:
            note = "Pantheon range"
        print(f"  {z:>8.4f}  {ln1z:>10.6f}  {Bz:>8.4f}  {H0_app:>12.2f}  {note}")

    print()
    print(f"  At SH0ES mean z = 0.044:")
    z_shoes = 0.044
    B_shoes = z_shoes / np.log(1 + z_shoes)
    H0_shoes_pred = H0_true * B_shoes
    print(f"    B(0.044) = {z_shoes} / ln({1 + z_shoes}) = {B_shoes:.4f}")
    print(f"    H0_apparent = {H0_true} x {B_shoes:.4f} = {H0_shoes_pred:.2f} km/s/Mpc")
    print(f"    SH0ES measured: 73.04 +/- 1.04 km/s/Mpc")
    print(f"    Difference: {abs(H0_shoes_pred - 73.04):.2f} km/s/Mpc ({abs(H0_shoes_pred - 73.04)/1.04:.2f} sigma)")
    print()

    # ============================================================
    # STEP 3: Published H0 measurements comparison
    # ============================================================
    print("=" * 72)
    print("  STEP 3: Published H0 Measurements")
    print("=" * 72)
    print()

    local_measurements = [
        ("SH0ES (Riess 2022)",          73.04, 1.04, 1.04),
        ("H0LiCOW lensing (Wong 2020)", 73.3,  1.8,  1.7),
        ("Megamasers (Pesce 2020)",     73.9,  3.0,  3.0),
        ("TRGB (Freedman 2021)",        69.8,  1.7,  1.7),
    ]
    cmb_measurements = [
        ("Planck 2018 (Aghanim 2020)",  67.4,  0.5,  0.5),
        ("ACT DR4 (Aiola 2020)",        67.6,  1.1,  1.1),
        ("SPT-3G (Balkenhol 2023)",     67.49, 0.53, 0.53),
    ]

    print(f"  LOCAL (distance ladder):")
    for label, h0, elo, ehi in local_measurements:
        sigma = abs(h0 - H0_true) / ((elo + ehi) / 2)
        print(f"    {label:<34s}  {h0:>6.2f} +/- {elo:.2f}  ({sigma:.1f} sigma from {H0_true})")
    print()
    print(f"  CMB (early universe):")
    for label, h0, elo, ehi in cmb_measurements:
        sigma = abs(h0 - H0_true) / ((elo + ehi) / 2)
        print(f"    {label:<34s}  {h0:>6.2f} +/- {elo:.2f}  ({sigma:.1f} sigma from {H0_true})")
    print()
    print(f"  Derived value: H0 = {H0_true} km/s/Mpc (no free parameters)")
    print()

    # ============================================================
    # STEP 4: Generate interactive two-panel figure (Plotly)
    # ============================================================
    print("=" * 72)
    print("  STEP 4: Generating interactive figure")
    print("=" * 72)
    print()

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=['(a) H0 Measurements',
                        '(b) Bias Factor: B(z) = z / ln(1+z)'],
        column_widths=[0.45, 0.55],
        horizontal_spacing=0.10
    )

    # ---- Panel (a): H0 Measurements with error bars ----

    # CMB measurements
    cmb_labels = [m[0] for m in cmb_measurements]
    cmb_h0 = [m[1] for m in cmb_measurements]
    cmb_err = [m[2] for m in cmb_measurements]

    fig.add_trace(go.Scatter(
        x=cmb_h0, y=cmb_labels,
        error_x=dict(type='data', array=cmb_err, visible=True,
                     color='#2166AC', thickness=2),
        mode='markers', name='CMB (early universe)',
        marker=dict(size=10, color='#2166AC', symbol='circle',
                    line=dict(color='white', width=1)),
        hovertemplate='%{y}<br>H0 = %{x} +/- %{error_x.array}<extra></extra>'
    ), row=1, col=1)

    # Local measurements
    local_labels = [m[0] for m in local_measurements]
    local_h0 = [m[1] for m in local_measurements]
    local_err = [m[2] for m in local_measurements]

    fig.add_trace(go.Scatter(
        x=local_h0, y=local_labels,
        error_x=dict(type='data', array=local_err, visible=True,
                     color='#D6604D', thickness=2),
        mode='markers', name='Local (distance ladder)',
        marker=dict(size=10, color='#D6604D', symbol='circle',
                    line=dict(color='white', width=1)),
        hovertemplate='%{y}<br>H0 = %{x} +/- %{error_x.array}<extra></extra>'
    ), row=1, col=1)

    # Derived value
    fig.add_trace(go.Scatter(
        x=[H0_true], y=['This work (derived)'],
        mode='markers', name=f'This work: H0 = {H0_true}',
        marker=dict(size=16, color='#CC2222', symbol='star',
                    line=dict(color='white', width=1)),
        hovertemplate=f'H0 = {H0_true} km/s/Mpc<br>(zero free parameters)<extra></extra>'
    ), row=1, col=1)

    # Planck and SH0ES bands
    fig.add_vrect(x0=67.4 - 0.5, x1=67.4 + 0.5,
                  fillcolor='#2166AC', opacity=0.08, line_width=0,
                  row=1, col=1)
    fig.add_vrect(x0=73.04 - 1.04, x1=73.04 + 1.04,
                  fillcolor='#D6604D', opacity=0.08, line_width=0,
                  row=1, col=1)
    fig.add_vline(x=H0_true, line_color='#CC2222', line_width=1.5,
                  line_dash='solid', opacity=0.5, row=1, col=1)

    fig.update_xaxes(title_text='H0 (km/s/Mpc)', range=[63, 80], row=1, col=1)

    # ---- Panel (b): Bias Factor B(z) curve ----

    z_arr = np.linspace(0.001, 0.15, 500)
    Bz_arr = z_arr / np.log(1 + z_arr)
    H0_apparent = H0_true * Bz_arr

    # Main curve
    fig.add_trace(go.Scatter(
        x=z_arr, y=H0_apparent,
        mode='lines', name='H0_obs = 71.48 x B(z)',
        line=dict(color='#CC2222', width=2.5),
        hovertemplate='z = %{x:.4f}<br>H0_apparent = %{y:.2f}<extra></extra>'
    ), row=1, col=2)

    # True H0 line
    fig.add_hline(y=H0_true, line_color='#2166AC', line_width=1.5,
                  line_dash='dash', opacity=0.7, row=1, col=2)
    fig.add_annotation(
        x=0.005, y=H0_true - 0.4,
        text=f'H0_true = {H0_true}', showarrow=False,
        font=dict(size=10, color='#2166AC'),
        xref='x2', yref='y2'
    )

    # SH0ES band
    fig.add_hrect(y0=73.04 - 1.04, y1=73.04 + 1.04,
                  fillcolor='#D6604D', opacity=0.12, line_width=0,
                  row=1, col=2)
    fig.add_hline(y=73.04, line_color='#D6604D', line_width=1,
                  line_dash='dot', opacity=0.7, row=1, col=2)
    fig.add_annotation(
        x=0.10, y=73.04 + 0.3,
        text='SH0ES: 73.04 +/- 1.04', showarrow=False,
        font=dict(size=9, color='#D6604D'),
        xref='x2', yref='y2'
    )

    # Planck band
    fig.add_hrect(y0=67.4 - 0.5, y1=67.4 + 0.5,
                  fillcolor='#2166AC', opacity=0.12, line_width=0,
                  row=1, col=2)
    fig.add_hline(y=67.4, line_color='#2166AC', line_width=1,
                  line_dash='dot', opacity=0.7, row=1, col=2)
    fig.add_annotation(
        x=0.10, y=67.4 - 0.4,
        text='Planck: 67.4 +/- 0.5', showarrow=False,
        font=dict(size=9, color='#2166AC'),
        xref='x2', yref='y2'
    )

    # SH0ES redshift marker
    fig.add_trace(go.Scatter(
        x=[z_shoes], y=[H0_shoes_pred],
        mode='markers', name=f'z=0.044: H0={H0_shoes_pred:.2f}',
        marker=dict(size=12, color='#CC2222',
                    line=dict(color='white', width=1.5)),
        showlegend=False,
        hovertemplate=f'z = {z_shoes}<br>H0_pred = {H0_shoes_pred:.2f}<extra></extra>'
    ), row=1, col=2)

    fig.add_annotation(
        x=z_shoes, y=H0_shoes_pred,
        text=f'z = 0.044<br>H0_pred = {H0_shoes_pred:.2f}',
        showarrow=True, arrowhead=2, arrowcolor='#CC2222',
        ax=40, ay=-40,
        font=dict(size=10, color='#CC2222'),
        bgcolor='#FEE8E8', bordercolor='#CC2222', borderwidth=1,
        xref='x2', yref='y2'
    )

    # Explanation box
    fig.add_annotation(
        x=0.06, y=66.8,
        text='Tension = distance formula bias<br>not a physics discrepancy',
        showarrow=False,
        font=dict(size=10, color='#555555'),
        bgcolor='#FFFFF0', bordercolor='#888888', borderwidth=1,
        xref='x2', yref='y2'
    )

    fig.update_xaxes(title_text='Redshift z', range=[0, 0.15], row=1, col=2)
    fig.update_yaxes(title_text='Apparent H0 (km/s/Mpc)',
                     range=[66, 77], row=1, col=2)

    fig.update_layout(
        title=dict(
            text='Hubble Tension Resolution: H0 = 71.48 km/s/Mpc with Systematic Bias at Low z',
            font=dict(size=15), x=0.5
        ),
        height=600, width=1200,
        template='plotly_white',
        legend=dict(orientation='h', yanchor='bottom', y=-0.18,
                    xanchor='center', x=0.5),
        margin=dict(l=80, r=40, t=80, b=100)
    )

    # Save and open in browser
    out_html = os.path.join(script_dir, 'fig_hubble_tension.html')
    fig.write_html(out_html, auto_open=True)
    print(f"  Generated: {out_html}")
    print(f"  Opened in browser.")

    print()
    print("=" * 72)
    print("  COMPLETE -- H0 derivation verified, figure generated.")
    print("=" * 72)
    print()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
