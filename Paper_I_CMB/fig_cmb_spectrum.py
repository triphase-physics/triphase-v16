"""
Paper I Figure: Planck CMB TT Power Spectrum with TriPhase Predictions
======================================================================
Companion figure for:
  "CMB Acoustic Peak Structure from Vacuum Electromagnetic Properties"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Generates a publication-quality figure showing:
  - Planck 2018 TT power spectrum (D_l vs multipole l)
  - Three-phase predicted peak positions (vertical lines)
  - Observed Gaussian-fitted peak positions with uncertainties
  - Zero free parameters: all predictions from vacuum constants

Output: fig_cmb_spectrum.pdf (for LaTeX inclusion)
        fig_cmb_spectrum.html (interactive browser version)

Data Source: Planck 2018 Release 3 (R3.01)
  https://irsa.ipac.caltech.edu/data/Planck/release_3/

Requires: numpy, matplotlib, plotly (optional, for HTML)
"""
import os
import sys
import traceback


def main():
    import numpy as np

    # ============================================================
    # Load Planck data
    # ============================================================
    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = script_dir

    data_file = "COM_PowerSpect_CMB-TT-full_R3.01.txt"
    subdirs = ["Planck_Data", "Notes", "", os.path.join("..", "..", "..", "Data", "CMB")]

    data_path = None
    for subdir in subdirs:
        candidate = os.path.join(script_dir, subdir, data_file)
        if os.path.exists(candidate):
            data_path = candidate
            break

    if data_path is None:
        print(f"ERROR: Cannot find {data_file}")
        print(f"Download from: https://irsa.ipac.caltech.edu/data/Planck/release_3/")
        return

    print(f"Loading: {os.path.basename(data_path)}")
    data = np.loadtxt(data_path, comments='#')
    ell = data[:, 0]
    Dl = data[:, 1]          # D_l = l(l+1)C_l / 2pi  [uK^2]
    err_minus = data[:, 2]
    err_plus = data[:, 3]

    # Symmetric error for fill_between
    err_sym = (err_minus + err_plus) / 2.0

    # Smooth with 21-point running average (Paper I method)
    kernel = 21
    Dl_smooth = np.convolve(Dl, np.ones(kernel) / kernel, mode='same')

    # Trim edges where convolution is unreliable
    trim = kernel // 2
    ell_plot = ell[trim:-trim]
    Dl_plot = Dl_smooth[trim:-trim]
    err_plot = err_sym[trim:-trim]
    Dl_raw = Dl[trim:-trim]

    # ============================================================
    # Paper I predictions and observations
    # ============================================================

    # Predictions from vacuum constants (zero free parameters)
    pred = {
        r'$\ell_1$': {'pos': 220.02, 'source': r'$\varepsilon_0$'},
        r'$\ell_2$': {'pos': 532.84, 'source': r'$\alpha, \mu_0$'},
        r'$\ell_3$': {'pos': 816.90, 'source': r'$Z_0$'},
    }

    # Planck 2018 observed (Gaussian-fitted)
    obs = {
        r'$\ell_1$': {'pos': 220.8, 'err': 3.5},
        r'$\ell_2$': {'pos': 533.2, 'err': 5.2},
        r'$\ell_3$': {'pos': 816.9, 'err': 2.8},
    }

    # ============================================================
    # Generate PDF figure with Matplotlib
    # ============================================================
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator

        print("Generating PDF figure...")

        fig, ax = plt.subplots(figsize=(10, 5.5))

        # Plot raw spectrum as faint points
        ax.scatter(ell_plot[::3], Dl_raw[::3], s=0.3, color='#cccccc',
                   alpha=0.5, zorder=1, rasterized=True)

        # Plot smoothed spectrum
        ax.plot(ell_plot, Dl_plot, color='#333333', linewidth=0.8,
                zorder=2, label='Planck 2018 TT (smoothed)')

        # Error band
        ax.fill_between(ell_plot, Dl_plot - err_plot, Dl_plot + err_plot,
                        color='#cccccc', alpha=0.3, zorder=1)

        # Predicted peak positions (solid red lines)
        colors_pred = ['#d62728', '#d62728', '#d62728']
        for i, (label, info) in enumerate(pred.items()):
            pos = info['pos']
            # Find D_l at this position for line height
            idx = np.argmin(np.abs(ell_plot - pos))
            y_top = Dl_plot[idx] * 1.12
            ax.axvline(pos, ymin=0, ymax=0.85, color='#d62728',
                       linewidth=1.5, linestyle='-', alpha=0.8, zorder=3)
            ax.annotate(
                f'{label} = {pos}\n({info["source"]})',
                xy=(pos, y_top), fontsize=9,
                ha='center', va='bottom', color='#d62728',
                fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='#d62728', alpha=0.9)
            )

        # Observed positions (blue dashed + error bars)
        for label, info in obs.items():
            pos = info['pos']
            err = info['err']
            idx = np.argmin(np.abs(ell_plot - pos))
            y_marker = Dl_plot[idx]
            ax.errorbar(pos, y_marker, xerr=err, fmt='o', color='#1f77b4',
                        markersize=5, capsize=4, capthick=1.2,
                        linewidth=1.2, zorder=5)

        # Axis labels and formatting
        ax.set_xlabel(r'Multipole $\ell$', fontsize=12)
        ax.set_ylabel(r'$\mathcal{D}_\ell = \ell(\ell+1)C_\ell / 2\pi$ [$\mu$K$^2$]',
                       fontsize=12)
        ax.set_title(
            r'Planck 2018 TT Power Spectrum with Three-Phase Predictions'
            '\n(zero free parameters)',
            fontsize=13, fontweight='bold')

        ax.set_xlim(2, 1200)
        ax.set_ylim(0, 7000)
        ax.xaxis.set_major_locator(MultipleLocator(200))
        ax.xaxis.set_minor_locator(MultipleLocator(50))

        # Legend
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        legend_elements = [
            Line2D([0], [0], color='#333333', linewidth=0.8,
                   label='Planck 2018 TT (21-pt smooth)'),
            Patch(facecolor='#cccccc', alpha=0.3, label='Measurement uncertainty'),
            Line2D([0], [0], color='#d62728', linewidth=1.5,
                   label='Predicted (vacuum constants)'),
            Line2D([0], [0], marker='o', color='#1f77b4', linewidth=0,
                   markersize=5, label='Observed (Gaussian fit)'),
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=8.5,
                  framealpha=0.9)

        plt.tight_layout()

        # Save PDF
        out_pdf = os.path.join(out_dir, "fig_cmb_spectrum.pdf")
        fig.savefig(out_pdf, format='pdf', dpi=300, bbox_inches='tight')
        print(f"Saved: {out_pdf}")

        # Also save PNG for quick preview
        out_png = os.path.join(out_dir, "fig_cmb_spectrum.png")
        fig.savefig(out_png, format='png', dpi=200, bbox_inches='tight')
        print(f"Saved: {out_png}")

        plt.close(fig)

    except ImportError:
        print("matplotlib not available -- skipping PDF generation")

    # ============================================================
    # Generate interactive HTML with Plotly
    # ============================================================
    try:
        import plotly.graph_objects as go

        print("Generating HTML figure...")

        fig2 = go.Figure()

        # Raw data as faint scatter
        fig2.add_trace(go.Scatter(
            x=ell_plot[::2], y=Dl_raw[::2],
            mode='markers',
            marker=dict(size=1.5, color='#cccccc', opacity=0.4),
            name='Planck 2018 TT (raw)',
            hoverinfo='skip'
        ))

        # Smoothed spectrum
        fig2.add_trace(go.Scatter(
            x=ell_plot, y=Dl_plot,
            mode='lines',
            line=dict(color='#333333', width=1.5),
            name='Planck 2018 TT (smoothed)'
        ))

        # Error band (upper)
        fig2.add_trace(go.Scatter(
            x=np.concatenate([ell_plot, ell_plot[::-1]]),
            y=np.concatenate([Dl_plot + err_plot, (Dl_plot - err_plot)[::-1]]),
            fill='toself',
            fillcolor='rgba(200,200,200,0.25)',
            line=dict(color='rgba(200,200,200,0)'),
            hoverinfo='skip',
            name='Uncertainty band'
        ))

        # Predicted peaks
        peak_labels = [r'ℓ₁', r'ℓ₂', r'ℓ₃']
        sources = ['ε₀', 'α, μ₀', 'Z₀']
        pred_pos = [220.02, 532.84, 816.90]
        obs_pos = [220.8, 533.2, 816.9]
        obs_err = [3.5, 5.2, 2.8]

        for i, (plabel, src, pp, op, oe) in enumerate(
                zip(peak_labels, sources, pred_pos, obs_pos, obs_err)):
            # Predicted vertical line
            fig2.add_vline(
                x=pp, line_width=2, line_dash="solid",
                line_color="#d62728", opacity=0.7,
                annotation_text=f"{plabel} = {pp} ({src})",
                annotation_position="top",
                annotation_font_size=11,
                annotation_font_color="#d62728"
            )

            # Observed marker
            idx = np.argmin(np.abs(ell_plot - op))
            fig2.add_trace(go.Scatter(
                x=[op], y=[Dl_plot[idx]],
                mode='markers',
                marker=dict(size=10, color='#1f77b4', symbol='circle'),
                error_x=dict(type='constant', value=oe, color='#1f77b4',
                             thickness=2, width=6),
                name=f'{plabel} obs = {op} ± {oe}',
                showlegend=(i == 0),
            ))

        fig2.update_layout(
            title=dict(
                text='Planck 2018 TT Power Spectrum with Three-Phase Predictions<br>'
                     '<sub>(zero free parameters — all positions from vacuum constants)</sub>',
                font=dict(size=16)
            ),
            xaxis_title='Multipole ℓ',
            yaxis_title='D_ℓ = ℓ(ℓ+1)C_ℓ / 2π [μK²]',
            xaxis=dict(range=[2, 1200]),
            yaxis=dict(range=[0, 7000]),
            template='plotly_white',
            width=1100,
            height=600,
            legend=dict(x=0.65, y=0.97, font=dict(size=10)),
        )

        out_html = os.path.join(out_dir, "fig_cmb_spectrum.html")
        fig2.write_html(out_html, auto_open=True)
        print(f"Saved: {out_html}")

    except ImportError:
        print("plotly not available -- skipping HTML generation")

    print()
    print("Done.")


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
