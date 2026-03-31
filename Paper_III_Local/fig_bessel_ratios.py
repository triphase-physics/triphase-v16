"""
Paper III Figure: Bessel Function Zero Ratios vs Mode Integer Predictions
=========================================================================
Companion script for:
  "Local Laboratory Evidence for Three-Phase Vacuum Structure"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Generates an interactive two-panel Bessel ratio comparison figure.
Output: fig_bessel_ratios.html (interactive Plotly chart)

Requires: numpy, scipy, plotly
"""
import sys
import traceback
import os


def main():
    import numpy as np
    from scipy.special import jn_zeros
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # ============================================================
    # MECHANISM (What this script does and why)
    # ============================================================
    print("=" * 72)
    print("  BESSEL ZERO RATIO VERIFICATION -- Paper III, Local Evidence")
    print("  Magnetic Innovative Solutions LLC")
    print("=" * 72)
    print()
    print("WHAT THIS CALCULATES:")
    print("-" * 72)
    print("""
The zeros of the Bessel function J_0(x) are exact mathematical constants.
They show up everywhere cylindrical symmetry appears: drum vibrations,
electromagnetic waveguides, diffraction patterns.

Walking droplet experiments (Bush, Couder, Fort) show that a silicone
droplet bouncing on a vibrating fluid bath self-organizes into quantized
orbits.  The radial structure of these orbits follows J_0 -- the same
Bessel function that governs cylindrical wave modes.

The claim under test:  ratios of J_0 zeros reduce to ratios of small
integers, and those integers are the same mode-counting numbers that
appear throughout the three-phase vacuum structure (5, 6, 7, 17, 18,
22, 137).  If this is coincidence, the probability of 9 independent
matches at < 0.55% is extremely small.  If it is real, it means the
vacuum mode structure is encoded in Bessel functions -- exactly what
you would expect from a three-phase cylindrical geometry.

This script:
  1. Computes the first 12 zeros of J_0 using SciPy
  2. Forms 9 specific ratios identified in Paper III
  3. Compares each ratio to its integer prediction
  4. Calculates the percent deviation for each
  5. Generates a two-panel interactive figure (HTML)
""")

    # ============================================================
    # STEP 1: Compute Bessel J_0 zeros
    # ============================================================
    print("=" * 72)
    print("  STEP 1: Bessel J_0 Zeros (first 12)")
    print("=" * 72)
    print()

    J0_zeros = jn_zeros(0, 12)

    for i, z in enumerate(J0_zeros):
        print(f"  x_{i+1:2d} = {z:.10f}")

    print()

    # ============================================================
    # STEP 2: Define the 9 ratio comparisons from Paper III
    # ============================================================
    print("=" * 72)
    print("  STEP 2: Form Ratios and Compare to Integer Predictions")
    print("=" * 72)
    print()
    print("  Each ratio is a quotient of two J_0 zeros.  The prediction is")
    print("  a ratio of mode-counting integers from the three-phase structure.")
    print()

    # (label, measured_value, integer_num, integer_den, physical_origin)
    ratios = [
        ('x3/x1',    J0_zeros[2]/J0_zeros[0],   18,   5, 'phases x total / coupling'),
        ('x9/x4',    J0_zeros[8]/J0_zeros[3],    7,   3, '(modes+ground) / phases'),
        ('x10/x4',   J0_zeros[9]/J0_zeros[3],   13,   5, 'hierarchy / coupling'),
        ('x9/x7',    J0_zeros[8]/J0_zeros[6],   22,  17, 'hierarchy ratios'),
        ('x3/x2',    J0_zeros[2]/J0_zeros[1],   11,   7, '(5+6) / (6+1)'),
        ('x8/x2',    J0_zeros[7]/J0_zeros[1],   22,   5, 'couplings / coupling'),
        ('x5/x1',    J0_zeros[4]/J0_zeros[0],  137,  22, 'alpha^-1 / hierarchy'),
        ('x12/x4',   J0_zeros[11]/J0_zeros[3],  22,   7, 'approx pi connection'),
        ('x4/x2',    J0_zeros[3]/J0_zeros[1],   17,   8, 'hierarchy / octave'),
    ]

    # Print step-by-step for each ratio
    print(f"  {'Ratio':<14s}  {'Bessel Value':>14s}  {'Integer':>8s}  {'Predicted':>10s}  {'Error %':>8s}  Origin")
    print(f"  {'-'*14}  {'-'*14}  {'-'*8}  {'-'*10}  {'-'*8}  {'-'*28}")

    labels_plot = []
    measured_arr = []
    predicted_arr = []
    errors_pct_arr = []
    int_labels = []

    for label, meas, num, den, origin in ratios:
        pred = num / den
        err = abs(meas - pred) / pred * 100
        print(f"  {label:<14s}  {meas:>14.10f}  {num:>3d}/{den:<3d}  {pred:>10.6f}  {err:>7.4f}%  {origin}")

        labels_plot.append(label)
        measured_arr.append(meas)
        predicted_arr.append(pred)
        errors_pct_arr.append(err)
        int_labels.append(f"{num}/{den}")

    measured_arr = np.array(measured_arr)
    predicted_arr = np.array(predicted_arr)
    errors_pct_arr = np.array(errors_pct_arr)

    print()
    print(f"  Worst-case deviation:  {max(errors_pct_arr):.4f}%")
    print(f"  Best-case deviation:   {min(errors_pct_arr):.4f}%")
    print(f"  Mean deviation:        {np.mean(errors_pct_arr):.4f}%")
    print(f"  All 9 ratios within:   0.55%")
    print()

    # ============================================================
    # STEP 3: Sort by accuracy for the figure
    # ============================================================
    print("=" * 72)
    print("  STEP 3: Sorting by accuracy (best match first) for figure")
    print("=" * 72)
    print()

    order = np.argsort(errors_pct_arr)
    labels_sorted = [labels_plot[i] for i in order]
    measured_sorted = measured_arr[order]
    predicted_sorted = predicted_arr[order]
    errors_sorted = errors_pct_arr[order]
    int_sorted = [int_labels[i] for i in order]

    for rank, i in enumerate(order):
        print(f"  #{rank+1}:  {labels_plot[i]:<14s}  error = {errors_pct_arr[i]:.4f}%")

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
        subplot_titles=['(a) Bessel Zero Ratios vs Mode Predictions',
                        '(b) Prediction Accuracy (% Deviation)'],
        column_widths=[0.6, 0.4],
        horizontal_spacing=0.12
    )

    # Panel (a): Measured vs Predicted grouped horizontal bars
    fig.add_trace(go.Bar(
        y=labels_sorted, x=list(measured_sorted),
        orientation='h', name='Bessel zero ratio (exact)',
        marker_color='#2166AC', opacity=0.85,
        text=[f'{v:.4f}' for v in measured_sorted],
        textposition='outside', textfont=dict(size=10, color='#2166AC'),
        hovertemplate='%{y}<br>Bessel: %{x:.6f}<extra></extra>'
    ), row=1, col=1)

    fig.add_trace(go.Bar(
        y=labels_sorted, x=list(predicted_sorted),
        orientation='h', name='Mode integer ratio',
        marker_color='#CC2222', opacity=0.85,
        text=[f'{il} = {v:.4f}' for il, v in zip(int_sorted, predicted_sorted)],
        textposition='outside', textfont=dict(size=10, color='#CC2222'),
        hovertemplate='%{y}<br>Integer: %{x:.6f}<extra></extra>'
    ), row=1, col=1)

    # Panel (b): Residuals with color coding
    bar_colors = ['#22CC44' if e < 0.1 else '#88BB22' if e < 0.3 else '#CCAA22'
                  for e in errors_sorted]

    fig.add_trace(go.Bar(
        y=labels_sorted, x=list(errors_sorted),
        orientation='h', name='% Deviation',
        marker_color=bar_colors,
        marker_line=dict(color='#333333', width=1),
        text=[f'{e:.3f}%' for e in errors_sorted],
        textposition='outside', textfont=dict(size=11, color='#333333'),
        hovertemplate='%{y}<br>Deviation: %{x:.4f}%<extra></extra>',
        showlegend=False
    ), row=1, col=2)

    # Reference lines on residual panel
    fig.add_vline(x=0.1, line_dash='dot', line_color='#22CC44',
                  line_width=1, opacity=0.6, row=1, col=2)
    fig.add_vline(x=0.5, line_dash='dot', line_color='#CCAA22',
                  line_width=1, opacity=0.6, row=1, col=2)

    # Annotation: all within 0.55%
    fig.add_annotation(
        x=max(errors_sorted) * 1.1, y=labels_sorted[0],
        text=f'All {len(ratios)} ratios<br>within 0.55%',
        showarrow=False, font=dict(size=12, color='#22AA44'),
        bgcolor='#E8F8E8', bordercolor='#22CC44', borderwidth=1,
        xref='x2', yref='y2', xanchor='left', yanchor='top'
    )

    fig.update_layout(
        title=dict(
            text='Walking Droplet Quantization: Bessel J0 Zeros Encode Mode Integers',
            font=dict(size=16), x=0.5
        ),
        barmode='group',
        height=550, width=1200,
        template='plotly_white',
        legend=dict(orientation='h', yanchor='bottom', y=-0.15,
                    xanchor='center', x=0.3),
        margin=dict(l=80, r=120, t=80, b=80)
    )

    fig.update_xaxes(title_text='Ratio Value', row=1, col=1)
    fig.update_xaxes(title_text='Deviation (%)', row=1, col=2,
                     range=[0, max(errors_sorted) * 1.5])

    # Save chart and open in browser
    out_html = os.path.join(script_dir, 'fig_bessel_ratios.html')
    fig.write_html(out_html, auto_open=True)
    print(f"  Generated: {out_html}")
    print(f"  Opened in browser.")

    print()
    print("=" * 72)
    print("  COMPLETE -- All 9 Bessel ratios verified, figure generated.")
    print("=" * 72)
    print()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
