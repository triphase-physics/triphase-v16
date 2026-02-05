"""
Paper I Calculation: [111] Crystallographic Projection (90 -> 120 deg)
=====================================================================
Companion script for:
  "CMB Acoustic Peak Structure from Vacuum Electromagnetic Properties"
  C.R. Fuccillo, Magnetic Innovative Solutions LLC

Verifies that three perpendicular axes viewed along the [111] body
diagonal project to 120-degree separation -- the geometric foundation
of three-phase structure.

Requires: numpy
"""
import sys
import traceback


def main():
    import numpy as np

    print("=" * 72)
    print("  [111] CRYSTALLOGRAPHIC PROJECTION -- Paper I")
    print("  Magnetic Innovative Solutions LLC")
    print("=" * 72)
    print()
    print("WHAT THIS CALCULATES:")
    print("-" * 72)
    print("""
    The three vacuum electromagnetic properties (epsilon_0, mu_0, Z_0)
    are functionally independent -- they live along orthogonal axes in
    a 3D parameter space (like x, y, z axes of a cube).

    When we observe physics, we project from 3D into 2D (observable
    space).  The natural projection direction is the body diagonal of
    the cube -- the [111] direction in crystallographic notation.

    Under [111] projection, three perpendicular 90-degree axes appear
    at 120-degree separation in the projected plane.  This is the same
    geometry as a three-phase electrical system (three voltages at 120
    degrees, like a Mercedes-Benz logo).

    This is pure geometry -- no physics assumptions.  If three
    independent quantities project to 2D, they MUST form 120-degree
    separation.  The script verifies this numerically.

    The physical consequence: vacuum properties that are orthogonal in
    parameter space produce three-phase oscillatory structure when
    projected into observables.
    """)

    # ============================================================
    # STEP 1: Define the orthogonal axes
    # ============================================================
    print("=" * 72)
    print("  STEP 1: Orthogonal Unit Vectors")
    print("=" * 72)
    print()

    x_hat = np.array([1, 0, 0])
    y_hat = np.array([0, 1, 0])
    z_hat = np.array([0, 0, 1])

    print(f"  Three perpendicular unit vectors (the cube edges):")
    print(f"    x_hat = [{x_hat[0]}, {x_hat[1]}, {x_hat[2]}]  (epsilon_0 direction)")
    print(f"    y_hat = [{y_hat[0]}, {y_hat[1]}, {y_hat[2]}]  (mu_0 direction)")
    print(f"    z_hat = [{z_hat[0]}, {z_hat[1]}, {z_hat[2]}]  (Z_0 direction)")
    print()
    print(f"  All pairwise angles = 90 degrees (by construction)")
    print()

    # ============================================================
    # STEP 2: Define the [111] viewing direction
    # ============================================================
    print("=" * 72)
    print("  STEP 2: [111] Viewing Direction")
    print("=" * 72)
    print()

    n_raw = np.array([1, 1, 1])
    n_hat = n_raw / np.linalg.norm(n_raw)

    print(f"  Body diagonal of the cube:")
    print(f"    [111] = [{n_raw[0]}, {n_raw[1]}, {n_raw[2]}]")
    print(f"    |[111]| = sqrt(1^2 + 1^2 + 1^2) = sqrt(3) = {np.sqrt(3):.10f}")
    print(f"    n_hat = [111] / sqrt(3) = [{n_hat[0]:.10f}, {n_hat[1]:.10f}, {n_hat[2]:.10f}]")
    print()
    print(f"  This is the direction you look from to see all three axes equally.")
    print()

    # ============================================================
    # STEP 3: Project each axis onto the plane perpendicular to [111]
    # ============================================================
    print("=" * 72)
    print("  STEP 3: Project Axes onto Plane Perpendicular to [111]")
    print("=" * 72)
    print()
    print(f"  Projection formula: v_proj = v - (v . n_hat) * n_hat")
    print(f"  This removes the component along the viewing direction.")
    print()


    def project(v, n):
        return v - np.dot(v, n) * n


    x_proj = project(x_hat, n_hat)
    y_proj = project(y_hat, n_hat)
    z_proj = project(z_hat, n_hat)

    # Show the dot products
    print(f"  x_hat . n_hat = {np.dot(x_hat, n_hat):.10f}")
    print(f"  y_hat . n_hat = {np.dot(y_hat, n_hat):.10f}")
    print(f"  z_hat . n_hat = {np.dot(z_hat, n_hat):.10f}")
    print()
    print(f"  Projected vectors:")
    print(f"    x_proj = [{x_proj[0]:+.10f}, {x_proj[1]:+.10f}, {x_proj[2]:+.10f}]")
    print(f"    y_proj = [{y_proj[0]:+.10f}, {y_proj[1]:+.10f}, {y_proj[2]:+.10f}]")
    print(f"    z_proj = [{z_proj[0]:+.10f}, {z_proj[1]:+.10f}, {z_proj[2]:+.10f}]")
    print()
    print(f"  Projected lengths (all equal by symmetry):")
    print(f"    |x_proj| = {np.linalg.norm(x_proj):.10f}")
    print(f"    |y_proj| = {np.linalg.norm(y_proj):.10f}")
    print(f"    |z_proj| = {np.linalg.norm(z_proj):.10f}")
    print(f"    Expected: sqrt(2/3) = {np.sqrt(2/3):.10f}")
    print()

    # ============================================================
    # STEP 4: Measure angles between projected vectors
    # ============================================================
    print("=" * 72)
    print("  STEP 4: Angles Between Projected Vectors")
    print("=" * 72)
    print()


    def angle_between(v1, v2):
        cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return np.degrees(np.arccos(np.clip(cos_theta, -1, 1)))


    angle_xy = angle_between(x_proj, y_proj)
    angle_yz = angle_between(y_proj, z_proj)
    angle_zx = angle_between(z_proj, x_proj)

    print(f"  cos(theta) = (v1 . v2) / (|v1| * |v2|)")
    print()
    print(f"  x_proj . y_proj = {np.dot(x_proj, y_proj):.10f}")
    print(f"  Angle(x, y) = arccos({np.dot(x_proj, y_proj) / (np.linalg.norm(x_proj) * np.linalg.norm(y_proj)):.10f})")
    print(f"              = {angle_xy:.6f} degrees")
    print()
    print(f"  Angle(y, z) = {angle_yz:.6f} degrees")
    print(f"  Angle(z, x) = {angle_zx:.6f} degrees")
    print()
    print(f"  Sum of angles = {angle_xy + angle_yz + angle_zx:.6f} degrees (should be 360)")
    print()

    # ============================================================
    # STEP 5: Verify result
    # ============================================================
    print("=" * 72)
    print("  STEP 5: Verification")
    print("=" * 72)
    print()

    max_error = max(abs(angle_xy - 120), abs(angle_yz - 120), abs(angle_zx - 120))

    print(f"  All three angles:")
    print(f"    Angle(x, y) = {angle_xy:.6f} deg  (target: 120.000000)")
    print(f"    Angle(y, z) = {angle_yz:.6f} deg  (target: 120.000000)")
    print(f"    Angle(z, x) = {angle_zx:.6f} deg  (target: 120.000000)")
    print()
    print(f"  Maximum deviation from 120: {max_error:.2e} degrees")
    print()
    print(f"  Physical interpretation:")
    print(f"    90-degree orthogonal  ->  120-degree three-phase")
    print(f"    Cube parameter space  ->  Three-phase oscillation")
    print(f"    (epsilon_0, mu_0, Z_0)  ->  Three-phase vacuum structure")
    print()

    passed = max_error < 0.0001
    status = "PASSED" if passed else "FAILED"
    print(f"  Result: {status}")
    print()

    print("=" * 72)
    print(f"  COMPLETE -- [111] projection verified: 90 deg -> 120 deg")
    print("=" * 72)
    print()


if __name__ == '__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
    finally:
        input("Press Enter to exit...")
