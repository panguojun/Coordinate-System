"""
Examples for High-Precision Curvature Computation
=================================================

This file demonstrates how to use the new curvature computation module
added in v2.3.0.

**Authors:** Pan Guojun
Date: 2025-10-27
**DOI:** https://doi.org/10.5281/zenodo.14435613
"""

import math
from coordinate_system import (
    # Surface classes
    Sphere,
    Torus,

    # High-precision curvature functions
    CurvatureCalculator,
    gaussian_curvature,
    mean_curvature,
    principal_curvatures,
    all_curvatures,
)


def example1_simple_gaussian_curvature():
    """Example 1: Compute Gaussian curvature of a sphere"""
    print("="*60)
    print("Example 1: Gaussian Curvature of a Sphere")
    print("="*60)

    # Create a sphere with radius 2.0
    sphere = Sphere(radius=2.0)

    # Compute at point (θ=π/4, φ=π/6)
    u, v = math.pi/4, math.pi/6

    # Simple interface
    K = gaussian_curvature(sphere, u, v)

    print(f"\nSphere parameters:")
    print(f"  Radius R = 2.0")
    print(f"  Point (θ, φ) = ({u:.4f}, {v:.4f})")
    print(f"\nComputed Gaussian curvature:")
    print(f"  K = {K:.10f}")
    print(f"\nTheoretical value:")
    print(f"  K = 1/R² = {1.0/(2.0**2):.10f}")
    print(f"\nRelative error:")
    error_pct = abs(K - 0.25) / 0.25 * 100
    print(f"  {error_pct:.6f}%")
    print()


def example2_all_curvatures():
    """Example 2: Compute all curvature quantities"""
    print("="*60)
    print("Example 2: All Curvatures of a Torus")
    print("="*60)

    # Create a torus
    torus = Torus(major_radius=2.0, minor_radius=1.0)

    # Test three positions: outer, inner, top
    positions = [
        (math.pi/4, 0.0, "Outer side"),
        (math.pi/4, math.pi, "Inner side"),
        (math.pi/4, math.pi/2, "Top"),
    ]

    for u, v, label in positions:
        print(f"\n{label} (u={u:.4f}, v={v:.4f}):")
        print("-" * 50)

        # Compute all curvatures
        curv = all_curvatures(torus, u, v)

        print(f"  Gaussian curvature K   = {curv['K']:.8f}")
        print(f"  Mean curvature H       = {curv['H']:.8f}")
        print(f"  Principal curvature k1 = {curv['k1']:.8f}")
        print(f"  Principal curvature k2 = {curv['k2']:.8f}")

        # Verify relationships
        K_check = curv['k1'] * curv['k2']
        H_check = (curv['k1'] + curv['k2']) / 2

        print(f"\n  Verification:")
        print(f"    K = k1×k2 = {K_check:.8f} (✓)" if abs(K_check - curv['K']) < 1e-6 else f"    K = k1×k2 = {K_check:.8f} (✗)")
        print(f"    H = (k1+k2)/2 = {H_check:.8f} (✓)" if abs(H_check - curv['H']) < 1e-6 else f"    H = (k1+k2)/2 = {H_check:.8f} (✗)")

        # Classify geometry
        if abs(curv['K']) < 1e-6:
            geom_type = "Flat (zero curvature)"
        elif curv['K'] > 1e-6:
            geom_type = "Elliptic (positive curvature)"
        else:
            geom_type = "Hyperbolic (negative curvature)"

        print(f"  Geometry: {geom_type}")

    print()


def example3_curvature_calculator():
    """Example 3: Using CurvatureCalculator class"""
    print("="*60)
    print("Example 3: CurvatureCalculator Class Usage")
    print("="*60)

    # Create a sphere
    sphere = Sphere(radius=2.0)

    # Create calculator with specific step size
    calc = CurvatureCalculator(sphere, step_size=1e-4)

    # Compute point
    u, v = math.pi/4, math.pi/6

    print(f"\nComputing curvatures at (u={u:.4f}, v={v:.4f}):")
    print("-" * 50)

    # Method 1: Compute individually
    K = calc.compute_gaussian_curvature(u, v)
    H = calc.compute_mean_curvature(u, v)
    k1, k2, dir1, dir2 = calc.compute_principal_curvatures(u, v)

    print(f"\nIndividual computations:")
    print(f"  K  = {K:.10f}")
    print(f"  H  = {H:.10f}")
    print(f"  k1 = {k1:.10f}")
    print(f"  k2 = {k2:.10f}")

    print(f"\nPrincipal directions (unit vectors):")
    print(f"  dir1 = [{dir1[0]:.6f}, {dir1[1]:.6f}, {dir1[2]:.6f}]")
    print(f"  dir2 = [{dir2[0]:.6f}, {dir2[1]:.6f}, {dir2[2]:.6f}]")

    # Method 2: Compute all at once (more efficient)
    all_curv = calc.compute_all_curvatures(u, v)

    print(f"\nAll-at-once computation:")
    print(f"  K  = {all_curv['K']:.10f}")
    print(f"  H  = {all_curv['H']:.10f}")
    print(f"  First fundamental form det(g) = {all_curv['g'][0,0] * all_curv['g'][1,1] - all_curv['g'][0,1]**2:.6f}")
    print()


def example4_convergence_analysis():
    """Example 4: Convergence analysis"""
    print("="*60)
    print("Example 4: Convergence Analysis")
    print("="*60)

    # Create a sphere
    sphere = Sphere(radius=2.0)

    # Create calculator
    calc = CurvatureCalculator(sphere)

    # Run convergence test
    u, v = math.pi/4, math.pi/6
    results = calc.convergence_analysis(u, v)

    print(f"\nTesting convergence at (u={u:.4f}, v={v:.4f}):")
    print(f"Theoretical K = {0.25:.10f}")
    print("\n" + "-" * 70)
    print(f"{'Step Size':<12} {'K computed':<18} {'Error':<15} {'Convergence Ratio':<20}")
    print("-" * 70)

    K_true = 0.25
    errors = []

    for r in results:
        error = abs(r['K'] - K_true)
        errors.append(error)

        ratio = ""
        if len(errors) >= 2:
            # Convergence ratio (should be ~16 for O(h⁴))
            ratio_val = errors[-2] / errors[-1] if errors[-1] > 1e-15 else 0
            ratio = f"{ratio_val:.2f}"

        print(f"{r['h']:<12.1e} {r['K']:<18.10f} {error:<15.2e} {ratio:<20}")

    print("\nNote: Convergence ratio ~16 indicates O(h⁴) accuracy")
    print("      (halving h reduces error by factor of 16)")
    print()


def example5_custom_surface():
    """Example 5: Custom surface definition"""
    print("="*60)
    print("Example 5: Custom Surface - Hyperbolic Paraboloid")
    print("="*60)

    from coordinate_system import Surface, vec3

    class HyperbolicParaboloid(Surface):
        """Saddle surface: z = x²/a² - y²/b²"""

        def __init__(self, a=1.0, b=1.0):
            super().__init__()
            self.a = a
            self.b = b

        def position(self, u, v):
            """u = x, v = y"""
            x, y = u, v
            z = (x*x)/(self.a*self.a) - (y*y)/(self.b*self.b)
            return vec3(x, y, z)

    # Create custom surface
    saddle = HyperbolicParaboloid(a=1.0, b=1.0)

    # Compute curvatures
    u, v = 0.5, 0.3

    calc = CurvatureCalculator(saddle, step_size=1e-4)
    curv = calc.compute_all_curvatures(u, v)

    print(f"\nHyperbolic paraboloid at (x={u}, y={v}):")
    print("-" * 50)
    print(f"  Gaussian curvature K   = {curv['K']:.8f}")
    print(f"  Mean curvature H       = {curv['H']:.8f}")
    print(f"  Principal curvatures:")
    print(f"    k1 = {curv['k1']:.8f} (bending upward)")
    print(f"    k2 = {curv['k2']:.8f} (bending downward)")

    # Theoretical value
    denom = 1 + 4*u*u + 4*v*v
    K_theory = -4 / (denom * denom)

    print(f"\nTheoretical K = {K_theory:.8f}")
    print(f"Relative error = {abs(curv['K'] - K_theory)/abs(K_theory)*100:.4f}%")

    print(f"\nNote: K < 0 confirms this is a saddle surface (hyperbolic geometry)")
    print()


def main():
    """Run all examples"""
    print("\n" + "="*60)
    print(" High-Precision Curvature Computation Examples")
    print(" coordinate_system package v2.3.0")
    print("="*60 + "\n")

    example1_simple_gaussian_curvature()
    input("Press Enter to continue to Example 2...")

    example2_all_curvatures()
    input("Press Enter to continue to Example 3...")

    example3_curvature_calculator()
    input("Press Enter to continue to Example 4...")

    example4_convergence_analysis()
    input("Press Enter to continue to Example 5...")

    example5_custom_surface()

    print("="*60)
    print(" All examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
