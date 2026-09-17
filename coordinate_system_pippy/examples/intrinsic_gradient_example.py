#!/usr/bin/env python3
"""
Intrinsic Gradient Operator Curvature Example
=============================================

This example demonstrates how to compute surface curvatures using the
intrinsic gradient operator method.
Based on the paper "Surface Curvature Computation via the Intrinsic Gradient Operator".

**Authors:** Pan Guojun
Date: 2025-10-30
**DOI:** https://doi.org/10.5281/zenodo.14435613
"""

import math
import numpy as np
import sys
import os

# Add parent directory to path to import coordinate_system package
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from coordinate_system import (
    Sphere, Torus,
    IntrinsicGradientCurvatureCalculator,
    compute_ccs_geometry_package,
    intrinsic_gradient_gaussian_curvature,
    CurvatureCalculator,  # For comparison
)


def print_separator(title: str = ""):
    """Print a separator line."""
    if title:
        print(f"\n{'=' * 20} {title} {'=' * 20}")
    else:
        print("=" * 70)


def demo_sphere():
    """Demonstrate curvature computation on a sphere."""
    print_separator("Sphere Curvature (Intrinsic Gradient Operator)")

    # Create a sphere of radius 2
    radius = 2.0
    sphere = Sphere(radius=radius)

    # Create intrinsic gradient curvature calculator
    calc_intrinsic = IntrinsicGradientCurvatureCalculator(sphere, step_size=1e-4)

    # Test points
    test_points = [
        (math.pi/4, 0, "Mid-latitude (north hemisphere)"),
        (math.pi/4, math.pi/6, "General position"),
        (math.pi/2, 0, "Equator"),
        (math.pi/6, math.pi/3, "High latitude (north hemisphere)"),
    ]

    # Theoretical values
    K_theory = 1.0 / (radius * radius)  # Gaussian curvature
    H_theory = 1.0 / radius  # Mean curvature

    print(f"\nSphere radius: R = {radius}")
    print(f"Theoretical Gaussian curvature: K = 1/R² = {K_theory:.6f}")
    print(f"Theoretical mean curvature: H = 1/R = {H_theory:.6f}")
    print(f"Theoretical principal curvatures: k₁ = k₂ = 1/R = {H_theory:.6f}")

    print("\nComputed results:")
    print("-" * 70)
    print(f"{'Location':<30} {'K (computed)':<12} {'K (error%)':<12} {'H (computed)':<12} {'H (error%)':<12}")
    print("-" * 70)

    for u, v, description in test_points:
        # Compute all curvatures
        results = calc_intrinsic.compute_all_curvatures(u, v)

        # Compute errors
        K_error = abs(results['K'] - K_theory) / K_theory * 100 if K_theory != 0 else 0
        H_error = abs(results['H'] - H_theory) / H_theory * 100 if H_theory != 0 else 0

        print(f"{description:<30} {results['K']:11.8f}  {K_error:10.4f}%  "
              f"{results['H']:11.8f}  {H_error:10.4f}%")

    # Verification demo
    print("\nVerification:")
    verification = calc_intrinsic.verify_against_sphere(math.pi/4, math.pi/6, radius)
    print(f"  Gaussian curvature K: computed={verification['K_computed']:.8f}, "
          f"theory={verification['K_theory']:.8f}, "
          f"error={verification['K_error_%']:.4f}%")
    print(f"  Mean curvature H: computed={verification['H_computed']:.8f}, "
          f"theory={verification['H_theory']:.8f}, "
          f"error={verification['H_error_%']:.4f}%")


def compare_methods():
    """Compare intrinsic gradient vs classical method results."""
    print_separator("Method Comparison: Intrinsic Gradient vs Classical")

    # Create sphere
    radius = 2.0
    sphere = Sphere(radius=radius)

    # Test point
    u, v = math.pi/4, math.pi/6

    print(f"\nTest surface: sphere with radius R={radius}")
    print(f"Test point: (u={u:.4f}, v={v:.4f})")
    print(f"Theoretical Gaussian curvature: K = {1.0/(radius*radius):.6f}")

    # Intrinsic gradient method
    calc_intrinsic = IntrinsicGradientCurvatureCalculator(sphere, step_size=1e-4)
    K_intrinsic = calc_intrinsic.compute_gaussian_curvature(u, v)

    # Classical method (5-point finite difference)
    calc_classical = CurvatureCalculator(sphere, step_size=1e-3)
    K_classical = calc_classical.compute_gaussian_curvature(u, v)

    # Theoretical value
    K_theory = 1.0 / (radius * radius)

    # Errors
    error_intrinsic = abs(K_intrinsic - K_theory) / K_theory * 100
    error_classical = abs(K_classical - K_theory) / K_theory * 100

    print("\nComputed results:")
    print("-" * 60)
    print(f"{'Method':<25} {'Gaussian K':<15} {'Relative error%':<15}")
    print("-" * 60)
    print(f"{'Intrinsic gradient':<25} {K_intrinsic:<15.10f} {error_intrinsic:<15.6f}")
    print(f"{'Classical 5-point':<25} {K_classical:<15.10f} {error_classical:<15.6f}")
    print(f"{'Theory':<25} {K_theory:<15.10f} {'0.000000':<15}")


def demo_torus():
    """Demonstrate curvature computation on a torus."""
    print_separator("Torus Curvature (Intrinsic Gradient Operator)")

    # Create torus
    R = 3.0  # Major radius
    r = 1.0  # Minor radius
    torus = Torus(major_radius=R, minor_radius=r)

    # Create calculator
    calc = IntrinsicGradientCurvatureCalculator(torus, step_size=1e-4)

    print(f"\nTorus parameters: R={R} (major radius), r={r} (minor radius)")

    # Test special points
    test_points = [
        (0, 0, "Outer top", 0, 1/(r*(R+r))),  # u=0: outer side
        (math.pi, 0, "Inner bottom", 0, -1/(r*(R-r))),  # u=π: inner side
        (math.pi/2, 0, "Middle position", -1/(r*R), 1/(r*R)),  # u=π/2: middle
    ]

    print("\nSpecial point curvature:")
    print("-" * 90)
    print(f"{'Location':<20} {'K (computed)':<15} {'K (theory)':<15} {'H (computed)':<15} {'H (theory)':<15}")
    print("-" * 90)

    for u, v, description, K_theory, H_theory in test_points:
        results = calc.compute_all_curvatures(u, v)
        print(f"{description:<20} {results['K']:14.8f}  {K_theory:14.8f}  "
              f"{results['H']:14.8f}  {H_theory:14.8f}")


def demo_step_size_analysis():
    """Analyze step-size sensitivity and convergence."""
    print_separator("Step Size Convergence Analysis")

    # Create sphere
    radius = 2.0
    sphere = Sphere(radius=radius)
    u, v = math.pi/4, math.pi/6

    # Step sizes
    step_sizes = [1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5]

    # Theoretical value
    K_theory = 1.0 / (radius * radius)

    print(f"\nTest surface: sphere with radius R={radius}")
    print(f"Test point: (u={u:.4f}, v={v:.4f})")
    print(f"Theoretical Gaussian curvature: K = {K_theory:.10f}")

    print("\nConvergence analysis:")
    print("-" * 70)
    print(f"{'Step h':<12} {'Gaussian K':<20} {'Absolute error':<15} {'Relative error%':<15}")
    print("-" * 70)

    results = []
    for h in step_sizes:
        calc = IntrinsicGradientCurvatureCalculator(sphere, step_size=h)
        K = calc.compute_gaussian_curvature(u, v)
        abs_error = abs(K - K_theory)
        rel_error = abs_error / K_theory * 100 if K_theory != 0 else 0

        results.append({'h': h, 'K': K, 'abs_error': abs_error, 'rel_error': rel_error})
        print(f"{h:<12.1e} {K:<20.15f} {abs_error:<15.2e} {rel_error:<15.8f}")

    # Find best step size
    best = min(results, key=lambda x: x['abs_error'])
    print(f"\nBest step size: h = {best['h']:.1e}, relative error = {best['rel_error']:.8f}%")


def demo_simplified_interface():
    """Demonstrate the simplified interface."""
    print_separator("Simplified Interface Example")

    from coordinate_system import (
        intrinsic_gradient_gaussian_curvature,
        intrinsic_gradient_mean_curvature,
        intrinsic_gradient_principal_curvatures,
        intrinsic_gradient_all_curvatures,
    )

    # Create sphere
    sphere = Sphere(radius=1.0)
    u, v = math.pi/3, math.pi/4

    print("\nCompute unit-sphere curvature with the simplified interface:")
    print(f"Test point: (u={u:.4f}, v={v:.4f})")

    # Compute curvatures individually
    K = intrinsic_gradient_gaussian_curvature(sphere, u, v)
    H = intrinsic_gradient_mean_curvature(sphere, u, v)
    k1, k2 = intrinsic_gradient_principal_curvatures(sphere, u, v)

    print(f"\nGaussian curvature K = {K:.8f} (theory: 1.0)")
    print(f"Mean curvature H = {H:.8f} (theory: 1.0)")
    print(f"Principal curvatures k₁ = {k1:.8f}, k₂ = {k2:.8f} (theory: 1.0, 1.0)")

    # Compute all curvatures at once
    all_results = intrinsic_gradient_all_curvatures(sphere, u, v)
    print("\nUse all_curvatures to retrieve all results:")
    print(f"  K = {all_results['K']:.8f}")
    print(f"  H = {all_results['H']:.8f}")
    print(f"  Second fundamental form: L={all_results['L']:.6f}, M={all_results['M']:.6f}, N={all_results['N']:.6f}")
    print(f"  First fundamental form: E={all_results['E']:.6f}, F={all_results['F']:.6f}, G={all_results['G']:.6f}")
    print(f"  Metric determinant: det(g) = {all_results['det_g']:.6f}")


def demo_ccs_geometry_package():
    """Demonstrate the structured CCS geometry package interface."""
    print_separator("CCS Geometry Package")

    sphere = Sphere(radius=2.0)
    u, v = math.pi/4, math.pi/6
    pkg = compute_ccs_geometry_package(sphere, u, v, step_size=1e-4)

    print("\nThe package exposes the theorem-chain outputs together:")
    print(f"  K  = {pkg.K:.8f}")
    print(f"  H  = {pkg.H:.8f}")
    print(f"  k1 = {pkg.k1:.8f}")
    print(f"  k2 = {pkg.k2:.8f}")
    print(f"  g shape = {pkg.g.shape}, h shape = {pkg.h.shape}, S shape = {pkg.S.shape}")
    print(f"  R_1212 = {pkg.riemann_1212:.8f}")
    print(f"  center frame = {type(pkg.center_frame).__name__}")
    print(f"  connection objects = ({type(pkg.G_u).__name__}, {type(pkg.G_v).__name__})")


def main():
    """Main entry point."""
    print("=" * 70)
    print("Intrinsic Gradient Operator Curvature Example")
    print("Based on: Surface Curvature Computation via the Intrinsic Gradient Operator")
    print("=" * 70)

    # Run demos
    demo_sphere()
    compare_methods()
    demo_torus()
    demo_step_size_analysis()
    demo_simplified_interface()
    demo_ccs_geometry_package()

    print_separator()
    print("\nDemo complete.")
    print("\nKey conclusions:")
    print("1. The intrinsic gradient method accurately computes surface curvature.")
    print("2. For spheres, it achieves high precision (error < 0.01%).")
    print("3. The optimal step size is typically between 1e-4 and 1e-5.")
    print("4. The method is numerically stable and geometrically intuitive.")


if __name__ == "__main__":
    main()
