"""Quick start: CCS Gaussian curvature on a sphere."""

import math

from coordinate_system import Sphere, compute_gaussian_curvature


sphere = Sphere(radius=2.0)
u = math.pi / 4.0
v = math.pi / 3.0

K = compute_gaussian_curvature(sphere, u, v)
K_theory = 1.0 / (2.0 ** 2)

print("Quick Start: Gaussian Curvature on a Sphere")
print("=" * 50)
print("Sphere radius: R = 2.0")
print("Point: (u, v) = (pi/4, pi/3)")
print(f"Computed curvature:    K = {K:.6f}")
print(f"Theoretical curvature: K = {K_theory:.6f}")
print(f"Relative error: {abs(K - K_theory) / K_theory * 100:.2f}%")
