"""Inspect the CCS theorem-chain geometry package."""

import math

from coordinate_system import Sphere, compute_ccs_geometry_package


sphere = Sphere(radius=2.0)
pkg = compute_ccs_geometry_package(sphere, math.pi / 4.0, math.pi / 3.0, step_size=1e-4)

print("Center frame:", pkg.center_frame)
print("K:", pkg.K)
print("H:", pkg.H)
print("k1, k2:", pkg.k1, pkg.k2)
print("g:\n", pkg.g)
print("h:\n", pkg.h)
print("S:\n", pkg.S)
print("R_1212:", pkg.riemann_1212)
