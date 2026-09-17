"""Define a custom surface and compute its CCS curvature package."""

from coordinate_system import Surface, compute_ccs_geometry_package, vec3


class Paraboloid(Surface):
    def position(self, u: float, v: float) -> vec3:
        return vec3(u, v, 0.25 * (u * u + v * v))

    def derivs(self, u: float, v: float):
        return (
            vec3(1.0, 0.0, 0.5 * u),
            vec3(0.0, 1.0, 0.5 * v),
            vec3(0.0, 0.0, 0.5),
            vec3(0.0, 0.0, 0.0),
            vec3(0.0, 0.0, 0.5),
        )


pkg = compute_ccs_geometry_package(Paraboloid(), 0.4, -0.2, step_size=1e-4)
print("K:", pkg.K)
print("H:", pkg.H)
print("g:\n", pkg.g)
print("h:\n", pkg.h)
