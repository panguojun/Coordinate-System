# coordinate_system API Reference

Most public APIs are exported from the package root:

```python
import coordinate_system as cs
from coordinate_system import vec3, quat, coord3
```

Interactive discovery:

```python
print(cs.__version__)
print(cs.__all__)
help(cs.coord3)
help(cs.compute_ccs_geometry_package)
```

## Core Objects

### `vec3`

3D vector.

Constructors:

```python
v = cs.vec3()
v = cs.vec3(x, y, z)
```

Fields: `x`, `y`, `z`.

Common methods: `dot`, `cross`, `cross_right`, `len`, `length`, `sqrlen`,
`len_squared`, `normalized`, `normalize`, `normcopy`, `project`, `reflect`,
`distance`, `mean`, `volum`, `hash`, `isINF`, `flipX`, `flipY`, `flipZ`.

Static helpers/constants: `vec3.min3`, `vec3.max3`, `vec3.rnd`, `vec3.lerp`,
`vec3.angle`, `vec3.ZERO`, `vec3.UX`, `vec3.UY`, `vec3.UZ`.

```python
a = cs.vec3(1, 0, 0)
b = cs.vec3(0, 1, 0)
print(a.cross(b))
```

### `vec2`

2D vector.

Constructors: `vec2()`, `vec2(x, y)`, `vec2(value)`.

Fields: `x`, `y`.

Common methods: `dot`, `cross`, `len`, `length`, `sqrlen`, `normalized`,
`normalize`, `angle`, `rot`, `rotcopy`, `roted`, `distance`, `lerp`, `isINF`,
`vec2.ang_len`.

### `quat`

Quaternion rotation.

Constructors:

```python
q = cs.quat()
q = cs.quat(w, x, y, z)
q = cs.quat(pitch, yaw, roll)
q = cs.quat(angle, axis)
q = cs.quat(from_vec, to_vec)
```

Common methods: `normalize`, `normalized`, `angle`, `axis`, `conj`,
`conjcopy`, `inverse`, `length`, `dot`, `angle_to`, `rotate`, `to_eulers`,
`to_angle_axis`, `xyz`, `set_angle`, `rotate_x`, `rotate_y`, `rotate_z`,
`is_finite`, `from_vectors`, `from_eulers`, `log`.

Static helpers: `quat.slerp`, `quat.nlerp`, `quat.from_euler`,
`quat.from_axis_angle`.

```python
import math
q = cs.quat.from_axis_angle(cs.vec3.UZ, math.pi / 2)
print(q.rotate(cs.vec3.UX))
```

### `coord3`

3D coordinate-system object. This is the software-level CCS object.

Fields: `o`/`p`, `x`, `y`, `z`, `ux`, `uy`, `uz`, `s`.

Constructors:

```python
C = cs.coord3()
C = cs.coord3(x, y, z)
C = cs.coord3(position)
C = cs.coord3(position, rotation)
C = cs.coord3(position, rotation, scale)
C = cs.coord3(origin, ux, uy, uz)
C = cs.coord3(origin, ux, uy, uz, scale)
```

Common methods: `Q`, `R`, `P`, `V`, `normalize`, `to_eulers`, `rot`,
`equal_dirs`, `dump`, `lie_cross`, `grad`, `inverse`, `inversed`, `reverse`,
`reversed`, `distance_to`, `rotation_distance_to`, `pose`, `pos`, `to_world`,
`to_local`, `VX`, `VY`, `VZ`, `X`, `Y`, `Z`, `ucoord`, `UC`, `VC`,
`compute_metric_det`, `metric`, `metric_det`.

Static helpers: `coord3.from_axes`, `coord3.from_angle`, `coord3.look_at`,
`coord3.from_forward`, `coord3.from_eulers`, `coord3.lerp`, `coord3.slerp`,
`coord3.identity`, `coord3.zero`, `coord3.from_position`,
`coord3.from_rotation`, `coord3.lie_bracket`.

```python
C = cs.coord3(cs.vec3(10, 0, 0), cs.quat(), cs.vec3(2, 2, 2))
world = C.to_world(cs.vec3(1, 0, 0))
print(world)
print(C.to_local(world))
```

## Constants and Utilities

Constants: `ZERO3`, `UNITX`, `UNITY`, `UNITZ`, `ONE3`, `ONE4`, `ONEC`.

Functions: `lerp`, `cross`, `cross_right`, `set_handedness`,
`get_handedness`, `is_left_handed`, `is_right_handed`.

## CCS Surface Geometry

### Surfaces

Classes: `Surface`, `Sphere`, `Torus`.

Custom surface:

```python
class MySurface(cs.Surface):
    def position(self, u, v):
        return cs.vec3(u, v, u * v)
```

Optional high-accuracy hook:

```python
def derivs(self, u, v):
    return r_u, r_v, r_uu, r_uv, r_vv
```

### Data Classes

`MetricTensor`: `from_surface`, `determinant`, `inverse`, `as_matrix`.

`GradientResult`: fields `dn`, `direction`.

`CCSGeometryPackage`: fields `center_frame`, `G_u`, `G_v`, `g`, `h`, `S`,
`K`, `H`, `k1`, `k2`, `normal`, `riemann_1212`; method `as_dict()`.

### Curvature Functions

- `compute_gaussian_curvature(surface, u, v, step_size=1e-3)`
- `compute_mean_curvature(surface, u, v, step_size=1e-3)`
- `compute_riemann_curvature(surface, u, v, step_size=1e-3)`
- `compute_curvature_tensor(surface, u, v, step_size=1e-3)`
- `compute_all_curvatures(surface, u, v, step_size=1e-3)`
- `compute_ccs_geometry_package(surface, u, v, step_size=1e-3)`

```python
import math
sphere = cs.Sphere(2.0)
pkg = cs.compute_ccs_geometry_package(sphere, math.pi / 4, math.pi / 3)
print(pkg.K, pkg.H)
print(pkg.g, pkg.h, pkg.S)
```

### Intrinsic Gradient

Classes: `IntrinsicGradientOperator`, `IntrinsicGradientCurvatureCalculator`,
`LieGroupCurvatureCalculator`, `CurvatureCalculator`.

`IntrinsicGradientOperator` methods: `calc_intrinsic_frame`, `compute_both`,
`compute_u`, `compute_v`, `compute_connection_matrices`.

Functions: `compute_intrinsic_gradient`, `compute_connection_matrices`.

Calculator methods: `compute_gaussian_curvature`, `compute_mean_curvature`,
`compute_riemann_curvature`, `compute_principal_curvatures`,
`compute_all_curvatures`.

### Classical and Compatibility Curvature

Classical functions: `gaussian_curvature_classical`,
`mean_curvature_classical`, `principal_curvatures_classical`,
`all_curvatures_classical`.

Default aliases: `gaussian_curvature`, `mean_curvature`,
`principal_curvatures`, `all_curvatures`.

Legacy aliases: `gaussian_curvature_lie`,
`intrinsic_gradient_gaussian_curvature`, `intrinsic_gradient_mean_curvature`,
`intrinsic_gradient_principal_curvatures`, `intrinsic_gradient_all_curvatures`.

Comparison/helper functions: `compare_methods`, `derivative_5pt`,
`derivative_2nd_5pt`, `richardson_extrapolation`.

### `CCS` Wrapper

```python
ccs = cs.CCS(step_size=1e-4)
ccs.geometry_package(surface, u, v)
ccs.gaussian(surface, u, v)
ccs.mean(surface, u, v)
ccs.riemann(surface, u, v)
ccs.connection_matrices(surface, u, v)
ccs.curvature_tensor(surface, u, v)
```

## Analytic Surface Constraints

Classes: `AnalyticParametricSurface`, `AnalyticSphere`, `AnalyticEllipsoid`,
`AnalyticTorus`, `AnalyticCylinder`, `SurfaceCoord`,
`IntersectionTraceResult`.

Functions: `tangent_basis`, `coord_from_axes`, `solve_feature_coord`,
`interpolate_surface_coords`, `trace_surface_intersection`,
`build_surface_from_spec`.

## Curves

Classes: `InterpolatedCurve`.

Functions: `frame_from_pose`, `polyline_arc_lengths`,
`resample_curve_equal_arc_length`, `frame_pair_hermite_curve`,
`estimate_arc_tangent_scale`, `frame_pair_arc_curve`, `mirror_curve`,
`generate_frenet_frames`, `frame_field_spline`, `frame_field_spline_c2`,
`reconstruct_curve_from_polygon`, `compute_curvature_profile`, `catmull_rom`,
`squad_interp`.

## Curve Intersections

Classes: `CurveIntersectionPoint`, `CurveIntersectionResult`.

Functions: `closest_points_on_segments`, `frame_at_curve_hit`,
`intersect_polyline_curves`, `sample_parametric_curve`,
`intersect_parametric_curves`, `intersect_curve_with_implicit_surface`.

## NURBS

Classes: `NURBSCurve`, `NURBSCurve2D`, `NURBSFrameCurve`, `NURBSSurface`,
`NURBSIntersectionFrame`.

Functions: `open_uniform_knot_vector`, `make_nurbs_curve`,
`make_nurbs_surface`, `sample_nurbs_frame_curve`, `nurbs_arc_length`,
`nurbs_curvature_profile`, `intersect_nurbs_curves`,
`intersect_nurbs_with_implicit_surface`, `solve_curve_intersections`,
`solve_curve_surface_intersections`, `solve_surface_intersection_curve`,
`correct_nurbs_surface_intersection`.

## Helicity Phase

Constants/classes: `TAU`, `HelicityPhaseTrajectory`,
`CompactBranchTrajectory`.

Functions: `wrap_phase`, `branch_index`, `branch_residual`,
`branch_transition_indices`, `cumulative_pose_exposure`,
`integrate_helicity_phase`, `integrate_compact_branch`.

## Spectral Geometry

Classes: `FourierFrame`, `FourierFrameSpectrum`, `IntrinsicGradient`,
`CurvatureFromFrame`, `BerryPhase`, `ChernNumber`, `SpectralDecomposition`,
`HeatKernel`, `FrequencyProjection`, `FrequencyBandState`.

Functions/constants: `spectral_transform`, `inverse_spectral_transform`,
`HBAR`, `GPU_AVAILABLE`.

## Complex Frames

Classes: `ComplexFrame`, `SU3Component`, `GaugeConnection`, `FieldStrength`,
`SymmetryBreakingPotential`.

## CoordScript

Exceptions: `CoordScriptError`, `CoordScriptNameError`,
`CoordScriptRequirementError`.

Result classes: `RequirementResult`, `AuditResult`, `CoordScriptResult`.

Functions: `parse_coordscript`, `eval_coordscript`, `eval_coordscript_file`.

```python
result = cs.eval_coordscript("""
p = vec3(1, 2, 3)
C = coord(p)
require_close(norm(origin(C)), norm(p), 1e-9)
""")
print(result.audit.passed)
```

## Recommended Learning Path

1. `examples/quickstart.py`
2. `examples/ccs_geometry_package.py`
3. `examples/custom_surface.py`
4. This API reference for the public API index
