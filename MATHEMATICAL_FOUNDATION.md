# Mathematical Foundation of Coordinate Systems

**A Comprehensive Guide to 3D Coordinate System Mathematics**

Author: PanGuoJun (romeosoft)
Version: 1.2.0
License: MIT

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Vectors](#2-vectors)
3. [Coordinate Systems](#3-coordinate-systems)
4. [Vector Operations](#4-vector-operations)
5. [Quaternions](#5-quaternions)
6. [Coordinate Transformations](#6-coordinate-transformations)
7. [Euler Angles](#7-euler-angles)
8. [Interpolation Methods](#8-interpolation-methods)
9. [Advanced Topics](#9-advanced-topics)
10. [Practical Applications](#10-practical-applications)

---

## 1. Introduction

### 1.1 What is a Coordinate System?

A **coordinate system** is a mathematical framework that allows us to describe positions and orientations in space. In 3D graphics, physics simulation, and robotics, we use coordinate systems to:

- Represent object positions and rotations
- Transform between different reference frames
- Perform hierarchical transformations (parent-child relationships)
- Implement camera systems and view transformations

### 1.2 The Fundamental Coordinate System Equation

The core concept of this library is expressed by the fundamental equation:

```
P_world = P_local * CoordinateFrame
```

Where:
- `P_local`: A point in local space
- `CoordinateFrame`: A 3D coordinate system (position + rotation + scale)
- `P_world`: The resulting point in world space

This simple equation enables all coordinate transformations in 3D space.

---

## 2. Vectors

### 2.1 Vector Definition

A 3D vector represents a direction and magnitude in space:

```
v = (x, y, z)
```

**Properties:**
- Direction: normalized vector (length = 1)
- Magnitude: length of the vector
- Components: x, y, z values

### 2.2 Vector Operations

#### Addition and Subtraction

```python
v1 + v2 = (x1 + x2, y1 + y2, z1 + z2)
v1 - v2 = (x1 - x2, y1 - y2, z1 - z2)
```

**Geometric meaning:** Vector addition follows the parallelogram rule.

#### Scalar Multiplication

```python
k * v = (k*x, k*y, k*z)
```

**Geometric meaning:** Scales the vector's magnitude without changing its direction.

#### Dot Product (Inner Product)

```python
v1 · v2 = x1*x2 + y1*y2 + z1*z2 = |v1| * |v2| * cos(θ)
```

**Properties:**
- Returns a scalar value
- Measures "similarity" of two vectors
- Used to calculate angles between vectors
- If dot product = 0, vectors are perpendicular

**Applications:**
- Angle calculation: `θ = arccos(v1·v2 / (|v1|*|v2|))`
- Projection: `proj = (v1·v2 / |v2|²) * v2`
- Determining if vectors point in same/opposite directions

#### Cross Product (Vector Product)

```python
v1 × v2 = (y1*z2 - z1*y2, z1*x2 - x1*z2, x1*y2 - y1*x2)
```

**Properties:**
- Returns a vector perpendicular to both input vectors
- Magnitude: `|v1 × v2| = |v1| * |v2| * sin(θ)`
- Direction: follows right-hand rule (or left-hand rule depending on coordinate system)
- Anti-commutative: `v1 × v2 = -(v2 × v1)`

**Applications:**
- Finding surface normals
- Calculating torque in physics
- Determining polygon orientation

### 2.3 Vector Length and Normalization

#### Length (Magnitude)

```python
|v| = √(x² + y² + z²)
```

#### Normalization

Convert a vector to unit length (magnitude = 1):

```python
v_normalized = v / |v| = (x/|v|, y/|v|, z/|v|)
```

**Important:** Normalized vectors represent pure direction, removing magnitude information.

---

## 3. Coordinate Systems

### 3.1 Coordinate System Components

A complete 3D coordinate system consists of:

1. **Origin (o)**: Position vector defining the coordinate system's center
2. **Basis vectors (ux, uy, uz)**: Three orthonormal vectors defining the axes
3. **Scale (s)**: Scale factors for each axis

Mathematical representation:

```
CoordSystem = {o, ux, uy, uz, s}
```

Where:
- `ux · uy = 0` (perpendicular)
- `uy · uz = 0` (perpendicular)
- `uz · ux = 0` (perpendicular)
- `|ux| = |uy| = |uz| = 1` (unit length)

### 3.2 World vs Local Coordinates

#### World Coordinates
The global reference frame, typically defined as:
```
World = {(0,0,0), (1,0,0), (0,1,0), (0,0,1), (1,1,1)}
```

#### Local Coordinates
Any coordinate system relative to world or another coordinate system.

### 3.3 Coordinate System Construction

#### From Position Only
```python
frame = coord3(x, y, z)  # Origin at (x,y,z), axes aligned with world
```

#### From Position and Rotation (Quaternion)
```python
frame = coord3(position, quaternion, scale)
```

#### From Three Basis Vectors
```python
frame = coord3.from_axes(ux, uy, uz)
```

---

## 4. Vector Operations in Detail

### 4.1 Vector Projection

Project vector `a` onto vector `b`:

```python
proj_b(a) = (a · b / |b|²) * b
```

**Geometric meaning:** The component of `a` that points in the direction of `b`.

### 4.2 Vector Reflection

Reflect vector `v` across a surface with normal `n`:

```python
v_reflected = v - 2 * (v · n) * n
```

**Application:** Ray tracing, physics simulations (bouncing objects).

### 4.3 Vector Distance

Distance between two points:

```python
distance = |v2 - v1| = √((x2-x1)² + (y2-y1)² + (z2-z1)²)
```

---

## 5. Quaternions

### 5.1 Quaternion Definition

A quaternion is a 4-component number used to represent 3D rotations:

```
q = w + xi + yj + zk
```

In component form: `q = (w, x, y, z)`

Where:
- `w`: scalar (real) part
- `(x, y, z)`: vector (imaginary) part
- `i² = j² = k² = ijk = -1`

### 5.2 Why Quaternions?

**Advantages over Euler angles:**
1. No gimbal lock
2. Smooth interpolation (slerp)
3. More efficient rotation composition
4. Numerically stable

**Advantages over rotation matrices:**
1. Compact (4 values vs 9 values)
2. Easy normalization
3. Natural interpolation

### 5.3 Quaternion Operations

#### Quaternion Multiplication

```python
q1 * q2 = (
    w1*w2 - x1*x2 - y1*y2 - z1*z2,
    w1*x2 + x1*w2 + y1*z2 - z1*y2,
    w1*y2 - x1*z2 + y1*w2 + z1*x2,
    w1*z2 + x1*y2 - y1*x2 + z1*w2
)
```

**Note:** Quaternion multiplication is NOT commutative: `q1 * q2 ≠ q2 * q1`

#### Quaternion Conjugate

```python
q* = (w, -x, -y, -z)
```

#### Quaternion Inverse

For unit quaternions: `q⁻¹ = q*`

For general quaternions: `q⁻¹ = q* / |q|²`

#### Quaternion Normalization

```python
q_normalized = q / |q|
where |q| = √(w² + x² + y² + z²)
```

### 5.4 Rotation Representation

#### From Angle-Axis to Quaternion

Given rotation angle `θ` around axis `n = (nx, ny, nz)`:

```python
q = (cos(θ/2), sin(θ/2)*nx, sin(θ/2)*ny, sin(θ/2)*nz)
```

**Example:** 90° rotation around Z-axis:
```python
θ = π/2, n = (0, 0, 1)
q = (cos(π/4), 0, 0, sin(π/4)) = (0.707, 0, 0, 0.707)
```

#### Rotating a Vector

To rotate vector `v` by quaternion `q`:

```python
v_rotated = q * v * q⁻¹
```

Where `v` is treated as quaternion `(0, vx, vy, vz)`.

### 5.5 Quaternion from Two Vectors

Find rotation that transforms `v1` to `v2`:

```python
axis = v1 × v2
angle = arccos(v1 · v2 / (|v1| * |v2|))
q = quaternion(angle, axis)
```

---

## 6. Coordinate Transformations

### 6.1 The Fundamental Transformation Equation

Transform a point from local space to world space:

```
P_world = P_local * Coord
```

Expanded form:

```
P_world = Coord.o + P_local.x * Coord.ux + P_local.y * Coord.uy + P_local.z * Coord.uz
```

### 6.2 Inverse Transformation

Transform from world space to local space:

```
P_local = P_world / Coord
```

Mathematically:

```
P_local.x = (P_world - Coord.o) · Coord.ux
P_local.y = (P_world - Coord.o) · Coord.uy
P_local.z = (P_world - Coord.o) · Coord.uz
```

### 6.3 Coordinate System Composition

Combine two coordinate transformations:

```
CoordResult = Coord1 * Coord2
```

**Geometric meaning:** Apply `Coord2` transformation, then apply `Coord1` transformation.

**Application:** Parent-child hierarchies in scene graphs.

### 6.4 Hierarchical Transformations

Example: Robot arm with multiple joints

```python
# Shoulder joint
shoulder = coord3(0, 0, 0, quat(angle1, vec3(0, 1, 0)))

# Elbow joint (relative to shoulder)
elbow_local = coord3(5, 0, 0, quat(angle2, vec3(0, 1, 0)))

# Elbow in world space
elbow_world = elbow_local * shoulder

# Hand (relative to elbow)
hand_local = coord3(5, 0, 0, quat(angle3, vec3(0, 1, 0)))

# Hand in world space
hand_world = hand_local * elbow_world
```

---

## 7. Euler Angles

### 7.1 Euler Angle Definition

Three sequential rotations around coordinate axes:

- **Pitch** (X-axis): Up/down rotation
- **Yaw** (Y-axis): Left/right rotation
- **Roll** (Z-axis): Tilt rotation

Convention used (intrinsic rotations):
```
Rotation = Yaw * Pitch * Roll
```

### 7.2 Euler Angles to Quaternion

```python
cy = cos(yaw * 0.5)
sy = sin(yaw * 0.5)
cp = cos(pitch * 0.5)
sp = sin(pitch * 0.5)
cr = cos(roll * 0.5)
sr = sin(roll * 0.5)

q.w = cr * cp * cy + sr * sp * sy
q.x = sr * cp * cy - cr * sp * sy
q.y = cr * sp * cy + sr * cp * sy
q.z = cr * cp * sy - sr * sp * cy
```

### 7.3 Gimbal Lock Problem

**Definition:** Loss of one degree of freedom when two rotation axes align.

**Example:** When pitch = ±90°, yaw and roll rotations become equivalent.

**Solution:** Use quaternions instead of Euler angles for intermediate calculations.

---

## 8. Interpolation Methods

### 8.1 Linear Interpolation (LERP)

Interpolate between two values linearly:

```python
result = (1 - t) * a + t * b
```

Where `t ∈ [0, 1]`:
- `t = 0` → result = a
- `t = 0.5` → result = midpoint
- `t = 1` → result = b

**Vector LERP:**
```python
v_lerp = (1-t) * v1 + t * v2
```

**Caution:** LERP of quaternions requires normalization and doesn't preserve constant angular velocity.

### 8.2 Spherical Linear Interpolation (SLERP)

For quaternions, SLERP provides smooth rotation interpolation:

```python
q_slerp = sin((1-t)*θ)/sin(θ) * q1 + sin(t*θ)/sin(θ) * q2
```

Where `θ = arccos(q1 · q2)`

**Properties:**
- Constant angular velocity
- Shortest path on quaternion sphere
- Numerically stable

**Implementation:**

```python
def slerp(q1, q2, t):
    dot = q1.dot(q2)

    # Ensure shortest path
    if dot < 0:
        q2 = -q2
        dot = -dot

    # If quaternions are very close, use LERP
    if dot > 0.9995:
        return normalize(lerp(q1, q2, t))

    # Calculate interpolation
    theta = arccos(dot)
    return (sin((1-t)*theta)/sin(theta))*q1 + (sin(t*theta)/sin(theta))*q2
```

### 8.3 Normalized Linear Interpolation (NLERP)

Fast approximation of SLERP:

```python
q_nlerp = normalize((1-t) * q1 + t * q2)
```

**Advantages:**
- Much faster than SLERP
- Good approximation for small angles

**Disadvantages:**
- Non-constant angular velocity
- Less accurate for large rotations

### 8.4 Coordinate System Interpolation

Interpolate between two coordinate systems:

```python
def coord_lerp(c1, c2, t):
    pos = lerp(c1.o, c2.o, t)           # Position: linear
    rot = slerp(c1.Q(), c2.Q(), t)      # Rotation: spherical
    scale = lerp(c1.s, c2.s, t)         # Scale: linear
    return coord3(pos, rot, scale)
```

---

## 9. Advanced Topics

### 9.1 Rotation Matrices

Convert quaternion to 3×3 rotation matrix:

```python
R = [
    [1-2(y²+z²),  2(xy-wz),    2(xz+wy)   ],
    [2(xy+wz),    1-2(x²+z²),  2(yz-wx)   ],
    [2(xz-wy),    2(yz+wx),    1-2(x²+y²) ]
]
```

### 9.2 Transformation Matrix (4×4)

Homogeneous transformation matrix:

```python
T = [
    [ux.x  uy.x  uz.x  o.x]
    [ux.y  uy.y  uz.y  o.y]
    [ux.z  uy.z  uz.z  o.z]
    [0     0     0     1  ]
]
```

This matrix combines rotation and translation in a single operation.

### 9.3 Inverse Transformation

For orthonormal coordinate system:

```python
T⁻¹ = [
    [ux.x  ux.y  ux.z  -o·ux]
    [uy.x  uy.y  uy.z  -o·uy]
    [uz.x  uz.y  uz.z  -o·uz]
    [0     0     0     1    ]
]
```

**Note:** Transpose of rotation part + transformed negative translation.

### 9.4 Coordinate System Derivatives

Rate of change of coordinate system (velocity):

```python
velocity = d(coord)/dt
angular_velocity = 2 * dq/dt * q⁻¹
```

Useful for physics simulation and motion control.

---

## 10. Practical Applications

### 10.1 Camera System

**Look-At Matrix:**

Create a camera coordinate system looking at a target:

```python
def create_look_at(eye, target, up=vec3(0, 1, 0)):
    # Forward direction
    forward = (target - eye).normcopy()

    # Right direction
    right = up.cross(forward).normcopy()

    # Corrected up direction
    up_corrected = forward.cross(right)

    # Create coordinate system
    return coord3.from_axes(right, up_corrected, forward)
```

**Camera transformations:**
```python
# View matrix (world to camera)
view_matrix = inverse(camera_coord)

# Point in camera space
point_camera = world_point / camera_coord
```

### 10.2 Physics Simulation

**Rigid Body Transformation:**

```python
class RigidBody:
    def __init__(self):
        self.coord = coord3()           # Position and orientation
        self.velocity = vec3(0, 0, 0)   # Linear velocity
        self.angular_vel = vec3(0, 0, 0) # Angular velocity (axis-angle)

    def update(self, dt):
        # Update position
        self.coord.o += self.velocity * dt

        # Update rotation
        if self.angular_vel.length() > 0:
            angle = self.angular_vel.length() * dt
            axis = self.angular_vel.normcopy()
            rotation = quat(angle, axis)
            self.coord.rot(rotation)
```

### 10.3 Skeletal Animation

**Bone Hierarchy:**

```python
class Bone:
    def __init__(self, name, parent=None):
        self.name = name
        self.parent = parent
        self.local_transform = coord3()  # Relative to parent
        self.world_transform = coord3()  # Absolute position

    def update_world_transform(self):
        if self.parent:
            self.world_transform = self.local_transform * self.parent.world_transform
        else:
            self.world_transform = self.local_transform
```

**Animation Blending:**

```python
def blend_animations(anim1, anim2, weight):
    """Blend two animation poses"""
    result = {}
    for bone_name in anim1.keys():
        pose1 = anim1[bone_name]
        pose2 = anim2[bone_name]

        # Interpolate position, rotation, scale
        result[bone_name] = coord_lerp(pose1, pose2, weight)

    return result
```

### 10.4 Scene Graph

**Hierarchical Scene Structure:**

```python
class SceneNode:
    def __init__(self, name):
        self.name = name
        self.local_coord = coord3()      # Local transformation
        self.world_coord = coord3()      # World transformation
        self.parent = None
        self.children = []

    def add_child(self, child):
        child.parent = self
        self.children.append(child)

    def update_transforms(self):
        """Update world transforms recursively"""
        if self.parent:
            self.world_coord = self.local_coord * self.parent.world_coord
        else:
            self.world_coord = self.local_coord

        for child in self.children:
            child.update_transforms()
```

---

## Mathematical Notation Summary

| Symbol | Meaning |
|--------|---------|
| `v, u` | Vectors |
| `q` | Quaternion |
| `c, C` | Coordinate system |
| `·` | Dot product |
| `×` | Cross product |
| `*` | Multiplication (scalar or composition) |
| `\|v\|` | Vector magnitude |
| `θ` | Angle |
| `t` | Interpolation parameter [0, 1] |
| `ux, uy, uz` | Coordinate system basis vectors |
| `o` | Origin/position |
| `s` | Scale vector |

---

## References and Further Reading

1. **Quaternions:** Ken Shoemake, "Animating rotation with quaternion curves", ACM SIGGRAPH 1985
2. **Coordinate Systems:** John Vince, "Mathematics for Computer Graphics"
3. **3D Transformations:** F. S. Hill Jr., "Computer Graphics Using OpenGL"
4. **Interpolation:** Erik B. Dam et al., "Quaternions, Interpolation and Animation"

---

## Key Takeaways

1. **Vectors** are the foundation - master dot product, cross product, and normalization
2. **Quaternions** avoid gimbal lock and provide smooth interpolation
3. **Coordinate systems** enable hierarchical transformations through multiplication
4. **The fundamental equation** `P_world = P_local * Coord` is all you need for transformations
5. **SLERP** provides the smoothest rotation interpolation
6. **Hierarchical transformations** enable complex scenes with parent-child relationships

---

**This mathematical foundation is implemented efficiently in the `coordinate_system` library, providing high-performance tools for 3D graphics, physics, and robotics applications.**

For code examples and API reference, see [README.md](README.md).

---

*Written by PanGuoJun (romeosoft) - 2025*
*License: MIT - Free to use and modify*
