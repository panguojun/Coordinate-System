![image](https://github.com/user-attachments/assets/066f8e50-1146-4ebb-b319-cabc8b3f7add)

# The Coordinate System (Coord) Framework

## History of Coordinate Systems
The concept of coordinate systems traces back to René Descartes, who sought to use geometry to describe celestial motion. However, his methods lacked the precision required for exact calculations. Long before Descartes, early civilizations already had notions of coordinate-like references—particularly the idea of a "world center."

During the Hellenistic period, Ptolemaic cosmology placed Earth at the center of the universe, while Copernicus later shifted this central reference to the Sun. The key difference between these models was not just the choice of origin but the mathematical framework they enabled. By repositioning the center at the Sun, scientists recognized the need for a dynamic, motion-based mathematical-physical system. This realization paved the way for calculus, equations of motion, and Newton's laws of inertia—cornerstones of modern science.

Thus, the choice of coordinate system profoundly influences both worldview and computational paradigms. Einstein’s relativity theory, for instance, can be viewed as a consequence of extending coordinates from flat Euclidean space to curved manifolds. Moreover, it appears that all precisely calculable problems ultimately reduce to coordinate transformations. Approximate methods—such as statistical approaches in quantum mechanics and thermodynamics—remain necessary where exact solutions are intractable (though modern techniques like Density Functional Theory (DFT) have achieved notable, if still imperfect, success).

## Mathematical Foundation
If we study differential geometry, at the level of differential elements, the coordinate system can be linearized. In this way, a concept of a certain universal coordinate system object can be formed, serving as a ruler or, alternatively, understood as the dimension of space.

The Coordinate System (or Frame), referred to here as the Coord object, is a mathematical construct rooted in group theory that defines a coordinate system in three-dimensional space. In physics, such a structure is commonly known as a reference frame, while in differential geometry, it is often called a frame field or moving frame, borrowing terminology from classical mechanics.

From a group-theoretic perspective, a coordinate system (or its simplified form, a coordinate) can be treated as an algebraic object capable of participating in group operations. The Coord object unifies these concepts, allowing both coordinate systems and individual coordinates to serve as elements in algebraic operations, such as multiplication and division.

By extending coordinate systems with arithmetic operations (addition, subtraction, multiplication, and division), the Coord object enables direct differential calculus, eliminating the need for cumbersome exterior calculus formulations. This approach provides an intuitive geometric interpretation of operations like vector division, which traditionally require complex tensor algebra. Moreover, it simplifies advanced differential geometry concepts such as connections (affine or Levi-Civita) and curvature tensors, offering a unified and geometrically intuitive framework for computations involving coordinate systems.

## Applications
World ↔ Local Coordinate Transformations: Seamlessly convert between global and local reference frames.

Multi-Node Hierarchies: Efficiently manage transformations in complex systems (e.g., robotics, computer graphics).

Differential Geometry & Physics: Streamline computations involving curvature, parallel transport, and dynamic reference frames.

The Coord object thus serves as a powerful abstraction, bridging algebraic operations, differential calculus, and geometric intuition in a computationally elegant manner.

## Structure of the Coordinate System

In C++, a coordinate system in three-dimensional space is defined by an origin, three directional axes, and three scaling components as follows:

```
struct coord {
    vec3 ux, uy, uz;   // Three basis vectors
    vec3 s;            // Scaling
    vec3 o;            // Origin
};
```

## Constructing a Coordinate System

A coordinate system can be constructed using three axes or Euler angles as follows:

```
coord C1(vec3 o);
coord C2(vec3 ux, vec3 uy, vec3 uz);
coord C3(vec3 o, vec3 s, quat q); 
```

## Multiplication and Division Operations

Multiplication and division operations are provided to transform vectors from one coordinate system to another and to project vectors from a parent coordinate system to a local one. For example, to transform a vector V1 from a local coordinate system C1 to a parent coordinate system C0, we can use the following operation:

```
V0 = V1 * C1
```

To project a vector V0 from a parent coordinate system C0 to a local coordinate system C1, we can use the following operation:

```
V1 = V0 / C1
```

## Common Scenarios

Coord can be applied in various scenarios, such as converting between world and local coordinate systems and using it in multi-node hierarchies. Here are some examples:

1. Convert a vector Vw from a world coordinate system to a local coordinate system C:

```
VL = Vw / C   
Vw = VL * C 
```

2. Convert between world and local coordinate systems:

```
C = C3 * C2 * C1
Vw = VL * C
VL = Vw / C
```

3. Use in multi-node hierarchies:

```
V1 = V4 * C4 * C3 * C2 
V4 = V1 / C2 / C3 / C4
```

4. Convert between parallel coordinate systems:

```
C0 { C1, C2 }
V2 = V1 * C1 / C2
```

5. More operations:

Scalar multiplication:

```
C * k = {C.o, C.s * k, C.u}
where: C.u = {C.ux, C.uy, C.uz}
```

Quaternion multiplication:

```
C0 = C1 * q1 
C1 = C0 / q1
q0 = q1 * C1
q1 = q0 / C1
```

Vector addition:

```
C2 = C1 + o
Where C2 = {C1.o + o, C1.v}, C1.v = {C1.ux*C1.s.x, C1.uy*C1.s.y, C1.uz*C1.s.z}
```

Coordinate Gradient:
```
G = C1 / C2 - I
Where C1 and C2 are coordinate systems on two points on a unit length distance.
```
## Coordinate System Differentiation

Coord can be used to differentiate coordinate systems in space in three ways:

Gradient: 

```
▽f = (u * df * Cuv) / Dxyz
Where:
    Cuv  = {u, v, 0}
    Dxyz = {ux * dx, uy * dy, uz * dz}
```

Divergence:

```
▽ ∙ F = dF / Dxyz ∙ Ic
Where: Ic = {ux, uy, uz}
```

Curl:

```
▽ x F = dF / Dxyz x Ic
```

## Differential Geometry Framework

### Connection Calculus
```cpp
// Finite connection between frames
G = C2 / C1 - I;

// Intrinsic connection (embedded surfaces)
G_intrinsic = C2 / C1 / c2 - I / c1;
```
Where:
- `C1,C2` are 3D coordinate frames along the surface
- `c1,c2` are mappings from intrinsic to global coordinates (e.g., cone development coordinates)

### Calculate the Space Curvature
Coord can transport vectors from a natural coordinate system to a curved coordinate system in a curved space. The curvature can be determined by comparing two paths projected onto the u and v curves, which is done using Gu and Gv. Gu and Gv represent the gradients of rotational changes of vectors along the u and v. By using a coordinate system, the spatial curvature can be calculated, and the curvature tensor in the u,v coordinate system is given by:

```
Ruv = Gu*Gv - Gv*Gu - G[u,v]

where:  Gu = C2 / C1 - I
        Connection vector: [u, v] (Lie bracket operation)
        W = Wu + Wv = [u, v]
        G[u,v] = Gu*Wu + Gv*Wv
```

## Lie-Theoretic Interpretation

The Coord framework naturally embeds Lie theory:
- **Group Multiplication**: Represents SE(3) action
- **Lie Algebra**: The tangent space at identity
- **Bracket Operation**: 
  ```cpp
  [C1, C2] = C1*C2 - C2*C1
  ```
- **Exponential Map**: From algebra to group

This provides a unified representation for:
- Rigid transformations (SE(3))
- Conformal transformations
- Gauge transformations

## Implementation and Usage

## Python Installation

To use the coordinate_system in Python(3.13), you can easily install it via pip:

```bash
pip install coordinate_system
```

```python
from coordinate_system import vec3,quat,coord3
a = coord3(0,0,1,0,45,0);
b = coord3(1,0,0,45,0,0);
a*=b;
print(a);
```
This will allow you to leverage the powerful features of the coord3 class in Python for your mathematical and computational needs.

## Conclusion

The Coord framework provides a unified language for:
- Geometric transformations
- Differential geometry
- Physical reference frames
- Lie-theoretic operations

By overloading algebraic operations, it creates a computational syntax that mirrors mathematical intuition while remaining efficient for computer implementation. This approach bridges the gap between abstract mathematics and practical computation, particularly in fields requiring rigorous treatment of coordinate systems and their transformations.

##  Paper online
https://zenodo.org/records/14435614
  
## Code Compilation and Usage
Regarding the compilation and usage of these codes, due to some codes being related to the company's confidentiality policy, I can only release a part of the codes. However, the key points are transparent. You can combine these coordinate system codes with your own vector library for use, or directly use the Python version (currently, it only supports the Windows version). I hope this can be helpful and inspiring to you.
