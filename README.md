# The Coordinate System
The Coord System is a mathematical object based on group theory, which defines a coordinate system in three-dimensional space. It includes an origin, three directional axes, and three scaling components that correspond to translation, rotation, and scaling transformations. Coord provides methods to construct coordinate systems, such as constructing them from three axes or Euler angles. Additionally, it offers multiplication and division operations to transform vectors from one coordinate system to another and to project vectors from a parent coordinate system to a local one. This document describes how Coord can be applied in common scenarios, such as converting between world and local coordinate systems and using it in multi-node hierarchies.

## Structure of the Coordinate System

In C++, a coordinate system in three-dimensional space is defined by an origin, three directional axes, and three scaling components as follows:

```
struct coord {
    vec3 ux, uy, uz;   // Three unit vectors
    vec3 scale;        // Scaling
    vec3 o;            // Origin
};
```

## Constructing a Coordinate System

A coordinate system can be constructed using three axes or Euler angles as follows:

```
coord C(vec3 ux, vec3 uy, vec3 uz);
coord C(vec3 ux, vec3 uy);
coord C(float angle, vec3 axis); 
coord C(float pitch, float yaw, float roll);
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
V2 = V5 * C5 * C4 * C3 
V3 = V2 / C3 / C4 / C5
```

4. Convert between parallel coordinate systems:

```
C0{ 
    C1, C2; 
}
V2 = V1 * C1 / C2
```

5. More operations:

Scalar multiplication:

```
C*k = C.scale * k
```

Quaternion multiplication:

```
C2 = C1 * q 
C1 = C2 / q
q2 = q1 * C
q1 = q2 / C
```

Vector addition:

```
C1 + o = C1’{C1.o + o, C1.v + C2.v}.normalized
```

## Coordinate System Differentiation

Coord can be used to differentiate coordinate systems in space in three ways:

Gradient: 

```
▽f = dF(U * df * Cuv) / dxyz
```

Divergence:

```
▽ ∙ F = dF / dxyz ∙ Ic
```

Curl:

```
▽ x F = dF / dxyz x Ic
```

## Applications in Differential Geometry

Coord can be used to transport vectors in a curved space from a natural coordinate system to a curved coordinate system. The coordinate system can be represented in exponential form as C = e^(V). The curvature can be calculated by comparing two paths projected onto the u and v curves, using Gu and Gv.

The spatial curvature can be computed using a coordinate system, and the Riemann curvature tensor in the u,v coordinate system is given by:

```
Ruv = Gu*Gv - Gv*Gu - G[uv]

where: Gu = UG - ONE
       UG = C2 / C1
       Connection vector: [U, V] (Lie bracket operation)
```
## Combination with Lie Groups and Lie Algebras

Coord can be combined with Lie groups and Lie algebras. The rotation matrix R is an element of the Lie group, and the multiplication operation of the coordinate system is equivalent to the rotation matrix. Therefore, the coordinate system C is an element of the Lie group, with multiplication as the operation and the identity element as ONE. The cross product of two vectors can be expressed using the Lie algebra bracket: [C1, C2] = C1 * C2 – C2 * C1.

## Summary

The text discusses the need for a more suitable language in the computer era to express mathematical problems in a way that is acceptable to both humans and computers. The solution proposed is to overload algebraic operations and redefine them in a coordinate system object, with the aim of achieving a more concise and computer-friendly mathematical language by extending algebraic expression.
