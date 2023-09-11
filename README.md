# The Coordinate System
The Coord System is a mathematical object based on group theory, which defines a coordinate system in three-dimensional space. It includes an origin, three directional axes, and three scaling components that correspond to translation, rotation, and scaling transformations. Coord provides methods to construct coordinate systems, such as constructing them from three axes or Euler angles. Additionally, it offers multiplication and division operations to transform vectors from one coordinate system to another and to project vectors from a parent coordinate system to a local one. This document describes how Coord can be applied in common scenarios, such as converting between world and local coordinate systems and using it in multi-node hierarchies.

## Structure of the Coordinate System

In C++, a coordinate system in three-dimensional space is defined by an origin, three directional axes, and three scaling components as follows:

```
struct coord {
    vec3 ux, uy, uz;   // Three unit vectors
    vec3 s;            // Scaling
    vec3 o;            // Origin
};
```

## Constructing a Coordinate System

A coordinate system can be constructed using three axes or Euler angles as follows:

```
coord C1(vec3 ux, vec3 uy, vec3 uz);
coord C2(float angle, vec3 axis); 
coord C3(float pitch, float yaw, float roll);
coord C4(vec3 o, vec3 s, quat q); 
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
C2 = C1 + o;
where C2 = {C1.o + o, C1.v}, C1.v = {C1.ux*C1.s.x, C1.uy*C1.s.y, C1.uz*C1.s.z}
```

Coordinate Gradient:
```
G = C1 / C2 - I
where C1 and C2 are coordinate systems on two points on a unit length distance.
```
## Coordinate System Differentiation

Coord can be used to differentiate coordinate systems in space in three ways:

Gradient: 

```
▽f = (U * df * Cuv) / dxyz
where: dxyz = {ux * dx, uy * dy, uz * dz}
```

Divergence:

```
▽ ∙ F = dF / dxyz ∙ Ic
where: Ic = {ux, uy, uz}
```

Curl:

```
▽ x F = dF / dxyz x Ic
```

## Applications in Differential Geometry

Coord can transport vectors from a natural coordinate system to a curved coordinate system in a curved space. The curvature can be determined by comparing two paths projected onto the u and v curves, which is done using Gu and Gv. Gu and Gv represent the gradients of rotational changes of vectors along the u and v. By using a coordinate system, the spatial curvature can be calculated, and the Riemann curvature tensor in the u,v coordinate system is given by:

```
Ruv = Gu*Gv - Gv*Gu - G[u,v]

where:  Gu = C2 / C1 - I
        Connection vector: [u, v] (Lie bracket operation)
        W = Wu + Wv = [u, v]
        G[u,v] = Gu^Wu * Gv^Wv
```
## Combination with Lie Groups and Lie Algebras

Coord can be combined with Lie groups and Lie algebras. The rotation matrix R is an element of the Lie group, and the multiplication operation of the coordinate system is equivalent to the rotation matrix. Therefore, the coordinate system C is an element of the Lie group, with multiplication as the operation and the identity element as ONE. The cross product of two vectors can be expressed using the Lie algebra bracket:  
```
[C1, C2] = C1 * C2 – C2 * C1
```
By using coordinate system objects to achieve a unified formal representation of operations such as rotation, translation, scaling, and differential gradients, Lie groups and Lie algebras can be operated on in a unified manner.

## Summary

The text discusses the need for a more suitable language in the computer era to express mathematical problems in a way that is acceptable to both humans and computers. The solution proposed is to overload algebraic operations and redefine them in a coordinate system object, with the aim of achieving a more concise and computer-friendly mathematical language by extending algebraic expression.
