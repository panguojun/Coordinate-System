# Coordinate System

## Introduction
Coordinate system transformations are typically achieved using matrices. However, matrices are not specifically designed for coordinate system changes, and their mathematical nature can lead to ambiguity. 
Tensors, while more tailored to coordinate system transformations, are often too abstract and challenging for both human comprehension and computer quantification. Therefore, there is a need for a concept that is specifically designed for coordinate system transformations. 
This article presents a straightforward and easily understandable coordinate system object and its corresponding algorithm.
## Definition
A coordinate system in three-dimensional space consists of an origin plus three orientation axes and three scaling components. Corresponding to the three transformations of displacement, rotation, and scaling, respectively.
We define a structure:
````
struct coord
{
    vec3 ux,uy,uz; 	// three axial unit vectors
    vec3 o; 		// define the origin position
    vec3 s; 	        // scale transformation
}
````
*Note that the position, rotation and scaling of the coordinate system are all defined under its parent coordinate system.*
## Multiplication operation: define a vector in a coordinate system.
```
friend vec3 operator * (crvec p, const coord& c)
{
    return c.ux * p.x + c.uy * p.y + c.uz * p.z;
}
coord operator * (coord& c)
{
    coord rc;
    rc.ux = ux.x * c.ux + uy.x * c.uy + uz.x * c.uz;
    rc.uy = ux.y * c.ux + uy.y * c.uy + uz.y * c.uz;
    rc.uz = ux.z * c.ux + uy.z * c.uy + uz.z * c.uz;
    return rc;
}
```
## Division operation: measure this vector in a coordinate system.
```
friend vec3 operator / (crvec v, const coord& c)
{
	return vec3(v.dot(c.ux), v.dot(c.uy), v.dot(c.uz));
}
coord operator / (const coord& c)
{
	coord_t rc;
	rc.ux = vec3(ux.dot(c.ux), ux.dot(c.uy), ux.dot(c.uz));
	rc.uy = vec3(uy.dot(c.ux), uy.dot(c.uy), uy.dot(c.uz));
	rc.uz = vec3(uz.dot(c.ux), uz.dot(c.uy), uz.dot(c.uz));
	return rc;
}
```
## Dot / Cross
```
real dot(crvec v)
{
	return v.dot(ux) + v.dot(uy) + v.dot(uz);
}
coord3 cross(const coord3& c)
{
	return coord3(
		vec3::UX * (uy.dot(c.uz) - uz.dot(c.uy)),
		vec3::UY * (uz.dot(c.ux) - ux.dot(c.uz)),
		vec3::UZ * (ux.dot(c.uy) - uy.dot(c.ux))
	);
}
```

## Sample 1: Centripetal force
### Polar coordinate transformation
````
void polecoord(coord& c1, real x, real y)
{
    c1.ux = vec3(cos(c1.o.y), sin(c1.o.y), 0);
    c1.uy = vec3(-c1.o.x * sin(c1.o.y), c1.o.x * cos(c1.o.y), 0);
}
````
### Kinematics settings
````
real w = 2.0f; 		// angular velocity
real dt = 0.001; 	// deta time
real dth = dt * w; 	// deta theta
````

### Observe a coordinate system C1 at point (r, ang) in global coordinates
````
real r = 1, ang = 180;
coord c1;
polecoord(c1, r, ang * PI / 180);
````
### Observe a coordinate system C2 at a point after dt next to C1
````
coord c2;
polecoord(c2, r, dth + ang * PI / 180);
````
### Velocity vector in polar coordinates (circular motion)
````
vec3 v(0.0, w, 0);
vec3 dv = (v * c2 - v * c1) ; // Instead, observe the velocity vector in the global coordinate system and differentiate the velocity
vec3 a = dv / dt; // the time derivative of the vector, the acceleration
````

### Print the result, which is consistent with the formula of w*w*r
````
PRINTVEC3(a); // wwr
````

## Sample 2: Maxwell Electromagnetism
### Electromagnetic potential
```
auto A = [](crvec p, real t)->vec3 {
	vec3 dp = p;
	return dp * sin(t); 
};
auto DXYZ_A = [A](crvec p, real t)->coord3 {
		coord3 c;
		real d = 0.001f;
		c.ux = ((A(p + vec3::UX * d, t) - A(p, t)) / d);
		c.uy = ((A(p + vec3::UY * d, t) - A(p, t)) / d);
		c.uz = ((A(p + vec3::UZ * d, t) - A(p, t)) / d);
		return c;
	};
auto DT_A = [A](crvec p, real t)->vec3 {real dt = 0.001f; return (A(p,t + dt) - A(p, t)) / dt; };
```
### Phase
```
auto Fai = [](crvec p)->real {
	return real(1.0f / p.len()); 
};
auto DXYZ_Fai = [Fai](crvec p)->vec3 {
	real d = 0.001f; 
	return vec3(
		(Fai(p + vec3::UX * d) - Fai(p)) / d,
		(Fai(p + vec3::UY * d) - Fai(p)) / d,
		(Fai(p + vec3::UZ * d) - Fai(p)) / d);
};
```
### Electromagnetic Field
```
{
	vec3 o = vec3::UX;
	real t = 1;
	vec3 B = c.cross(DXYZ_A(o, t));
	vec3 E = -(DXYZ_Fai(o) / c) - DT_A(o, t);

	PRINTVEC3(B);
	PRINTVEC3(E);
}
```
## Sample 3: Curvature
Ruv = Gu*Gv - Gv*Gu * Gu^Wu * Gv^Wv;  W = (U + V*Gu) - (V + U*Gv)
```
coord3 curvature()
{
	vec3 q = vec3(2, 1, 0);
	vec3 q11 = q;
	vec3 q21 = q + vec3(deta_d, 0, 0);
	vec3 q12 = q + vec3(0, deta_d, 0);
	vec3 q22 = q + vec3(deta_d, deta_d, 0);
	
	coord3 c11;
	coord_at(c11, q11);
	vec3 p11 = q11 * c11;

	coord3 c21;
	coord_at(c21, q21);
	vec3 p21 = q21 * c21;

	coord3 c12;
	coord_at(c12, q12);
	vec3 p12= q12 * c12;

	coord3 grad1 = c21 / c11; grad1.norm(false);
	coord3 grad2 = c12 / c11; grad2.norm(false);

	vec3 v(1,1,0);
	vec3 deta = v * grad1 * grad2 - v * grad2 * grad1;
	PRINTVEC3(v * grad1); PRINTVEC3(v * grad2);
	PRINTVEC3(deta);
	return  grad1 * grad2 - grad2 * grad1; // just related to curvature
}
```
