# Coordinate System Transformation

### *Introduction*
*Coordinate system transformation is often performed using a matrix, but as a mathematical object, the matrix is not customized for coordinate system changes, and the matrix is too mathematical and the meaning is ambiguous.*
*Tensors are too abstract not only difficult to grasp but also difficult to quantify by computer, and a concept specially designed for coordinate system transformation is required.*
*This article introduces a simple and easy-to-understand coordinate system object and its algorithm.*
## Definition
A coordinate system in three-dimensional space consists of an origin plus three orientation axes and three scaling components. Corresponding to the three transformations of displacement, rotation, and scaling, respectively.
We define a structure:
````
struct coord
{
    vec3 o; 		// define the origin position
    vec3 ux,uy,uz; 	// three axial unit vectors
    vec3 scale; 	// scale transformation
}
````
*Note that the position, rotation and scaling of the coordinate system are all defined under its parent coordinate system.*
## Multiplication operation: define a vector in a coordinate system.
```
vec3 operator * (crvec p)
{
    return ux * v.x + uy * v.y + uz * v.z + o;
}
friend vec3 operator * (crvec p, const coord& c)
{
    return c.ux * p.x + c.uy * p.y + c.uz * p.z + c.o;
}
coord operator * (coord& c)
{
    coord rc;
    rc.ux = ux * c.ux.x + uy * c.ux.y + uz * c.ux.z;
    rc.uy = ux * c.uy.x + uy * c.uy.y + uz * c.uy.z;
    rc.uz = ux * c.uz.x + uy * c.uz.y + uz * c.uz.z;
    rc.o = o + ux * c.o.x + uy * c.o.y + uz * c.o.z;
    return rc;
}
```
## Division operation: measure this vector in a coordinate system.
```
friend vec3 operator / (crvec p, const coord& c)
{
	vec3 v = p - c.o;
	return vec3(v.dot(c.ux), v.dot(c.uy), v.dot(c.uz));
}
coord operator / (const coord& c)
{
	coord_t rc;
	rc.ux = vec3(ux.dot(c.ux) / c.scl.x, ux.dot(c.uy) / c.scl.y, ux.dot(c.uz) / c.scl.z);
	rc.uy = vec3(uy.dot(c.ux) / c.scl.x, uy.dot(c.uy) / c.scl.y, uy.dot(c.uz) / c.scl.z);
	rc.uz = vec3(uz.dot(c.ux) / c.scl.x, uz.dot(c.uy) / c.scl.y, uz.dot(c.uz) / c.scl.z);
	rc.o -= c.o;
	return rc;
}
```
## Dot / Cross
```
vec3 dot(crvec v)
{
	return vec3(v.dot(ux) * scl.x, v.dot(uy) * scl.y, v.dot(uz) * scl.z);
}
vec3 cross(const coord3& c)
{
	return vec3(
		uy.dot(c.UZ()) - uz.dot(c.UY()) +
		uz.dot(c.UX()) - ux.dot(c.UZ()) +
		ux.dot(c.UY()) - uy.dot(c.UX())
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
real w = 2.0f; // angular velocity
real dt = 0.001; // deta time
real dth = dt * w; // deta theta
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
vec3 dv = (c2* v - c1 * v) ; // Instead, observe the velocity vector in the global coordinate system and differentiate the velocity
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
## Sample 2: Curvature
### Curvature, I'm not quite sure if this curvature algorithm is correct, if so it can greatly simplify tensor calculations.
``
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
	return  grad1 * grad2 - grad2 * grad1;
}
``
