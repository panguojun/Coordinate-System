# Coordinate System Transformation

### *Introduction*
*Coordinate system transformation is often performed using a matrix, but as a mathematical object, the matrix is not customized for coordinate system changes, and the matrix is too mathematical and the meaning is ambiguous.*
*Tensors are too abstract not only difficult to grasp but also difficult to quantify by computer, and a concept specially designed for coordinate system transformation is required.*
*This article introduces a simple and easy-to-understand coordinate system object and its algorithm.*
## Definition
A coordinate system in three-dimensional space consists of an origin plus three orientation axes and three scaling components. Corresponding to the three transformations of displacement, rotation, and scaling, respectively.
We define a structure:
````
struct coord_t
{
    vec3 o; // define the origin position
    vec3 ux,uy,uz; // three axial unit vectors
    vec3 scale; // scale transformation
}
````
*Note that the position, rotation and scaling of the coordinate system are all defined under its parent coordinate system.*
## Multiplication operation
```
// C1*C2*C3* ... *Cloc * Vloc （transfrom)
vec3 operator * (crvec v)
{
    return ux * v.x + uy * v.y + uz * v.z + o;
}
// V * C1 * C2 ...
friend vec3 operator * (crvec v, const coord_t& c)
{
    return c.ux * v.x + c.uy * v.y + c.uz * v.z + c.o;
}
coord_t operator * (coord_t& c)
{
    coord_t rc;
    rc.ux = ux * c.ux.x + uy * c.ux.y + uz * c.ux.z;
    rc.uy = ux * c.uy.x + uy * c.uy.y + uz * c.uy.z;
    rc.ux = ux * c.uz.x + uy * c.uz.y + uz * c.uz.z;
    rc.o += ux * c.o.x + uy * c.o.y + uz * c.o.z;
    return rc;
}
```
## Division operation
```
// Vworld/C1/C2/C3/ ... /Cloc（projection)
friend vec3 operator / (crvec v, const coord_t& c)
{
	vec3 dv = v - c.o;
	return vec3(dv.dot(c.ux), dv.dot(c.uy), dv.dot(c.uz));
}
coord_t operator / (const coord_t& c)
{
	coord_t rc;
	rc.ux = vec3(ux.dot(c.ux), ux.dot(c.uy), ux.dot(c.uz));
	rc.uy = vec3(uy.dot(c.ux), uy.dot(c.uy), uy.dot(c.uz));
	rc.uz = vec3(uz.dot(c.ux), uz.dot(c.uy), uz.dot(c.uz));
	rc.o -= c.o;
	return rc;
}
```
## Dot / Cross
```
vec3 dot(crvec v)
{
	vec3 dv = v - o;
	return vec3(dv.dot(ux) * scl.x, dv.dot(uy) * scl.y, dv.dot(uz) * scl.z);
}
vec3 cross(crvec v)
{
	vec3 dv = v - o;
	return dv.cross(ux) * scl.x + dv.cross(uy) * scl.y + dv.cross(uz) * scl.z;
}
```

## Sample 1: Centripetal force
### Polar coordinate transformation
````
void polecoord(coord_t& c1, real x, real y)
{
    c1.o.x = x;
    c1.o.y = y;
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
coord_t c1;
polecoord(c1, r, ang * PI / 180);
````
### Observe a coordinate system C2 at a point after dt next to C1
````
coord_t c2;
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
