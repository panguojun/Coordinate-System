# Coordinate system transformation

*Introduction:*
*Coordinate system transformation is often performed using a matrix, but as a mathematical object, the matrix is not customized for coordinate system changes, and the matrix is too mathematical and the meaning is ambiguous. *
*Tensors are too abstract not only difficult to grasp but also difficult to quantify by computer, and a concept specially designed for coordinate system transformation is required. *
*This article introduces a simple and easy-to-understand coordinate system object and its algorithm*
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

## Sample:
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
### Sports settings
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
