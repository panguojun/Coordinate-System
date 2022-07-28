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
*Note that the position, rotation and scaling of the coordinate system are all defined under its parent coordinate system. *
