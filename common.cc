// Using procedural algorithms should not be a numerical tool of mathematics,
// Program algorithms should be calculated in it's own form.

/** I think the physical formula is made up of several parts:
			
			Gauge * EigenNumbers * Dimension

 1) Coordinate System is a simple form of Gauge.
 2) Quaternions are typical EigenNumbers.
 3) Dimension is related to measurement and is the bridge between mathematics and physics.
*/

#define real float
#define ROT(p, a) p=cos(a)*p+sin(a)*vec2(p.y, -p.x)
#define deta_d 0.0001
#define deta_t 0.01

// ---------------------------------------------------------
// Coordinate System:
// A coordinate system in three-dimensional space 
// consists of an origin plus three orientation axes 
// ---------------------------------------------------------
struct coord3
{
   vec3 ux,uy,uz; // three axial unit vectors
};
coord3 create_coord(vec3 _ux, vec3 _uy, vec3 _uz)
{
    coord3 c;c.ux = _ux;c.uy = _uy; c.uz = _uz;
    return c;
}
// div: measure this vector in a coordinate system
vec3 coord_div (vec3 p, coord3 c)
{
    vec3 v = p;
    return vec3(dot(v,c.ux), dot(v,c.uy), dot(v,c.uz));
}
// mul: measure this vector in a coordinate system
vec3 coord_mul (vec3 p, coord3 c)
{
	return c.ux * (p.x) + c.uy * (p.y) + c.uz * (p.z);
}
// eigen vector
vec3 coord_eigenvec(coord3 c)
{
	return (c.ux + c.uy + c.uz);
}
// dot
real coord_dot(vec3 v, coord3 c)
{
	return dot(v, coord_eigenvec(c));
}
// cross
vec3 coord_cross(coord3 a,coord3 b)
{
    return vec3(
        dot(a.uy,b.uz) - dot(a.uz,b.uy),
        dot(a.uz,b.ux) - dot(a.ux,b.uz),
        dot(a.ux,b.uy) - dot(a.uy,b.ux)
    );
}

coord3 default_coord;		// Default Coordinate System

// ---------------------------------------------------------
// DIV, Curl
// ---------------------------------------------------------
#define DIV(Fai, p, t) \
        vec3((Fai(p + vec3(deta_d,0.0,0.0), t) - Fai(p, t)) / deta_d,\
        (Fai(p + vec3(0.0,deta_d,0.0), t) - Fai(p, t)) / deta_d, \
        (Fai(p + vec3(0.0,0.0,deta_d), t) - Fai(p, t)) / deta_d)

#define CURL(A, p, t) \
    create_coord( \
    (A(p + vec3(1.0,0.0,0.0) * deta_d, t) - A(p, t)) / deta_d, \
    (A(p + vec3(0.0,1.0,0.0) * deta_d, t) - A(p, t)) / deta_d, \
    (A(p + vec3(0.0,0.0,1.0) * deta_d, t) - A(p, t)) / deta_d)
    
#define DT(A, p, t) (A(p,t + deta_t) - A(p, t)) / deta_t


// ---------------------------------------------------------
// Quaternion:
// Quaternions are real mathematics numbers, 
// meaning there is number theory behind them.
// ---------------------------------------------------------
struct quaternion
{
   real w,x,y,z;
};
quaternion angax_q(real ang, vec3 ax)
{
    quaternion q;
    real halfang = 0.5 * ang;
    real fsin = sin(halfang);
    q.w = cos(halfang);
    q.x = fsin * ax.x;
    q.y = fsin * ax.y;
    q.z = fsin * ax.z;
    return q;
}
quaternion qmul (quaternion q1, quaternion q2)
{
	quaternion q;
    
	q.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
	q.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
	q.y = q1.w * q2.y + q1.y * q2.w + q1.z * q2.x - q1.x * q2.z;
	q.z = q1.w * q2.z + q1.z * q2.w + q1.x * q2.y - q1.y * q2.x;
	return q;
}
quaternion qrot(quaternion q, real ang, vec3 ax)
{
	quaternion qt = angax_q(ang, ax);
	return qmul(q, qt);
}

// ---------------------------------------------------------
// Dimension or Const:
// Dimensions and constants define our physical world!
// ---------------------------------------------------------
#define PI 3.1415926535
#define MET 1.0
#define SEC 0.5
