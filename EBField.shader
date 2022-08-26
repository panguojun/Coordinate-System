// Using procedural algorithms should not be a numerical tool of mathematics,
// Program algorithms should be calculated in it's own form.

/** I think the physical formula is made up of several parts:
			
			Gauge * EigenNumbers * Dimension

 1) Coordinate System is a simple form of Gauge.
 2) Quaternions are typical EigenNumbers.
 3) Dimension is related to measurement and is the bridge between mathematics and physics.
*/
// ---------------------------------------------------------
// Electromagnetic Field:
// I define an electromagnetic force field in the Y-axis direction, 
// its scalar part exhibits an inverse square distribution near the Zero Point.
// Note: This field does not exist in reality!

// ---------------------------------------------------------
// I want to display the electromagnetic fields together, 
// with the colors representing the field directions.
// ---------------------------------------------------------
#define real float
#define resolution_3d	0.02 
#define depth_of_field	100
#define draw_strength_E	0.002
#define draw_strength_B	0.004

#define ROT(p, a) p=cos(a)*p+sin(a)*vec2(p.y, -p.x)

// ---------------------------------------------------------
// Coordinate System:
// A coordinate system in three-dimensional space 
// consists of an origin plus three orientation axes 
// ---------------------------------------------------------
struct coord3
{
   vec3 ux,uy,uz; // three axial unit vectors
};
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

// =========================================================
// Electromagnetic Feild:
// I prefer the equations in Maxwell's original form:
// Q = Fai + |A>
// =========================================================
// The Scalar Part
float Fai(vec3 p, float t)
{
    float r = length(p);
    if(r < 0.001)
        r = 0.001;
	return (1.0 / (r*r)) * (cos(t*0.25)); 
}
// The Vector Part
vec3 A(vec3 p, float t) {
    float r = length(p.xz);
    if(r < 0.001)
        r = 0.001;
    return vec3(0.0,1.0,0.0) * (10.0*sin(t*0.25) / (r)); 
}

// ---------------------------------------------------------
// Scalar Divergence
vec3 DXYZ_Fai(vec3 p, float t)
{
    float d = 0.0001; 
    return vec3(
        (Fai(p + vec3(d,0.0,0.0), t) - Fai(p, t)) / d,
        (Fai(p + vec3(0.0,d,0.0), t) - Fai(p, t)) / d,
        (Fai(p + vec3(0.0,0.0,d), t) - Fai(p, t)) / d);
}
// Vector Curl
coord3 DXYZ_A(vec3 p, float t){
    coord3 c;
    float d = 0.0001;
    c.ux = ((A(p + vec3(1.0,0.0,0.0) * d, t) - A(p, t)) / d);
    c.uy = ((A(p + vec3(0.0,1.0,0.0) * d, t) - A(p, t)) / d);
    c.uz = ((A(p + vec3(0.0,0.0,1.0) * d, t) - A(p, t)) / d);
    return c;
}
// Time
vec3 DT_A(vec3 p, float t)
{
    float dt = 0.01;
    return (A(p,t + dt) - A(p, t)) / dt; 
}

coord3 default_coord;		// Default Coordinate System

vec3 get_Efield(vec3 p, real t)
{  
	return (DXYZ_Fai(p, t)) - DT_A(p, t); 
}
vec3 get_Bfield(vec3 p, real t)
{  
	return coord_cross(default_coord,DXYZ_A(p, t));
}

// ---------------------------------------------------------
// View
// ---------------------------------------------------------
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 pp = (-iResolution.xy + 2.0 * fragCoord.xy) / iResolution.y;
    float eyer = 0.5;
    float eyea = -((iMouse.x) / iResolution.x) * PI * 2.0;
    float eyef = ((iMouse.y / iResolution.y) - 0.24) * PI * 2.0;
    
	vec3 cam = vec3(
        eyer * cos(eyea) * sin(eyef),
        eyer * cos(eyef),
        eyer * sin(eyea) * sin(eyef));
         
    ROT(cam.xz, (0.25) * (iTime + 15.0)); // auto rotation
    
	vec3 front = normalize(- cam);
	vec3 left = normalize(cross(normalize(vec3(0.0,1,-0.001)), front));
	vec3 up = normalize(cross(front, left));
	vec3 v = normalize(front + left*pp.x + up*pp.y);
   
    vec3 p = cam;
	
	// Default Coordinate System
    default_coord.ux = vec3(1.0,0.0,0.0);
    default_coord.uy = vec3(0.0,1.0,0.0);
    default_coord.uz = vec3(0.0,0.0,1.0);
    
    float dt = resolution_3d;
    vec3 cor1 = vec3(0.0);
    vec3 cor2 = vec3(0.0);
	
    float t = iTime * SEC;
    for(int i = 0; i < depth_of_field; i ++)
    {
		cor1 += get_Bfield(p, t)*draw_strength_B;
		cor2 += get_Efield(p, t)*draw_strength_E;
        
        p += v * dt;
    }
    
    fragColor = vec4(
		sin(cor1+ cor2),	// MAKE it pretty!
		1.0);
}
