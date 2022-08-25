// =========================================================
// Extraordinary scientific theories must seem extraordinary
// =========================================================
// Electromagnetic Field:
// I define an electromagnetic force field in the Y-axis direction, 
// and its scalar part exhibits an inverse square distribution near the Zero Point
// Note: This field does not exist in reality!

// ---------------------------------------------------------
// I want to display the electromagnetic fields together, 
// with the colors representing the field directions.
// ---------------------------------------------------------
#define real float
#define resolution_3d	0.01 
#define depth_of_field	180
#define draw_strength	0.0001

// ---------------------------------------------------------
// Quaternion:
// Quaternions are real numbers, 
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
// Coordinate System:
// A coordinate system in three-dimensional space 
// consists of an origin plus three orientation axes 
// 
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
real coord_dot(vec3 v, coord3 c)
{
	return dot(v, coord_eigenvec(c));
}
vec3 coord_cross(coord3 a,coord3 b)
{
    return vec3(
        dot(a.uy,b.uz) - dot(a.uz,b.uy),
        dot(a.uz,b.ux) - dot(a.ux,b.uz),
        dot(a.ux,b.uy) - dot(a.uy,b.ux)
    );
}

// ---------------------------------------------------------
// Dimension or Const
// Dimensions and constants define our physical world!
// ---------------------------------------------------------
#define PI 3.1415926535
#define MET 1.0
#define SEC 1.0

// ---------------------------------------------------------
// Electromagnetic Feild
// I prefer the equations in Maxwell's original form:
// Q = Fai + |A>
// ---------------------------------------------------------
// The Scalar Part
float Fai(vec3 p, float t)
{
    float r = length(p) * 1.;
	return (1.0 / (r*r)) * (cos(t*0.25)); 
}
// The Vector Part
vec3 A(vec3 p, float t) {
    float r = length(p.xz) * 1.;
    return vec3(0.0,1.0,0.0) * (1.0 / (r * r)*sin(t*0.25)); 
}
// Scalar Divergence
vec3 DXYZ_Fai(vec3 p, float t)
{
    float d = 0.0001; 
    return vec3(
        (Fai(p + vec3(d,0.0,0.0), t) - Fai(p, t)) / d,
        (Fai(p + vec3(0.0,d,0.0), t) - Fai(p, t)) / d,
        (Fai(p + vec3(0.0,0.0,d), t) - Fai(p, t)) / d);
}
// Vectors Curl
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
    float dt = 0.0001;
    return (A(p,t + dt) - A(p, t)) / dt; 
}

// ---------------------------------------------------------
// 3D View
// ---------------------------------------------------------
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 pp = (-iResolution.xy + 2.0*fragCoord.xy) / iResolution.y;
    float eyer = 1.0;
    float eyea = -((iMouse.x) / iResolution.x) * PI * 2.0;
    float eyef = ((iMouse.y / iResolution.y)-0.24) * PI * 2.0;
    
	vec3 cam = vec3(
        eyer * cos(eyea) * sin(eyef),
        eyer * cos(eyef),
        eyer * sin(eyea) * sin(eyef));
    
	vec3 front = normalize(- cam);
	vec3 left = normalize(cross(normalize(vec3(0.0,1,-0.01)), front));
	vec3 up = normalize(cross(front, left));
	vec3 v = normalize(front*1.5 + left*pp.x + up*pp.y);
    
    vec3 p = cam;
    
    coord3 c;
    c.ux = vec3(1.0,0.0,0.0);
    c.uy = vec3(0.0,1.0,0.0);
    c.uz = vec3(0.0,0.0,1.0);
    
    float dt = resolution_3d;
    vec3 cor1 = vec3(0.0);
    vec3 cor2 = vec3(0.0);
	
    float t = iTime * SEC;
    for(int i = 0; i < depth_of_field; i ++)
    {
        float r = length(p);
        {  
            vec3 o = p;
            vec3 B = coord_cross(c,DXYZ_A(o, t));
            vec3 E = (DXYZ_Fai(o, t)) - DT_A(o, t); 
            cor1 += (B)*draw_strength;
            cor2 += (E)*draw_strength;
        }
        p += v * dt;
    }
    
    fragColor = vec4(
		sin(cor1*4.0)+
		sin(cor2*4.0),
		1.0);
}