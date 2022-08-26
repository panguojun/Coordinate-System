// =========================================================
// Electromagnetic Field And 3D Image
// =========================================================
#define resolution_3d	0.02 
#define depth_of_field	100
#define draw_strength_E	0.002
#define draw_strength_B	0.004

// ---------------------------------------------------------
// Electromagnetic Field:
// I define an electromagnetic force field in the Y-axis direction, 
// its scalar part exhibits an inverse square distribution near the Zero Point.
// Note: This field does not exist in reality!

// I prefer the equations in Maxwell's original form:
// Q = Fai + |A>
// ---------------------------------------------------------
// Eigen Electromagnetic Field Define

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
// E/B Field in default coordinate system
// ---------------------------------------------------------
vec3 get_Efield(vec3 p, real t)
{  
	return DIV(Fai, p, t) - DT(A, p, t); 
}
vec3 get_Bfield(vec3 p, real t)
{  
	return coord_cross(default_coord, CURL(A, p, t));
}

// ---------------------------------------------------------
// 3D Image
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
	
    float t = 0.5 * SEC;
    for(int i = 0; i < depth_of_field; i ++)
    {
		cor1 += (get_Efield(p - vec3(-0.25,0.,0.), t)-0.5*get_Efield(p - vec3(0.25,0.,0.), t))*draw_strength_B;
        
        p += v * dt;
    }
    
    fragColor = vec4(
		sin(cor1),	// MAKE it pretty!
		1.0);
}