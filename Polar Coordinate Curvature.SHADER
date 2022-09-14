
/*
                          .-''--.
                         _`>   `\.-'<
                      _.'     _     '._
                    .'   _.='   '=._   '.
                    >_   / /_\ /_\ \   _<
                      / (  \o/\\o/  ) \
                      >._\ .-,_)-. /_.<
                          /__/ \__\
                            '---'  
            
          “If I can’t picture it, I can’t understand it.”
*/
// Changes in the angular extent of the polar coordinate system produce the curvature of the pole.
// Instead of tensors, I use a coordinate system of my own design to calculate the curvature.
// Although it is likely to fail(99.9% to fail :-D), but if it succeeds, tensor operations can be greatly simplified, 
// and may allow us to calculate General Relativity directly!
// So why not have a try?

#define real float
#define ROT(p, a) p=cos(a)*p+sin(a)*vec2(p.y, -p.x)
#define pi 3.1415926535
// ---------------------------------------------------------
// Coordinate System:
// A coordinate system in two-dimensional space 
// consists of an origin plus two orientation axes 
// ---------------------------------------------------------
struct coord2
{
   vec2 ux,uy;   // two axial unit vectors
   float s;      // scaling[opt]
   vec2 o;       // origin
};
// mul: define a vector in a coordinate system or transform a vector
vec2 coord_mul (vec2 p, coord2 c)
{
	return c.ux * (p.x*c.s) + c.uy * (p.y*c.s) + c.o;
}
// 
coord2 coord_mul2 (coord2 c1, coord2 c2)
{
    coord2 rc;
    rc.ux = c1.ux * c2.ux.x + c1.uy * c2.ux.y;
    rc.uy = c1.ux * c2.uy.x + c1.uy * c2.uy.y;
    rc.s = c1.s * c2.s;
    //rc.o = c1.o + c1.ux * c2.o.x + c1.uy * c2.o.y;
    return rc;
}
// div: measure a vector in a coordinate system
vec2 coord_div (vec2 p, coord2 c)
{
    vec2 v = p - c.o;
    return vec2(dot(v,normalize(c.ux))/c.s, dot(v,normalize(c.uy))/c.s);
}
// grad: I'm not sure if gradient is the right word
coord2 coord_grad (coord2 c1, coord2 c2)
{
    coord2 rc;
    rc.ux = vec2(dot(c1.ux, normalize(c2.ux)), dot(c1.ux, normalize(c2.uy)));
    rc.uy = vec2(dot(c1.uy, normalize(c2.ux)), dot(c1.uy, normalize(c2.uy)));
    rc.o = c1.o - c2.o;
    rc.s = c1.s/c2.s;
    return rc;
}
// sub
coord2 coord_sub (coord2 c1, coord2 c2)
{
    coord2 rc;
    rc.ux = c1.ux - c2.ux;
    rc.uy = c1.uy - c2.uy;
    rc.o = c1.o - c2.o;
    rc.s = 1.;
    return rc;
}
// Eigen vector just like matrix
vec2 eigenvec(coord2 c)
{
    return (c.ux + c.uy) * c.s;
}

// Polar Coordinate Differential
void polar_coord(out coord2 c1, float x, float y)
{
    //y /= 2.0*pi; // change the angle
	c1.ux = vec2((cos(y)), (sin(y)));
	c1.uy = vec2(x*(-sin(y)), x*cos(y));
    c1.s = 1.;
    c1.o = vec2(0.0,0.0);
}
// ---------------------------------------
// mainImage
// ---------------------------------------
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (fragCoord/iResolution.xy - vec2(0.5,0.5)) * 2.0;

    float deta = 0.01; // delta
    coord2 c11;
    polar_coord(c11, uv.x, uv.y); // uv.x = r, uv.y = theta
    
    coord2 c21;
    polar_coord(c21, uv.x + deta, uv.y);
    coord2 grad12 = coord_grad(c21, c11);
    
    coord2 c12;
    polar_coord(c12, uv.x, uv.y + deta);
    coord2 grad21 = coord_grad(c12, c11);
    
    coord2 g12 = coord_mul2(grad12,grad21);
    coord2 g21 = coord_mul2(grad21,grad12);
    
    coord2 R = coord_sub(g12, g21); // related to curvature
    
    vec3 col = vec3(1.-cos(eigenvec(R)/(deta)), 0.0);// green for curvature.
    
    // Show Polar Coordinate Differential in the background.
    vec2 deta_polar_coord = sin(coord_mul(uv, c11) * 40.0);//sin(coord_div(uv, c11) * 40.0);
    col += vec3(deta_polar_coord.x, 0.0, deta_polar_coord.y)*0.2;
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
