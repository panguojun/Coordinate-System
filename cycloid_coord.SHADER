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
// The cycloid is related to the fastest curve, 
// but it is not convenient for us to calculate the integral problem here.

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
// mul: define a vector in a coordinate system
vec2 coord_mul (vec2 p, coord2 c)
{
	return c.ux * (p.x*c.s) + c.uy * (p.y*c.s) + c.o;
}
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
coord2 coord_grad (coord2 c1, coord2 c2)
{
    coord2 rc;
    rc.ux = vec2(dot(c1.ux, normalize(c2.ux)), dot(c1.ux, normalize(c2.uy)));
    rc.uy = vec2(dot(c1.uy, normalize(c2.ux)), dot(c1.uy, normalize(c2.uy)));
    rc.o = c1.o - c2.o;
    rc.s = c1.s/c2.s;
    return rc;
}
coord2 coord_sub (coord2 c1, coord2 c2)
{
    coord2 rc;
    rc.ux = c1.ux - c2.ux;
    rc.uy = c1.uy - c2.uy;
    rc.o = c1.o - c2.o;
    rc.s = 1.;
    return rc;
}
vec2 eigenvec(coord2 c)
{
    return (c.ux + c.uy) * c.s;
}
// cycloid coordinate differential
void cycloid_coord(out coord2 c1, float x, float y)
{
    //y /= 2.0;
	c1.ux = vec2((cos(y)), (sin(y)));
	c1.uy = vec2(x*(-sin(y)), x*cos(y));
    c1.s = 1.;
    c1.o = vec2(0.5,0.5);
}
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (fragCoord/iResolution.xy - vec2(0.5,0.5)) * 2.0;
    
    float deta = 0.01;
    coord2 c11;
    cycloid_coord(c11, uv.x, uv.y);
    vec2 duv = coord_div(uv, c11);
    
    coord2 c21;
    cycloid_coord(c21, uv.x + deta, uv.y);
    coord2 grad12 = coord_grad(c21, c11);
    
    coord2 c12;
    cycloid_coord(c12, uv.x, uv.y + deta);
    coord2 grad21 = coord_grad(c12, c11);
    
    coord2 g12 = coord_mul2(grad12,grad21);
    coord2 g21 = coord_mul2(grad21,grad12);
    
    coord2 R = coord_sub(g12, g21);
    //eigenvec(c11)
    vec3 col = vec3(sin(duv/deta), 0.25);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
