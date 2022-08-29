vec3 mul_coord(const vec3& v)
{
    vec3 cv; 
    cv.x = x;
    cv.y = 2.0 * y;
    return cv;
}
vec3 div_coord(const vec3& cv)
{
    vec3 v; 
    v.x = cx;
    v.y = cy / 2.0;
    return v;
}
void curvature()
{
	coord3 c;
	vec3 v = vec3(1,0,0);
	vec3 cv = v * c;
	
}
