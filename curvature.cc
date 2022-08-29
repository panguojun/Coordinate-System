void coord(vec3& cv, vec3& v)
{
    cv.x = x;
    cv.y = 2.0 * y;
}

void curvature()
{
	coord c;
	auto A = [](crvec p, real t)->vec3 {
		vec3 dp = p;
		return dp * sin(t); 
	};
	
}
