void curvature()
{
	coord c;
	
	auto A = [](crvec p, real t)->vec3 {
		vec3 dp = p;
		return dp * sin(t); 
	};
	
}
