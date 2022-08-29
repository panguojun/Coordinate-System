void curvature()
{
	vec3 v = vec3(1,0,0);
	coord3 c11(v);
	coord3 c21(c11.o + vec3::UX * deta_d);
	coord3 c12(c11.o + vec3::UY * deta_d);
	coord3 c22_21(c21.o + vec3::UY * deta_d);
	coord3 c22_12(c12.o + vec3::UX * deta_d);
	vec3 v11 = v * c11;
	vec3 v12 = v * c12;
	vec3 v21 = v * c21;
	vec3 v22 = v * c22;
}
