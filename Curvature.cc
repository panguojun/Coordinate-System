void coord_at(coord3& c, vec3 q)
{
	q.y /= 2.0;
	c.ux = vec3(cos(q.y), sin(q.y), 0);
	c.uy = vec3(-q.x * sin(q.y), q.x * cos(q.y), 0);
}
void curvature()
{
	vec3 q = vec3(2, 1, 0);
	vec3 q11 = q;
	vec3 q21 = q + vec3(deta_d, 0, 0);
	vec3 q12 = q + vec3(0, deta_d, 0);
	vec3 q22 = q + vec3(deta_d, deta_d, 0);
	
	coord3 c11;
	coord_at(c11, q11);
	vec3 p11 = q11 * c11;

	coord3 c21;
	coord_at(c21, q21);
	vec3 p21 = q21 * c21;

	coord3 c12;
	coord_at(c12, q12);
	vec3 p12= q12 * c12;

	coord3 grad1 = c21 / c11; grad1.dump();
	coord3 grad2 = c12 / c11; grad2.dump();

	vec3 v(1,1,0);
	vec3 deta = v * grad1 * grad2 - v * grad2 * grad1;
	PRINTVEC3(v * grad1); PRINTVEC3(v * grad2);
	PRINTVEC3(deta);
}
