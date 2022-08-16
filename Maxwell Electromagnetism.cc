void Maxwell_Electromagnetism()
{
	coord_t c;
	
	auto A = [](crvec p, real t)->vec3 {
		vec3 dp = p;
		return dp * sin(t); 
	};
	auto DXYZ_A = [A](crvec p, real t)->vec3 {
		real d = 0.001f;
		return 
			(A(p + vec3::UX * d, t) - A(p, t)) / d +
			(A(p + vec3::UY * d, t) - A(p, t)) / d +
			(A(p + vec3::UZ * d, t) - A(p, t)) / d;
	};
	auto DT_A = [A](crvec p, real t)->vec3 {real dt = 0.001f; return (A(p,t + dt) - A(p, t)) / dt; };

	auto Fai = [](crvec p)->real {
		return real(1.0f / p.len()); 
	};
	auto DXYZ_Fai = [Fai](crvec p)->vec3 {
		real d = 0.001f; 
		return vec3(
			(Fai(p + vec3::UX * d) - Fai(p)) / d,
			(Fai(p + vec3::UY * d) - Fai(p)) / d,
			(Fai(p + vec3::UZ * d) - Fai(p)) / d);
	};
	{
		vec3 o = vec3::UX;
		vec3 deta = vec3::UX * 0.0;
		real t = 1;
		vec3 B = c.cross(DXYZ_A(o, t));
		vec3 E = -(DXYZ_Fai(o + deta) / c) - DT_A(o, t);

		PRINTVEC3(B);
		PRINTVEC3(E);
	}
	{
		vec3 o = vec3::UX;
		vec3 deta = vec3::UX * 0.5;
		real t = 1;
		vec3 B = c.cross(DXYZ_A(o, t));
		vec3 E = -(DXYZ_Fai(o + deta) / c) - DT_A(o, t);

		PRINTVEC3(B);
		PRINTVEC3(E);
	}
}
