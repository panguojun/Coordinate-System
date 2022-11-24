/**
*				融合函数
*				包括非线性融合，线性插值等
*/
// **********************************************************************
// shape blend
// **********************************************************************
int blendi(int h1, int h2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	if (power != 1.0)
		alpha = pow(alpha, power);

	return int(h1 * (1.0f - alpha) + h2 * alpha);
}
real blend(real h1, real h2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	if (power != 1.0)
		alpha = pow(alpha, power);

	return h1 * (1.0f - alpha) + h2 * alpha;
}
vector2 blend(const vector2& v1, const vector2& v2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	if (power != 1.0)
		alpha = pow(alpha, power);

	return v1 * (1 - alpha) + v2 * alpha;
}
vector3 blend(const vector3& v1, const vector3& v2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	if (power != 1.0)
		alpha = pow(alpha, power);

	return v1 * (1 - alpha) + v2 * alpha;
}
vectorn blend(const vectorn& v1, const vectorn& v2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	if (power != 1.0)
		alpha = pow(alpha, power);

	return v1 * (1.0f - alpha) + v2 * alpha;
}
// ------------------------------------------------
real blend2(real h1, real h2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;

	if (alpha < 0.5f)
	{
		alpha = (0.5f - alpha) * 2.0f;
		if (power != 1.0)
			alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;
	}
	else
	{
		alpha = (alpha - 0.5f) * 2.0f;
		if (power != 1.0)
			alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;
	}
}
// ------------------------------------------------
real blendn(real h1, real h2, real alpha, int n, real power)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;

	/*
	alpha = alpha * n - int(alpha * n);
	if (alpha < 0.5)
	{
		alpha = (0.5 - alpha) / 0.5;

	}
	else
	{
		alpha = (alpha - 0.5) / 0.5;
	}*/
	alpha = pow(fabs(cos(alpha * n * PI)), power);

	return h2 * (1 - alpha) + h1 * alpha;
}
// ------------------------------------------------
vector2 blend2(vector2 h1, vector2 h2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;

	if (alpha < 0.5f)
	{
		alpha = (0.5f - alpha) * 2.0f;
		if (power != 1.0)
			alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;
	}
	else
	{
		alpha = (alpha - 0.5f) * 2.0f;
		if (power != 1.0)
			alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;
	}
}
// ------------------------------------------------
vector3 blend2(vector3 h1, vector3 h2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;

	if (alpha < 0.5f)
	{
		alpha = (0.5f - alpha) * 2.0f;
		if (power != 1.0)
			alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;
	}
	else
	{
		alpha = (alpha - 0.5f) * 2.0f;
		if (power != 1.0)
			alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;
	}
}

real blend3(real h1, real h2, real alpha, real mid = 0.5, real power = 1)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;

	if (alpha < mid)
	{
		alpha = (mid - alpha) / mid;
		alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;
	}
	else
	{
		alpha = (alpha - mid) / (1 - mid);
		alpha = pow(alpha, power);
		return h2 * (1 - alpha) + h1 * alpha;
	}
}

real blend4(real h1, real h2, real alpha, real power = 1)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;

	real mid = 0.5;
	if (alpha < mid)
	{
		alpha = alpha / mid;
		alpha = pow(alpha, 1.0f / power);
		return h1 * (1 - alpha) + h2 * alpha;
	}
	else
	{
		alpha = (1 - alpha) / (1 - mid);
		alpha = pow(alpha, 1.0f / power);
		return h1 * (1 - alpha) + h2 * alpha;
	}
}

real blend5(real h1, real h2, real alpha, real mid = 0.5, real power = 1)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;

	if (alpha < mid)
	{
		alpha = alpha / mid;
		alpha = pow(alpha, 1.0f / power);
		return h1 * (1 - alpha) + h2 * alpha;
	}
	else
	{
		alpha = (1 - alpha) / (1 - mid);
		alpha = pow(alpha, 1.0f / power);
		return h1 * (1 - alpha) + h2 * alpha;
	}
}
inline real clamp(real v, real lo, real hi)
{
	ASSERT(!(hi < lo));
	return (v < lo) ? lo : (hi < v) ? hi : v;
}
inline int clampi(int v, int lo, int hi)
{
	ASSERT(!(hi < lo));
	return (v < lo) ? lo : (hi < v) ? hi : v;
}
// ------------------------------------------------------------------------------------------------
vec3 bezier2(vec3 cp[3], real t)
{
	return cp[0] * ((1 - t) * (1 - t)) + cp[1] * (2 * t * (1 - t)) + cp[2] * (t * t);
}
// ------------------------------------------------------------------------------------------------
vec3 bezier3(vec3 cp[4], real t)
{
	real s = 1 - t;
	return  cp[0] * (s * s * s) +
		cp[1] * (3 * t * s * s) +
		cp[2] * (3 * t * t * s) +
		cp[3] * (t * t * t);
}
// ------------------------------------------------------------------------------------------------
vec2 bezier2(vec2 cp[3], real t)
{
	return cp[0] * ((1 - t) * (1 - t)) + cp[1] * (2 * t * (1 - t)) + cp[2] * (t * t);
}
// ------------------------------------------------------------------------------------------------
vec2 bezier2(vec2 p1, vec2 p2, vec2 p3, real t)
{
	return p1 * ((1 - t) * (1 - t)) + p2 * (2 * t * (1 - t)) + p3 * (t * t);
}
// ------------------------------------------------------------------------------------------------
vec2 bezier3(vec2 cp[4], real t)
{
	real s = 1 - t;
	return  cp[0] * (s * s * s) +
		cp[1] * (3 * t * s * s) +
		cp[2] * (3 * t * t * s) +
		cp[3] * (t * t * t);
}
// ------------------------------------------------------------------------------------------------
vec bezier2(vec p1, vec p2, vec p3, real t)
{
	return p1 * ((1 - t) * (1 - t)) + p2 * (2 * t * (1 - t)) + p3 * (t * t);
}
// ------------------------------------------------------------------------------------------------
real bezier2(real cp[3], real t)
{
	return cp[0] * ((1 - t) * (1 - t)) + cp[1] * (2 * t * (1 - t)) + cp[2] * (t * t);
}
// ------------------------------------------------------------------------------------------------
real bezier3(real cp[4], real t)
{
	real s = 1 - t;
	return  cp[0] * (s * s * s) +
		cp[1] * (3 * t * s * s) +
		cp[2] * (3 * t * t * s) +
		cp[3] * (t * t * t);
}
// -----------------------------------------------------------------------
inline real roundblend(real h1, real h2, real alpha)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	//alpha = sqrt(1 - (1-alpha)*(1-alpha));
	alpha = sqrt(1 - (alpha) * (alpha));
	return h2 * (1.0f - alpha) + h1 * alpha;

}
// -----------------------------------------------------------------------
inline real roundblend2(real h1, real h2, real alpha, real power = 1)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;

	if (power != 1)
		alpha = pow(alpha, power);

	if (alpha < 0.5)
	{
		alpha = alpha * 2;
		alpha = sqrt(1 - (1 - alpha) * (1 - alpha));

		return h1 * (1 - alpha) + h2 * alpha;
	}
	else
	{
		alpha = (alpha - 0.5) * 2;
		alpha = sqrt(1 - (alpha) * (alpha));

		return h1 * (1 - alpha) + h2 * alpha;
	}
}
// -----------------------------------------------------------------------
// 三角函数插值
inline real BlendSin(real h1, real h2, real alpha)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;

	alpha = sin(alpha * PI / 2);

	return h1 * (1 - alpha) + h2 * alpha;
}
// -----------------------------------------------------------------------
// 傅立叶级数
inline real FT(real angle, real t[] = 0, real dt = 0)
{
	if (t == 0)
	{
		static real s_t0[] = { rrnd(0, 1), rrnd(0, 1), rrnd(0, 1), rrnd(0, 1) };
		t = s_t0;
	}

	real yy = 0;
	yy += 1 * sin(1 * angle + (t[0] + dt) * PI);
	yy += 0.5 * sin(2 * angle + (t[1] + dt) * PI);
	yy += 0.25 * sin(4 * angle + (t[2] + dt) * PI);
	yy += 0.125 * sin(8 * angle + (t[3] + dt) * PI);

	return yy;
}
// -----------------------------------------------------------------------
inline real FTU(real ang, real t[] = 0, real dt = 0)
{
	real ft = FT(ang, t, dt);
	real max1 = (1 + 0.5 + 0.25 + 0.125);
	real min1 = -max1;
	return (ft - min1) / (max1 - min1);
}
// -----------------------------------------------------------------------
inline real BlendFT(real h1, real h2, real alpha, real t[] = 0, real dt = 0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;

	alpha = FTU(alpha, t, dt);

	return h1 * (1 - alpha) + h2 * alpha;
}
// -----------------------------------------------------------------------
inline real FT2D1(real anglex, real angley, real tx[] = 0, real ty[] = 0, real dtx = 0, real dty = 0)
{
	if (tx == 0)
	{
		static real s_t0x[] = { rrnd(0, 1), rrnd(0, 1), rrnd(0, 1), rrnd(0, 1) };
		tx = s_t0x;
	}
	if (ty == 0)
	{
		static real s_t0y[] = { rrnd(0, 1), rrnd(0, 1), rrnd(0, 1), rrnd(0, 1) };
		ty = s_t0y;
	}

	real yy = 0;
	yy += 1 * sin(1 * anglex + (tx[0] + dtx) * PI) + 1 * sin(1 * angley + (ty[0] + dty) * PI);
	yy += 0.5 * sin(2 * anglex + (tx[1] + dtx) * PI) + 0.5 * sin(2 * angley + (ty[1] + dty) * PI);
	yy += 0.25 * sin(4 * anglex + (tx[2] + dtx) * PI) + 0.25 * sin(4 * angley + (ty[2] + dty) * PI);
	yy += 0.125 * sin(8 * anglex + (tx[3] + dtx) * PI) + 0.125 * sin(8 * angley + (ty[3] + dty) * PI);

	return yy;
}
// ------------------------------------------------------------------------------------------------
real FFT(real ax, real ay, real rndmap[1024][1024])
{
	real ft = 0;
	for (int i = 0; i < 8; i++)
	{
		real dz = blend(ax, ay, rndmap[100][i], 2);
		{
			real A = 2 * rndmap[0][i];
			real F = blend(20, 50, rndmap[1][i]);
			real T0 = 2 * PI * rndmap[2][i];
			ft += A * (0.5 + 0.5 * sin(T0 + F * dz * PI * 2));
		}
	}
	return ft / (8);
}

// -----------------------------------------------------------------------
inline real FT2D(real anglex, real angley, real tx[] = 0, real ty[] = 0, real dtx = 0, real dty = 0)
{
	static real s_fx[] = { rrnd(0, 8), rrnd(0, 4), rrnd(0, 2), rrnd(0, 1) };
	static real s_fy[] = { rrnd(0, 8), rrnd(0, 4), rrnd(0, 2), rrnd(0, 1) };
	if (tx == 0)
	{
		static real s_t0x[] = { rrnd(0, 1), rrnd(0, 1), rrnd(0, 1), rrnd(0, 1) };
		tx = s_t0x;
	}
	if (ty == 0)
	{
		static real s_t0y[] = { rrnd(0, 1), rrnd(0, 1), rrnd(0, 1), rrnd(0, 1) };
		ty = s_t0y;
	}

	real yy = 0;
	yy += 1 * sin(s_fx[0] * anglex + (tx[0] + dtx) * PI + s_fy[0] * angley + (ty[0] + dty) * PI);
	yy += .5 * sin(s_fx[1] * anglex + (tx[1] + dtx) * PI + s_fy[1] * angley + (ty[1] + dty) * PI);
	yy += .25 * sin(s_fx[2] * anglex + (tx[2] + dtx) * PI + s_fy[2] * angley + (ty[2] + dty) * PI);
	yy += .125 * sin(s_fx[3] * anglex + (tx[3] + dtx) * PI + s_fy[3] * angley + (ty[3] + dty) * PI);

	return yy;
}

// ------------------------------------------------
// 2d
real blend2d(real h1, real h2, real alphaX, real alphaY)
{
	int size;
	real alpha;
	alphaX < 0 ? alphaX = 0 : 0;
	alphaX > 1 ? alphaX = 1 : 0;
	alphaY < 0 ? alphaY = 0 : 0;
	alphaY > 1 ? alphaY = 1 : 0;
	size = IMAGESCALE;
	alpha = fheightmap[(int)(alphaX * size + 0.5) % size][(int)(alphaY * size + 0.5) % size];
	alpha = (alpha - fmin1) / (fmax1 - fmin1);

	return h1 * (1 - alpha) + h2 * alpha;
}
// **********************************************************************
// color blend
// **********************************************************************
inline int blendcor(int c1, int c2, real alpha, real power = 1.0)
{
	//alpha = alpha != 1 ? abs(alpha) - (int)(alpha) : 1;	
	alpha = alpha > 1 ? 1 : alpha;
	alpha = alpha < 0 ? 0 : alpha;

	if (power != 1.0)
		alpha = pow(alpha, power);
	return _RGB(GetRValue(c2) * alpha + GetRValue(c1) * (1 - alpha),
		GetGValue(c2) * alpha + GetGValue(c1) * (1 - alpha),
		GetBValue(c2) * alpha + GetBValue(c1) * (1 - alpha)
	);

}
inline int blendcorRGBA(int c1, int c2, real alpha, real power = 1.0)
{
	//alpha = alpha != 1 ? abs(alpha) - (int)(alpha) : 1;
	alpha = alpha > 1 ? 1 : alpha;
	alpha = alpha < 0 ? 0 : alpha;

	if (power != 1.0)
		alpha = pow(alpha, power);
	return RGBA(GetRValue(c2) * alpha + GetRValue(c1) * (1 - alpha),
		GetGValue(c2) * alpha + GetGValue(c1) * (1 - alpha),
		GetBValue(c2) * alpha + GetBValue(c1) * (1 - alpha),
		GetAValue(c2) * alpha + GetAValue(c1) * (1 - alpha)
	);
}

// ----------------------------------------------------------------------
// 2d
inline int blendcor2d(int c1, int c2, real alphaX, real alphaY, real power = 1.0)
{
	real alpha;
	int size = IMAGESCALE;
	alphaX < 0 ? alphaX = 0 : 0;
	alphaX > 1 ? alphaX = 1 : 0;
	alphaY < 0 ? alphaY = 0 : 0;
	alphaY > 1 ? alphaY = 1 : 0;

	alpha = fheightmap[(int)(alphaX * size) % size][(int)(alphaY * size) % size];
	alpha = (alpha - fmin1) / (fmax1 - fmin1);

	return blendcor(c1, c2, alpha, power);
}

// color add

inline int addcor(int c1, int c2)
{
	int r = GetRValue(c1) + GetRValue(c2);
	int g = GetGValue(c1) + GetGValue(c2);
	int b = GetBValue(c1) + GetBValue(c2);
	return RGB(r > 255 ? 255 : r,
		g > 255 ? 255 : g,
		b > 255 ? 255 : b);
}

// ----------------------------------------------------------------------
// lerp 
// ----------------------------------------------------------------------
inline real lerp(real h1, real h2, real alpha, real power = 1.0)
{
	if (power != 1.0)
		alpha = pow(alpha, power);

	return h1 * (1.0f - alpha) + h2 * alpha;
}
inline vector2 lerp(const vector2& v1, const vector2& v2, real alpha, real power = 1.0)
{
	if (power != 1.0)
		alpha = pow(alpha, power);

	return v1 * (1 - alpha) + v2 * alpha;
}
inline vector3 lerp(const vector3& v1, const vector3& v2, real alpha, real power = 1.0)
{
	if (power != 1.0)
		alpha = pow(alpha, power);

	return v1 * (1 - alpha) + v2 * alpha;
}
inline vectorn lerp(const vectorn& v1, const vectorn& v2, real alpha, real power = 1.0)
{
	if (power != 1.0)
		alpha = pow(alpha, power);

	return v1 * (1.0f - alpha) + v2 * alpha;
}
// ----------------------------------------------------------------------
// 模板实现 
// ----------------------------------------------------------------------
template<typename T>
inline T blend(const T& v1, const T& v2, real alpha, real power = 1.0)
{
	alpha = alpha < 0 ? 0 : alpha;
	alpha = alpha > 1 ? 1 : alpha;
	if (power != 1.0)
		alpha = pow(alpha, power);

	return v1 * (1.0f - alpha) + v2 * alpha;
}
template<typename T>
inline T lerp(const T& v1, const T& v2, real alpha, real power = 1.0)
{
	if (power != 1.0)
		alpha = pow(alpha, power);

	return v1 * (1.0f - alpha) + v2 * alpha;
}