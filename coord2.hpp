/************************************************************************************************
*				[Coordinate System]
*				   by Guojun Pan
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
*   The coordinate system class is separately encapsulated by me for
*   simplifying coordinate transformation and deriving many algorithms,
*   which can solve some problems related to coordinate system transformation.
*   The operation of the coordinate system is similar to Lie group.
*   The coordinate system consists of three parts: C = M (position) + S (scaling) * R (rotation).
*
*  *  *  *  *  *  *  *  *  *  *  *  Detailed Explanation  *  *  *  *  *  *  *  *  *  *  *  *  *  *
*   The coordinate system transformation is divided into three steps:
*   			projection (/), translation (^), and restoration (*).
*
*   The symbol of the coordinate system itself is C. The transformation between coordinate systems
*   can be written as G = C2 / C1 - I, where G means gradient.
*           		oper(/)  =  C1 * C2^-1
*           		oper(\)  =  C1^-1 * C2
*
*   Specifically:
*   Define an intrinsic coordinate system (assuming it is a flat space, and the vector can move freely
*   without changing) under V. Observing V in a curved coordinate system, V is different at different points.
*   Therefore, the coordinate system is related to the position.
*   Take vectors V1 and V2 at adjacent points (1) and (2) respectively,
*   corresponding to coordinate systems C1 and C2. Then:
*           		V  = V1 * C1 = V2 * C2 =>
*           		V2 = V1 * C1 / C2, let R12 = C1 / C2 =>
*           		V2 = V1 * R12
*
*   The coordinate system can be used to calculate spatial curvature. In the u,v coordinate system,
*   the Riemann curvature tensor is:
*           	Ruv  = 	Gu*Gv - Gv*Gu - G[u,v]
*           	where:  Gu = C2 / C1 - I
*                   	Connection vector: W = [U, V] (Lie bracket operation)
*                   	G[u,v] = Gu*Wu + Gv*Wv
*/
// *******************************************************************
//  |_
// UC     2d Rotation Coordinate System
// *******************************************************************
struct ucoord2
{
	static const ucoord2 ZERO;
	static const ucoord2 ONE;

	vec2 ux = vec2::UX;		// basis 单位化基向量
	vec2 uy = vec2::UY;

	ucoord2() {}
	ucoord2(crvec2 _ux, crvec2 _uy)
	{
		ux = _ux;
		uy = _uy;
	}
	ucoord2(crvec2 _ux)
	{
		ux = _ux;
		uy = ux.rotcopy(PI / 2);
	}
	ucoord2(real ang)
	{
		ux.rot(ang);
		uy.rot(ang);
	}
	void rot(real ang)
	{
		ux.rot(ang);
		uy.rot(ang);
	}
	bool is_same_dirs(const ucoord2& c) const
	{
		return ux == c.ux && uy == c.uy;
	}
	// 在坐标系下定义一个向量
	friend vec2 operator * (crvec2 p, const ucoord2& c)
	{
		return c.ux * (p.x) + c.uy * (p.y);
	}
	ucoord2 operator * (const ucoord2& c) const
	{
		ucoord2 rc;
		rc.ux = ux.x * c.ux + ux.y * c.uy;
		rc.uy = uy.x * c.ux + uy.y * c.uy;
		return rc;
	}
	// 向量向坐标系投影
	friend vec2 operator / (crvec2 v, const ucoord2& c)
	{
#ifdef Parallel_Projection
		{// 对于非正交情况
			return vec2(pl_dot(v, c.ux, c.uy), pl_dot(v, c.uy, c.ux));
		}
#endif
		return vec2(v.dot(c.ux), v.dot(c.uy));
	}
	// oper(/) = C1 * C2^-1
	ucoord2 operator / (const ucoord2& c) const
	{
		ucoord2 rc;
#ifdef Parallel_Projection
		{// 对于非正交情况
			rc.ux = vec2(pl_dot(ux, c.ux, c.uy) / c.s.x, pl_dot(ux, c.uy, c.ux) / c.s.y);
			rc.uy = vec2(pl_dot(uy, c.ux, c.uy) / c.s.x, pl_dot(uy, c.uy, c.ux) / c.s.y);
		}
#else
		rc.ux = vec2(ux.dot(c.ux), ux.dot(c.uy));
		rc.uy = vec2(uy.dot(c.ux), uy.dot(c.uy));
#endif
		return rc;
	}
	// oper(//) = C1^-1 * C2
	ucoord2 operator % (const ucoord2& c) const
	{
		return (*this).reversed() * c;
	}
	// 倒置
	void reverse()
	{
		(*this) = ONE / (*this);
	}
	ucoord2 reversed() const
	{
		return ONE / (*this);
	}
	// 梯度坐标系
	static ucoord2 grad(const ucoord2& c1, const ucoord2& c2)
	{
		return c1.reversed() * c2;
	}
	real dot(crvec2 v) const
	{
		return v.dot(ux) + v.dot(uy);
	}
	void dump(const std::string& name = "") const
	{
		PRINT("----" << name << "---");
		PRINTVEC2(ux);
		PRINTVEC2(uy);
	}
};
#ifdef PMDLL
const ucoord2 ucoord2::ZERO = { 0 };
const ucoord2 ucoord2::ONE = ucoord2();
#endif

// *******************************************************************
//  |_
// C     2d Coordinate System
// *******************************************************************
struct coord2 : ucoord2
{
	static const coord2 ZERO;
	static const coord2 ONE;

	vec2 s = vec2::ONE;		// 缩放
	vec2 o;				// 原点

	coord2() {}
	coord2(const ucoord2& c) : ucoord2(c)
	{
	}
	coord2(const ucoord2& c, crvec2 _s) : ucoord2(c)
	{
		s = _s;
	}
	coord2(const coord2& c) : ucoord2(c.ux, c.uy)
	{
		s = c.s;
		o = c.o;
	}
	coord2(crvec2 _ux, crvec2 _uy) : ucoord2(_ux, _uy)
	{
	}
	coord2(crvec2 p)
	{
		o = p;
	}
	coord2(real ang)
	{
		vec2 z = vector2::ang_len(ang, 1);
		ux = complex_mul(ux, z);
		uy = complex_mul(uy, z);
	}
	coord2(real ang, real _r)
	{
		vec2 z = vector2::ang_len(ang, 1);
		ux = complex_mul(ux, z);
		uy = complex_mul(uy, z);
		s *= _r;
	}
	coord2(crvec2 p, real ang)
	{
		o = p;
		vec2 z = vector2::ang_len(ang, 1);
		ux = complex_mul(ux, z);
		uy = complex_mul(uy, z);
	}
	vec2 VX() const { return ux * s.x; }
	vec2 VY() const { return uy * s.y; }

	bool is_same_dirs(const coord2& c) const
	{
		return ux == c.ux && uy == c.uy && o == c.o && s == c.s;
	}
	coord2 operator + (const coord2& c) const
	{
		coord2 rc;
		rc.ux = VX() + c.VX();
		rc.uy = VY() + c.VY();
		rc.norm();
		rc.o = o + c.o;
		return rc;
	}
	void operator += (const coord2& c)
	{
		*this = *this + c;
	}
	coord2 operator + (const vec2& v) const
	{
		coord2 c; c.o = v; c.s = vec2::ZERO;
		return (*this) + c;
	}
	void operator += (const vec2& v)
	{
		*this = *this + v;
	}
	coord2 operator - (const coord2& c) const
	{
		coord2 rc;
		rc.ux = VX() - c.VX();
		rc.uy = VY() - c.VY();
		rc.norm();
		rc.o = o - c.o;
		return rc;
	}
	void operator -= (const coord2& c)
	{
		*this = *this - c;
	}
	coord2 operator - (const vec2& v) const
	{
		coord2 c; c.o = v; c.s = vec2::ZERO;
		return (*this) - c;
	}
	void operator -= (const vec2& v)
	{
		*this = *this - v;
	}
	// 在坐标系下定义一个向量
	friend vec2 operator * (crvec2 p, const coord2& c)
	{
		return c.ux * (c.s.x * p.x) + c.uy * (c.s.y * p.y) + c.o;
	}
	coord2 operator * (crvec2 v) const
	{
		return  (*this) * coord2(vec2::UX * v.x, vec2::UY * v.y);
	}
	void operator *= (crvec2 v)
	{
		(*this) *= coord2(vec2::UX * v.x, vec2::UY * v.y);
	}
	coord2 operator * (const coord2& c) const
	{
		coord2 rc;
		rc.ux = ux.x * c.ux + ux.y * c.uy;
		rc.uy = uy.x * c.ux + uy.y * c.uy;
		rc.s = s * c.s;
		rc.o = c.o + o.x * c.s.x * c.ux + o.y * c.s.y * c.uy;
		return rc;
	}
	void operator *= (const coord2& c)
	{
		(*this) = (*this) * c;
	}
	coord2 operator * (real s) const
	{
		coord2 c = *this;
		//{// C*S 缩放乘法
		//	c.s.x *= s; c.s.y *= s;
		//}
		{// C*S 移动乘法
			c.o.x *= s; c.o.y *= s;
		}
		return c;
	}
	void operator *= (real s)
	{
		*this = (*this) * s;
	}
#ifdef Parallel_Projection
	// 非正交坐标系下平行投影 Parallel projection
	static real pl_dot(crvec2 v, crvec2 ax1, crvec2 ax2)
	{
		real co = ax1.dot(ax2);
		real si = sqrt(1 - co * co);
		real sc = (co / si);
		return (v.dot(ax1) - v.cross(ax1) * sc);
	}
#endif
	// 向量向坐标系投影
	friend vec2 operator / (crvec2 p, const coord2& c)
	{
		vec2 v = p - c.o;
#ifdef Parallel_Projection
		{// 对于非正交情况
			return vec2(pl_dot(v, c.ux, c.uy) / c.s.x, pl_dot(v, c.uy, c.ux) / c.s.y);
		}
#endif
		return vec2(v.dot(c.ux) / c.s.x, v.dot(c.uy) / c.s.y);
	}
	// oper(/) = C1 * C2^-1
	coord2 operator / (const coord2& c) const
	{
		coord2 rc;
#ifdef Parallel_Projection
		{// 对于非正交情况
			rc.ux = vec2(pl_dot(ux, c.ux, c.uy) / c.s.x, pl_dot(ux, c.uy, c.ux) / c.s.y);
			rc.uy = vec2(pl_dot(uy, c.ux, c.uy) / c.s.x, pl_dot(uy, c.uy, c.ux) / c.s.y);
		}
#else
		rc.ux = vec2(ux.dot(c.ux) / c.s.x, ux.dot(c.uy) / c.s.y);
		rc.uy = vec2(uy.dot(c.ux) / c.s.x, uy.dot(c.uy) / c.s.y);
#endif
		rc.s = s / c.s;
		rc.o = o - c.o;
		rc.o = vec2(rc.o.dot(c.ux) / c.s.x, rc.o.dot(c.uy) / c.s.y);
		return rc;
	}
	coord2 operator / (crvec2 v) const
	{
		return (*this) / coord2(ux * v.x, uy * v.y);
	}
	// oper(//) = C1^-1 * C2
	coord2 operator % (const coord2& c) const
	{
		return (*this).reversed() * c;
	}
	coord2 operator ^ (real f) const
	{
		real ang = ux.angle() * f;
		real rad = ::exp(::log(ux.len()) * f);
		return coord2(ang, rad);
	}
	void norm(bool bscl = true)
	{
#define ISZERO(a) (fabs(a) < 1e-10)
		s.x = ux.len(); if (!ISZERO(s.x)) ux /= s.x;
		s.y = uy.len(); if (!ISZERO(s.y)) uy /= s.y;

		if (!bscl)
			s = vec2::ONE;
	}
	// 倒置
	void reverse()
	{
		(*this) = ONE / (*this);
	}
	coord2 reversed() const
	{
		return ONE / (*this);
	}
	real dot(crvec2 v) const
	{
		return v.dot(ux) * s.x + v.dot(uy) * s.y;
	}
	// 梯度
	static coord2 grad(const coord2& c1, const coord2& c2)
	{
		return c1.reversed() * c2 - ONE;
	}
	// 位置
	inline vec2 pos() const
	{
		return o;
	}
	// 角度
	inline real angle() const
	{
		return ux.angle();
	}
	// 旋转
	void rot(real angle)
	{
		vec2 z = vector2::ang_len(angle, 1);
		ux = complex_mul(ux, z);
		uy = complex_mul(uy, z);
	}
	coord2 rotcopy(real angle) const
	{
		coord2 c = (*this);
		vec2 z = vector2::ang_len(angle, 1);
		c.ux = complex_mul(c.ux, z);
		c.uy = complex_mul(c.uy, z);
		return c;
	}
	void rot2dir(crvec2 _dir)
	{
		ux = _dir; uy = _dir.rotcopy(PI / 2);
	}
	std::string serialise() const
	{
		return o.serialise() + "," + std::to_string(angle());
	}
	void dump(const std::string& name = "") const
	{
		PRINT("----" << name << "---");
		PRINTVEC2(ux);
		PRINTVEC2(uy);
		PRINTVEC2(s);
		PRINTVEC2(o);
	}
};
#ifdef PMDLL
const coord2 coord2::ZERO = { ucoord2::ZERO, vec2::ZERO };
const coord2 coord2::ONE = coord2();
#endif
