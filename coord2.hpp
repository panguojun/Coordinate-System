/****************************************************************************************************
*						[Coordinate System (Coordinate Frame)]
*
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
*   The Coordinate System class is specifically encapsulated to simplify coordinate transformations
*   and derive geometric algorithms, capable of solving various problems related to coordinate system
*   transformations. The coordinate system operations exhibit a Lie group-like structure.
*   The coordinate system consists of three components: C = M (position) + S (scaling) * R (rotation).
*
*  *  *  *  *  *  *  *  *  *  *  *  Detailed Explanation  *  *  *  *  *  *  *  *  *  *  *  *  *  *
*   Coordinate system transformation is divided into three steps:
*					projection (/), translation (^), and restoration (*).
*
*   The coordinate system itself is denoted as C. Transformations between coordinate systems
*   can be expressed as G = C2 / C1 - I, where G represents the geometric gradient.
*					oper(/)  =  C1 * C2^-1
*					oper(\)  =  C1^-1 * C2
*
*   Specifically:
*   Define a vector V in an intrinsic coordinate system (assuming a flat space where vectors can
*   move freely without change). When observing V in a curved coordinate system, V appears different
*   at different points. Therefore, the coordinate system is position-dependent.
*
*   Take vectors V1 and V2 at adjacent points (1) and (2) respectively,
*   corresponding to coordinate systems C1 and C2. Then:
*					V  = V1 * C1 = V2 * C2 =>
*					V2 = V1 * C1 / C2, let R12 = C1 / C2 =>
*					V2 = V1 * R12
*
*   Based on the frame field combination operator theory proposed in this paper,
*   the direct geometric information extraction formula is:
*           G = (c₂·c₁⁻¹)/C₂ - I/C₁
*   where c is the intrinsic frame field and C is the embedding frame field.
*
*   The coordinate system can be used to compute spatial curvature. In the u,v coordinate system,
*   the curvature tensor is:
*					Ruv  =  Gu·Gv - Gv·Gu - G[u,v]
*   where:
*				 Gu = (c(u+du,v)·c⁻¹(u,v))/C(u+du,v) - I/C(u,v)
*				 Gv = (c(u,v+dv)·c⁻¹(u,v))/C(u,v+dv) - I/C(u,v)
*				 Connection vector: W = [U, V] (Lie bracket operation)
*				 G[u,v] = Gu·Wu + Gv·Wv
*
*   Compared with traditional methods, this framework avoids the complex Christoffel symbol
*   computation chain and directly extracts geometric invariants through frame field combinations,
*   offering higher computational efficiency and geometric intuitiveness.
*/

// **************************************************************************************************
//  |_
// UC     2d Rotation Coordinate System
// **************************************************************************************************
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
		return vec2(v.dot(c.ux), v.dot(c.uy));
	}
	// oper(/) = C1 * C2^-1
	ucoord2 operator / (const ucoord2& c) const
	{
		ucoord2 rc;
		rc.ux = vec2(ux.dot(c.ux), ux.dot(c.uy));
		rc.uy = vec2(uy.dot(c.ux), uy.dot(c.uy));
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

	// 角度
	real angle() const
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
	ucoord2 rotcopy(real angle) const
	{
		ucoord2 c = (*this);
		vec2 z = vector2::ang_len(angle, 1);
		c.ux = complex_mul(c.ux, z);
		c.uy = complex_mul(c.uy, z);
		return c;
	}
	void rot2dir(crvec2 _dir)
	{
		ux = _dir; uy = _dir.rotcopy(PI / 2);
	}
	void dump(const std::string& name = "") const
	{
		PRINT("----" << name << "---");
		PRINTVEC2(ux);
		PRINTVEC2(uy);
	}
};
const ucoord2 ucoord2::ZERO = { 0 };
const ucoord2 ucoord2::ONE = ucoord2();

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
	coord2(crvec2 p, crvec2 _s, real ang)
	{
		o = p;
		s = _s;

		vec2 z = vector2::ang_len(ang, 1);
		ux = complex_mul(ux, z);
		uy = complex_mul(uy, z);
	}
	operator vec2 () const
	{
		return o;
	}
	vec2 VX() const { return ux * s.x; }
	vec2 VY() const { return uy * s.y; }
	coord2 VC() const
	{
		return { VX(), VY() };
	}
	ucoord2 UC() const
	{
		return {ux, uy};
	}
	void UC(const ucoord2& uc)
	{
		ux = uc.ux;
		uy = uc.uy;
	}
	bool equal_dirs(const coord2& c) const
	{
		return ux == c.ux && uy == c.uy && o == c.o && s == c.s;
	}
	bool operator == (const coord2& c) const
	{
		return o == c.o && s == c.s && equal_dirs(c);
	}
	coord2 operator + (const coord2& c) const
	{
		coord2 rc;
		rc.ux = VX() + c.VX();
		rc.uy = VY() + c.VY();
		// rc.norm();
		rc.o = o + c.o;
		return rc;
	}
	void operator += (const coord2& c)
	{
		*this = *this + c;
	}
	coord2 operator + (const vec2& v) const
	{
		coord2 c = (*this);
		c.o += v;
		return c;
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
		// rc.norm();
		rc.o = o - c.o;
		return rc;
	}
	void operator -= (const coord2& c)
	{
		*this = *this - c;
	}
	coord2 operator - (const vec2& v) const
	{
		coord2 c = (*this);
		c.o -= v;
		return c;
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
#ifdef	NON_UNIFORM_SCALE
		rc.ux = (ux.x * s.x) * (c.ux * c.s.x) + (ux.y * s.x) * (c.uy * c.s.y);
		rc.uy = (uy.x * s.y) * (c.ux * c.s.x) + (uy.y * s.y) * (c.uy * c.s.y);
		rc.norm();
#else
		rc = ucoord2::operator*(c);
		rc.s = s * c.s;
#endif
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
		{// C*S 移动乘法
			c.o *= s;
		}
		return c;
	}
	void operator *= (real s)
	{
		*this = (*this) * s;
	}
	// 向量向坐标系投影
	friend vec2 operator / (crvec2 p, const coord2& c)
	{
		vec2 v = p - c.o;
		v /= c.s;
		return vec2(v.dot(c.ux), v.dot(c.uy));
	}
	// oper(/) = C1 * C2^-1
	coord2 operator / (const coord2& c) const
	{
		coord2 rc;
#ifdef	NON_UNIFORM_SCALE
		vec2 vx = VX();
		vec2 vy = VY();

		vec2 cvx = c.ux / c.s.x;
		vec2 cvy = c.uy / c.s.y;

		rc.ux = vec2(vx.dot(cvx), vx.dot(cvy));
		rc.uy = vec2(vy.dot(cvx), vy.dot(cvy));

		rc.norm();
#else
		rc = ucoord2::operator/(c);
		rc.s = s / c.s;
#endif
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
	vec2 pos() const
	{
		return o;
	}
	string serialise() const
	{
		return o.serialise() + "," + std::to_string(angle());
	}
	void dump(const string& name = "") const
	{
		PRINT("----" << name << "---");
		PRINTVEC2(o);
		PRINTVEC2(s);
		PRINTVEC2(ux);
		PRINTVEC2(uy);
	}
};
const coord2 coord2::ZERO = { ucoord2::ZERO, vec2::ZERO };
const coord2 coord2::ONE = coord2();



