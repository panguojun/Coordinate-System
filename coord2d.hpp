/*********************************************************************
*						【坐标系】
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
* 	坐标系类是我单独封装，用于简化坐标变换，衍生出许多算法，能解决一些
* 	坐标系变换相关的问题。
* 	坐标系的运算跟李群很相似。
*	坐标系由三个部分组成：C = M(位置） + S（缩放） * R（旋转）
*
*  *  *  *  *  *  *  *  *  *  详解  *  *  *  *  *  *  *  *  *  *  *  *
*	坐标系变换分为投影（*), 平移（^), 还原（*）三个步骤，以平移最精深：
*	坐标系本体符号 C，坐标系之间的变换可以写成G = C1//C2,GRAD梯度的意思
*			oper(/)  = C1 * C2^-1
*			oper(//) = C1^-1 * C2, oper(//) = grad()
*	具体来说：
*	定义一个内禀坐标系(假设它是平直空间，向量可以随意移动而不变)下V,在弯
*	曲坐标系下观察V，不同点上V是不同的，故而坐标系跟位置有关，取相邻两点
*	（1),(2)点处有向量V1,V2，对应坐标系C1,C2，那么：
*			V = V1 * C1 = V2 * C2 =>
*			V2 = V1 * C1 / C2, 令 G12 = C1 / C2 =>
*			V2 = V1 * G12
*
*	可以使用坐标系计算空间曲率，在u,v坐标系下黎曼曲率张量为：
*			Ruv = Gu*Gv - Gv*Gu - G[u, v]
*			其中：Gu = UG - ONE
*			     UG = C2 / C1
*			     [U, V] : 李括号运算
*/

//#define	Parallel_Projection		 // 非正交坐标系下平行投影
// *******************************************************************
//  |_
// UC     2d Rotation Coordinate System
// *******************************************************************
struct ucoord2
{
	static const ucoord2 ZERO;
	static const ucoord2 ONE;

	vec2 ux = vec2::UX;		// 方向
	vec2 uy = vec2::UY;

	ucoord2() {}
	ucoord2(crvec2 _ux, crvec2 _uy, crvec2 _uz)
	{
		ux = _ux; uy = _uy;
	}
	ucoord2(crvec2 _ux, crvec2 _uy)
	{
		ux = _ux;
		uy = _uy;
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
	coord2(const coord2& c)
	{
		ux = c.ux; uy = c.uy;
		s = c.s;
		o = c.o;
	}
	coord2(crvec2 _ux, crvec2 _uy, crvec2 _uz)
	{
		ux = _ux; uy = _uy;
	}
	coord2(crvec2 _ux, crvec2 _uy)
	{
		ux = _ux;
		uy = _uy;
	}
	coord2(crvec2 p)
	{
		o = p;
	}
	coord2(real ang)
	{
		ux.rot(ang);
		uy.rot(ang);
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
		return (*this) * coord2(ux * v.x, uy * v.y);
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
	// 梯度坐标系
	static coord2 grad(const coord2& c1, const coord2& c2)
	{
		return c1.reversed() * c2 - ONE;
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
const coord2 coord2::ZERO = { 0 };
const coord2 coord2::ONE = coord2();
#endif