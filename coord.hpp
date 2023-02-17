/*********************************************************************
*				坐标系
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
* 	坐标系类是我单独封装，用于简化坐标变换，衍生出许多算法，能解决一些
* 	坐标系变换相关的问题。
* 	坐标系的运算跟李群很相似。
*	坐标系由三个部分组成：C = M(位置） + S（缩放） * R（旋转）
*  *  *  *  *  *  *  *  *  *  详解  *  *  *  *  *  *  *  *  *  *  *  *
*	坐标系变换分为投影（*), 平移（^), 还原（*）三个步骤，以平移最精深：
*	坐标系本体符号 C，坐标系之间的变换可以写成G = C1//C2,GRAD梯度的意思
*			oper(/)  = C1 * C2^-1
*			oper(//) = C1^-1 * C2, oper(//) = gradcoord()
*	坐标系的李括号: [C1,C2] = C1*C2 - C2*C1
*	具体来说：
*	定义一个内禀坐标系(假设它是平直空间，向量可以随意移动而不变)下V,在弯
*	曲坐标系下观察V，不同点上V是不同的，故而坐标系跟位置有关，取相邻两点
*	（1),(2)点处有向量V1,V2，对应坐标系C1,C2，那么：
*			V = V1 * C1 = V2 * C2 =>
*			V2 = V1 * C1 / C2, 令 G12 = C1 / C2 =>
*			V2 = V1 * G12
*
*	在弯曲坐标系下内禀坐标系x,y轴的平行线投影得到的u,v曲线上G12分别在两个
*	方向上对应Gu,Gv, 从(u1,v1)到(u2,v2) 计算两个路径的差别再加上修正项可
*	得曲率公式为：
*			Ruv = Gu*Gv - Gv*Gu * Gu^Wu * Gv^Wv
*			W = (U + V*Gu) - (V + U*Gv)
*/

//#define	Parallel_Projection		 // 非正交坐标系下平行投影
// ******************************************************************
//  |/_
// UC     3d Rotation Coordinate System
// ******************************************************************
struct ucoord3
{
	static const ucoord3 ZERO;
	static const ucoord3 ONE;

	vec3 ux = vec3::UX;		// 方向
	vec3 uy = vec3::UY;
	vec3 uz = vec3::UZ;

	ucoord3() {}
	ucoord3(const ucoord3& c)
	{
		ux = c.ux; uy = c.uy; uz = c.uz;
	}
	ucoord3(const vec3& _ux, const vec3& _uy, const vec3& _uz)
	{
		ux = _ux; uy = _uy; uz = _uz;
	}
	ucoord3(const vec3& _ux, const vec3& _uy)
	{
		ux = _ux; uy = _uy; uz = ux.cross(uy);
	}
	ucoord3(real ang, const vec3& ax)
	{
		ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);
	}
	ucoord3(real pit, real yaw, real rol)
	{
		ux.rot(pit, vec3::UX);
		uy.rot(yaw, vec3::UY);
		uz.rot(rol, vec3::UZ);
	}
	ucoord3(const quaternion& q)
	{
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}
	// 旋转差
	void fromvecsR(const vec3& v1, const vec3& v2)
	{
		quaternion q;
		q.fromvectors(v1, v2);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}
	void fromaxvecs(const vec3& ax, const vec3& v1, const vec3& v2)
	{
		vec3 pv1 = v1.crossdot(ax);
		vec3 pv2 = v2.crossdot(ax);
		real ang = acos(pv1.dot(pv2));
		quaternion q;q.fromangleaxis(ang, ax);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;;
	}
	quaternion toquat() const
	{
		vec3 pyr = coord2eulers();
		quaternion q;
		q.fromeuler(pyr.x, pyr.y, pyr.z);
		return q;
	}
	inline bool is_same_dirs(const ucoord3& c) const
	{
		return ux == c.ux && uy == c.uy && uz == c.uz;
	}
	bool operator == (const ucoord3& c) const
	{
		return is_same_dirs(c);
	}
	bool operator != (const ucoord3& c) const
	{
		return !is_same_dirs(c);
	}

	// 乘法：在坐标系下定义一个向量，或者向量向父空间还原
	friend vec3 operator * (const vec3& p, const ucoord3& c)
	{
		return c.ux * (p.x) + c.uy * (p.y) + c.uz * (p.z);
	}
	ucoord3 operator * (const ucoord3& c) const
	{// Cchild * Cparent * ...
		ucoord3 rc;
		rc.ux = ux.x * c.ux + ux.y * c.uy + ux.z * c.uz;
		rc.uy = uy.x * c.ux + uy.y * c.uy + uy.z * c.uz;
		rc.uz = uz.x * c.ux + uz.y * c.uy + uz.z * c.uz;

		return rc;
	}
	void operator *= (const ucoord3& c)
	{
		*this = (*this) * c;
	}
	ucoord3 operator * (const quaternion& q) const
	{
		ucoord3 rc = *this;
		rc.ux = q * ux;
		rc.uy = q * uy;
		rc.uz = q * uz;
		return rc;
	}
	void operator *= (const quaternion& q)
	{
		*this = (*this) * q;
	}

	// 除法：向量向坐标系投影
#ifdef Parallel_Projection
	// 非正交坐标系下平行投影 Parallel projection
	static real pl_prj(const vec3& v, const vec3& ax1, const vec3& ax2)
	{
		vec3 ax = ax1.cross(ax2); ax.norm();
		real co = ax1.dot(ax2);
		real si = sqrt(1 - co * co);
		real sc = (co / si);
		return (v.dot(ax1) - v.cross(ax1).dot(ax) * sc);
	}

#define PL_PRJ3(v) vec3( \
				pl_prj(v-c.uz*v.dot(c.uz), c.ux, c.uy) / c.s.x, \
				pl_prj(v-c.ux*v.dot(c.ux), c.uy, c.uz) / c.s.y, \
				pl_prj(v-c.uy*v.dot(c.uy), c.uz, c.ux) / c.s.z)
#endif
	friend vec3 operator / (const vec3& v, const ucoord3& c)
	{
#ifdef Parallel_Projection
		{// 对于非正交情况
			return vec3(
				pl_prj(v - c.uz * v.dot(c.uz), c.ux, c.uy) / c.s.x,
				pl_prj(v - c.ux * v.dot(c.ux), c.uy, c.uz) / c.s.y,
				pl_prj(v - c.uy * v.dot(c.uy), c.uz, c.ux) / c.s.z);
		}
#endif
		return vec3(v.dot(c.ux), v.dot(c.uy), v.dot(c.uz));
	}
	friend void operator /= (vec3& v, const ucoord3& c)
	{
		v = v / c;
	}
	// oper(/) = C1 * C2^-1
	ucoord3 operator / (const ucoord3& c) const
	{
		ucoord3 rc;
#ifdef Parallel_Projection
		{// 对于非正交情况
			rc.ux = PL_PRJ3(ux);
			rc.uy = PL_PRJ3(uy);
			rc.uz = PL_PRJ3(uz);
		}
#else
		rc.ux = vec3(ux.dot(c.ux), ux.dot(c.uy), ux.dot(c.uz));
		rc.uy = vec3(uy.dot(c.ux), uy.dot(c.uy), uy.dot(c.uz));
		rc.uz = vec3(uz.dot(c.ux), uz.dot(c.uy), uz.dot(c.uz));
#endif
		return rc;
	}
	void operator /= (const ucoord3& c)
	{
		*this = (*this) / c;
	}
	// oper(//) = C1^-1 * C2
	ucoord3 operator % (const ucoord3& c) const
	{
		return (*this).reversed() * c;
	}
	// oper(^)
	// 相空间的乘法运算,Ce^(th*v)
	// 如C表示某向量A在两点间的旋转，
	// 融合向量0<v<1,c=C^v; v=0时c=ONE,v=1时c=C
	ucoord3 operator ^ (const vec3& v) const
	{
		ucoord3 c = *this;
		c.ux = lerp(vec3::UX, c.ux, v.x); c.ux.norm();
		c.uy = lerp(vec3::UY, c.uy, v.y); c.uy.norm();
		c.uz = lerp(vec3::UZ, c.uz, v.z); c.uz.norm();
		
		return c;
	}
	// 转置(坐标轴交换）
	void transpose()
	{
		vec3 _ux = vec3(ux.x, uy.x, uz.x);
		vec3 _uy = vec3(ux.y, uy.y, uz.y);
		vec3 _uz = vec3(ux.z, uy.z, uz.z);
		ux = _ux;uy = _uy;uz = _uz;
	}
	ucoord3 transposed()
	{
		ucoord3 c = (*this);
		c.ux = vec3(ux.x, uy.x, uz.x);
		c.uy = vec3(ux.y, uy.y, uz.y);
		c.uz = vec3(ux.z, uy.z, uz.z);
		return c;
	}
	// 倒置
	void reverse()
	{
		(*this) = ONE / (*this);
	}
	ucoord3 reversed() const
	{
		return ONE / (*this);
	}
	// 翻转
	void flipX()
	{
		ux = -ux;
	}
	void flipY()
	{
		uy = -uy;
	}
	void flipZ()
	{
		uz = -uz;
	}
	void rot(real ang, const vec3& ax)
	{
		ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);
	}
	vec3 dir() const
	{
		return (ux + uy + uz).normlized();
	}
	// 本征向量（坐标系作为旋转变换时候的特征）
	vec3 eigenvec() const
	{
		return toquat().axis();
	}
	real dot(const vec3& v) const
	{
		return v.dot(ux) + v.dot(uy) + v.dot(uz);
	}
	real dot(const ucoord3& c) const
	{
		return c.ux.dot(ux) + c.uy.dot(uy) + c.uz.dot(uz);
	}
	// 由电磁场计算引出的叉乘
	ucoord3 cross(const ucoord3& c) const
	{
		return ucoord3(
			vec3::UX * (uy.dot(c.uz) - uz.dot(c.uy)),
			vec3::UY * (uz.dot(c.ux) - ux.dot(c.uz)),
			vec3::UZ * (ux.dot(c.uy) - uy.dot(c.ux))
		);
	}
	// v1 x v2 = v1 * (C x v2)
	ucoord3 cross(const vec3& v) const
	{
		return ucoord3(
			ux.cross(v),
			uy.cross(v),
			uz.cross(v)
		);
	}
	// 坐标系到欧拉角
	vec3 coord2eulers() const
	{
		const ucoord3& rm = *this;
		float sy = sqrt(rm.ux.x * rm.ux.x + rm.uy.x * rm.uy.x);
		bool singular = sy < 1e-6;

		float x, y, z;
		if (!singular)
		{
			x = atan2(rm.uz.y, rm.uz.z);
			y = atan2(-rm.uz.x, sy);
			z = atan2(rm.uy.x, rm.ux.x);
		}
		else
		{
			x = atan2(-rm.uy.z, rm.uy.y);
			y = atan2(-rm.uz.x, sy);
			z = 0;
		}
		//PRINT("rx: " << x * 180 / PI << ", ry: " << y * 180 / PI << ", rz: " << z * 180 / PI);
		//PRINT("rx: " << x << ", ry: " << y  << ", rz: " << z);
		return vec3(x, y, z);
	}

	// 梯度坐标系 = 梯度 X 切空间
	// 相当于一阶坐标系的导数
	static ucoord3 gradcoord(const ucoord3& c1, const ucoord3& c2)
	{
		return c1.reversed() * c2;
	}
	void dump(const std::string& name = "") const
	{
		PRINT("----" << name << "---");
		PRINTV3(ux);
		PRINTV3(uy);
		PRINTV3(uz);
	}
};
const ucoord3 ucoord3::ZERO = { 0 };
const ucoord3 ucoord3::ONE = ucoord3();

// ******************************************************************
//  |/_
// C     3d Coordinate System
// ******************************************************************
struct coord3 : ucoord3
{
	static const coord3 ZERO;
	static const coord3 ONE;

	vec3 s = vec3::ONE;		// 缩放
	vec3 o;					// 原点

	coord3() {}
	coord3(const coord3& c)
	{
		ux = c.ux; uy = c.uy; uz = c.uz;
		s = c.s;
		o = c.o;
	}
	coord3(const vec3& _ux, const vec3& _uy, const vec3& _uz)
	{
		ux = _ux; uy = _uy; uz = _uz;
	}
	coord3(const vec3& _ux, const vec3& _uy)
	{
		ux = _ux; uy = _uy; uz = ux.cross(uy);
	}
	explicit coord3(const vec3& _p)
	{
		o = _p;
	}
	coord3(real ang, const vec3& ax)
	{
		ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);
	}
	coord3(real pit, real yaw, real rol)
	{
		ux.rot(pit, vec3::UX);
		uy.rot(yaw, vec3::UY);
		uz.rot(rol, vec3::UZ);
	}
	coord3(const quaternion& q)
	{
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}
	// 移动差
	void fromvecsT(const vec3& v1, const vec3& v2)
	{
		o = v2 - v1;
	}

	vec3 VX() const { return ux * s.x; }
	vec3 VY() const { return uy * s.y; }
	vec3 VZ() const { return uz * s.z; }

	vec3 X() const { return ux * s.x + vec3::UX * o.x; }
	vec3 Y() const { return uy * s.y + vec3::UY * o.y; }
	vec3 Z() const { return uz * s.z + vec3::UZ * o.z; }

	// 归一化的正交坐标系
	coord3 ucoord() const
	{
		coord3 c = *this;
		c.norm(false);
		c.o = vec3::ZERO;
		return c;
	}
	// 向量坐标系 = 方向 X 缩放
	inline coord3 vcoord()
	{
		coord3 c = *this;
		c.o = vec3::ZERO;
		return c;
	}
	// 位置
	inline vec3 pos()
	{
		return o;
	}
	// 总向量
	inline vec3 sumvec()
	{
		return o + VX() + VY() + VZ();
	}
	// 向量
	vec3 tovec() const
	{
		return ux * s.x + uy * s.y + uz * s.z;
	}
	// 旋转
	quaternion toquat() const
	{
		coord3 c = ucoord();
		vec3 pyr = c.coord2eulers();
		quaternion q;
		q.fromeuler(pyr.x, pyr.y, pyr.z);
		return q;
	}
	inline bool is_same_dirs(const coord3& c) const
	{
		return ux == c.ux && uy == c.uy && uz == c.uz;
	}
	bool operator == (const coord3& c) const
	{
		return o == c.o && s == c.s && is_same_dirs(c);
	}
	bool operator != (const coord3& c) const
	{
		return o != c.o || s != c.s || !is_same_dirs(c);
	}

	// +/- 运算
	coord3 operator + (const coord3& c) const
	{
		coord3 rc = *this;
		rc.ux = VX() + c.VX();
		rc.uy = VY() + c.VY();
		rc.uz = VZ() + c.VZ();
		rc.norm();
		rc.o = o + c.o;
		return rc;
	}
	void operator += (const coord3& c)
	{
		*this = (*this) + c;
	}
	coord3 operator + (const vec3& v) const
	{// C+V 移动
		coord3 rc = (*this);
		rc.o = o + v;
		return rc;
	}
	void operator += (const vec3& v)
	{
		*this = *this + v;
	}
	friend vec3 operator + (const vec3& p, const coord3& c)
	{
		return p + c.o;
	}
	friend void operator += (vec3& p, const coord3& c)
	{
		p = p + c;
	}
	friend vec3 operator - (const vec3& p, const coord3& c)
	{
		return p - c.o;
	}
	friend void operator -= (vec3& p, const coord3& c)
	{
		p = p - c;
	}
	coord3 operator - (const coord3& c) const
	{
		coord3 rc = (*this);
		rc.ux = VX() - c.VX();
		rc.uy = VY() - c.VY();
		rc.uz = VZ() - c.VZ();
		rc.norm();
		rc.o = o - c.o;
		return rc;
	}
	coord3 operator - (const vec3& v) const
	{
		coord3 rc = (*this);
		rc.o = o - v;
		return rc;
	}
	void operator -= (const vec3& v)
	{
		*this = *this - v;
	}

	// 乘法：在坐标系下定义一个向量
	friend vec3 operator * (const vec3& p, const coord3& c)
	{
		return c.ux * (c.s.x * p.x) + c.uy * (c.s.y * p.y) + c.uz * (c.s.z * p.z) + c.o;
	}
	friend void operator *= (vec3& p, const coord3& c)
	{
		p = p * c;
	}
	coord3 operator * (const vec3& v) const
	{// C*V 缩放乘法
		coord3 c = *this;
		c.s.x *= v.x; c.s.y *= v.y; c.s.z *= v.z;
		c.o.x *= v.x; c.o.y *= v.y; c.o.z *= v.z;
		return c;
	}
	void operator *= (const vec3& v)
	{
		*this = (*this) * v;
	}
	coord3 operator * (real s) const
	{// C*S 缩放乘法
		coord3 c = *this;
		c.s.x *= s; c.s.y *= s; c.s.z *= s;
		return c;
	}
	void operator *= (real s)
	{
		*this = (*this) * s;
	}
	coord3 operator * (const coord3& c) const
	{// Cchild * Cparent * ...
		coord3 rc;
		rc.ux = ux.x * c.ux + ux.y * c.uy + ux.z * c.uz;
		rc.uy = uy.x * c.ux + uy.y * c.uy + uy.z * c.uz;
		rc.uz = uz.x * c.ux + uz.y * c.uy + uz.z * c.uz;

		rc.s = s * c.s;
		rc.o = c.o + o.x * c.s.x * c.ux  + o.y * c.s.y * c.uy + o.z * c.s.z * c.uz;
		return rc;
	}
	void operator *= (const coord3& c)
	{
		*this = (*this) * c;
	}

	// 除法：向量向坐标系投影 注意：要保证ux,uy,uz是单位向量！
#ifdef Parallel_Projection
	// 非正交坐标系下平行投影 Parallel projection
	static real pl_prj(const vec3& v, const vec3& ax1, const vec3& ax2)
	{
		vec3 ax = ax1.cross(ax2); ax.norm();
		real co = ax1.dot(ax2);
		real si = sqrt(1 - co * co);
		real sc = (co / si);
		return (v.dot(ax1) - v.cross(ax1).dot(ax) * sc);
	}

#define PL_PRJ3(v) vec3( \
				pl_prj(v-c.uz*v.dot(c.uz), c.ux, c.uy) / c.s.x, \
				pl_prj(v-c.ux*v.dot(c.ux), c.uy, c.uz) / c.s.y, \
				pl_prj(v-c.uy*v.dot(c.uy), c.uz, c.ux) / c.s.z)
#endif
	friend vec3 operator / (const vec3& p, const coord3& c)
	{
		vec3 v = p - c.o;
#ifdef Parallel_Projection
		{// 对于非正交情况
			return vec3(
				pl_prj(v - c.uz * v.dot(c.uz), c.ux, c.uy) / c.s.x,
				pl_prj(v - c.ux * v.dot(c.ux), c.uy, c.uz) / c.s.y,
				pl_prj(v - c.uy * v.dot(c.uy), c.uz, c.ux) / c.s.z);
		}
#endif
		return vec3(v.dot(c.ux) / c.s.x, v.dot(c.uy) / c.s.y, v.dot(c.uz) / c.s.z);
	}
	friend void operator /= (vec3& p, const coord3& c)
	{
		p = p / c;
	}
	coord3 operator / (const vec3& v) const
	{// C/V 缩放除法
		coord3 c = *this;
		c.s.x /= v.x; c.s.y /= v.y; c.s.z /= v.z;
		c.o.x /= v.x; c.o.y /= v.y; c.o.z /= v.z;
		return c;
	}
	void operator /= (const vec3& v)
	{
		*this = (*this) / v;
	}

	coord3 operator / (real s) const
	{// C/S 缩放除法
		coord3 c = *this;
		c.s.x /= s; c.s.y /= s; c.s.z /= s;
		return c;
	}
	void operator /= (real s)
	{
		*this = (*this) / s;
	}
	// oper(/) = C1 * C2^-1
	coord3 operator / (const coord3& c) const
	{
		coord3 rc;
#ifdef Parallel_Projection
		{// 对于非正交情况
			rc.ux = PL_PRJ3(ux);
			rc.uy = PL_PRJ3(uy);
			rc.uz = PL_PRJ3(uz);
		}
#else
		rc.ux = vec3(ux.dot(c.ux), ux.dot(c.uy), ux.dot(c.uz));
		rc.uy = vec3(uy.dot(c.ux), uy.dot(c.uy), uy.dot(c.uz));
		rc.uz = vec3(uz.dot(c.ux), uz.dot(c.uy), uz.dot(c.uz));
#endif
		rc.s = s / c.s;
		rc.o = o - c.o;
		rc.o = vec3(rc.o.dot(c.ux) / c.s.x, rc.o.dot(c.uy) / c.s.y, rc.o.dot(c.uz) / c.s.z);
		return rc;
	}
	void operator /= (const coord3& c)
	{
		*this = (*this) / c;
	}
	// oper(//) = C1^-1 * C2
	coord3 operator % (const coord3& c) const
	{
		return (*this).reversed() * c;
	}

	// 归一化
	void norm(bool bscl = true)
	{
#define ISZERO(a) (fabs(a) < 1e-10)
		s.x = ux.len(); if (!ISZERO(s.x)) ux /= s.x;
		s.y = uy.len(); if (!ISZERO(s.y)) uy /= s.y;
		s.z = uz.len(); if (!ISZERO(s.z)) uz /= s.z;
		if (!bscl)
			s = vec3::ONE;
	}
	// 倒置
	void reverse()
	{
		(*this) = ONE / (*this);
	}
	coord3 reversed() const
	{
		return ONE / (*this);
	}
	real dot(const vec3& v) const
	{
		return v.dot(ux) * s.x + v.dot(uy) * s.y + v.dot(uz) * s.z;
	}
	real dot(const coord3& c) const
	{
		return c.VX().dot(VX()) + c.VY().dot(VY()) + c.VZ().dot(VZ());
	}
	// 由李符号引出的叉乘，更加符合群论
	coord3 lie_cross(const coord3& c) const
	{
		return (*this) * c - c * (*this);
	}
	// 由电磁场计算引出的叉乘
	coord3 cross(const coord3& c) const
	{
		vec3 vx = VX();
		vec3 vy = VY();
		vec3 vz = VZ();

		vec3 cvx = c.VX();
		vec3 cvy = c.VY();
		vec3 cvz = c.VZ();

		return coord3(
			vec3::UX * (vy.dot(cvz) - vz.dot(cvy)),
			vec3::UY * (vz.dot(cvx) - vx.dot(cvz)),
			vec3::UZ * (vx.dot(cvy) - vy.dot(cvx))
		);
	}
	// v1 x v2 = v1 * (C x v2)
	coord3 cross(const vec3& v) const
	{
		vec3 vx = VX();
		vec3 vy = VY();
		vec3 vz = VZ();

		return coord3(
			vx.cross(v),
			vy.cross(v),
			vz.cross(v)
		);
	}
	void dump(const std::string& name = "") const
	{
		PRINT("----" << name << "---");
		PRINTV3(ux);
		PRINTV3(uy);
		PRINTV3(uz);
		PRINTV3(s);
		PRINTV3(o);
	}
};
const coord3 coord3::ZERO = { 0 };
const coord3 coord3::ONE = coord3();
