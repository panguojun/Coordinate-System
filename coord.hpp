/************************************************************************************************
*				[Coordinate System]
* 
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
*   The coordinate system class is separately encapsulated by me for
*   simplifying coordinate transformation and deriving many algorithms,
*   which can solve some problems related to coordinate system transformation.
*   The operation of the coordinate system is similar to Lie group.
*   The coordinate system consists of three parts: C = M (position) + S (scaling) * R (rotation).
*
*  *  *  *  *  *  *  *  *  *  *  *  Detailed Explanation  *  *  *  *  *  *  *  *  *  *  *  *  *  *	
*   The coordinate system transformation is divided into three steps: 
*   		projection (/), translation (^), and restoration (*).
*
*   The symbol of the coordinate system itself is C. The transformation between coordinate systems
*   can be written as G = C2 / C1 - I, where G means gradient.
*           	oper(/)  =  C1 * C2^-1
*           	oper(\)  =  C1^-1 * C2
*
*   Specifically:
*   Define an intrinsic coordinate system (assuming it is a flat space, and the vector can move freely
*   without changing) under V. Observing V in a curved coordinate system, V is different at different points.
*   Therefore, the coordinate system is related to the position.
*   Take vectors V1 and V2 at adjacent points (1) and (2) respectively,
*   corresponding to coordinate systems C1 and C2. Then:
*           	V  = V1 * C1 = V2 * C2 =>
*           	V2 = V1 * C1 / C2, let G12 = C1 / C2 =>
*           	V2 = V1 * G12
*
*   The coordinate system can be used to calculate spatial curvature. In the u,v coordinate system,
*   the Riemann curvature tensor is:
*           Ruv = Gu*Gv - Gv*Gu - G[u,v]
*           where:  Gu = C2 / C1 - I
*                   Connection vector: W = [U, V] (Lie bracket operation)
*                   G[u,v] = GuˆWu * GvˆWv 
*/

//#define	Parallel_Projection		 // 非正交坐标系下平行投影
// ********************************************************************************************
//  |/_
// UC     3d Rotation Coordinate System
// ********************************************************************************************
struct ucoord3
{
	static const ucoord3 ZERO;
	static const ucoord3 ONE;

	vec3 ux = vec3::UX;		// base 单位化向量
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
	// 轴，向量1，2
	void fromaxvecs(const vec3& ax, const vec3& v1, const vec3& v2)
	{
		vec3 pv1 = v1.crossdot(ax);
		vec3 pv2 = v2.crossdot(ax);
		real ang = acos(pv1.dot(pv2));
		quaternion q; q.ang_axis(ang, ax);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;;
	}
	// uy方向 推测ux,uz
	void fromuy(const vec3& _uy)
	{
		quat q; q.fromvectors(uy, _uy);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}
	void fromquat(const quaternion& q)
	{
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}
	void frompyr(real pit, real yaw, real rol)
	{
		quaternion q(pit, yaw, rol);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}
	inline void frompyr(crvec pyr)
	{
		fromquat(quaternion(pyr.x, pyr.y, pyr.z));
	}
	inline vec3 topyr()
	{
		return coord2eulers();
	}
	inline bool same_dirs(const ucoord3& c) const
	{
		return ux == c.ux && uy == c.uy && uz == c.uz;
	}
	inline bool operator == (const ucoord3& c) const
	{
		return same_dirs(c);
	}
	inline bool operator != (const ucoord3& c) const
	{
		return !same_dirs(c);
	}

	// 乘法：在坐标系下定义一个向量，或者向量向父空间还原
	inline friend vec3 operator * (const vec3& p, const ucoord3& c)
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
	inline friend quaternion operator * (const quaternion& q, const ucoord3& c)
	{
		return q * c.toquat();
	}
	ucoord3 operator * (const quaternion& q) const
	{
		ucoord3 rc;
		rc.ux = q * ux;
		rc.uy = q * uy;
		rc.uz = q * uz;
		return rc;
	}
	inline friend void operator *= (vec3& p, const ucoord3& c)
	{
		p = p * c;
	}
	inline void operator *= (const ucoord3& c)
	{
		*this = (*this) * c;
	}
	inline void operator *= (const quaternion& q)
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
				pl_prj(v - c.uz * v.dot(c.uz), c.ux, c.uy) / c.s.x, \
				pl_prj(v - c.ux * v.dot(c.ux), c.uy, c.uz) / c.s.y, \
				pl_prj(v - c.uy * v.dot(c.uy), c.uz, c.ux) / c.s.z)
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
	inline friend void operator /= (vec3& v, const ucoord3& c)
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
	inline void operator /= (const ucoord3& c)
	{
		*this = (*this) / c;
	}
	inline friend quaternion operator / (const quaternion& q, const ucoord3& c)
	{
		return q * c.toquat().conjcopy();
	}
	inline ucoord3 operator / (const quaternion& q) const
	{
		return (*this) * q.conjcopy();
	}
	inline void operator /= (const quaternion& q)
	{
		*this = (*this) / q;
	}
	// oper(//) = C1^-1 * C2
	inline ucoord3 operator % (const ucoord3& c) const
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
		c.ux = vec3::lerp(vec3::UX, c.ux, v.x); c.ux.norm();
		c.uy = vec3::lerp(vec3::UY, c.uy, v.y); c.uy.norm();
		c.uz = vec3::lerp(vec3::UZ, c.uz, v.z); c.uz.norm();

		return c;
	}
	inline void operator ^= (const vec3& v)
	{
		(*this) = (*this) ^ v;
	}
	ucoord3 operator ^ (real f) const
	{
		/*
		ucoord3 c = *this;
		c.ux = lerp(vec3::UX, c.ux, f); c.ux.norm();
		c.uy = lerp(vec3::UY, c.uy, f); c.uy.norm();
		c.uz = lerp(vec3::UZ, c.uz, f); c.uz.norm();
		*/
		// 四元数法
		return ucoord3((*this).toquat() ^ f);
	}
	inline void operator ^= (real f)
	{
		(*this) = (*this) ^ f;
	}
	// 转置(坐标轴交换）
	void transpose()
	{
		vec3 _ux = vec3(ux.x, uy.x, uz.x);
		vec3 _uy = vec3(ux.y, uy.y, uz.y);
		vec3 _uz = vec3(ux.z, uy.z, uz.z);
		ux = _ux; uy = _uy; uz = _uz;
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
	inline void reverse()
	{
		(*this) = ONE / (*this);
	}
	inline ucoord3 reversed() const
	{
		return ONE / (*this);
	}
	// 翻转
	inline void flipX()
	{
		ux = -ux;
	}
	inline void flipY()
	{
		uy = -uy;
	}
	inline void flipZ()
	{
		uz = -uz;
	}
	void rot(real ang, const vec3& ax)
	{
		ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);
	}
	inline vec3 dir() const
	{
		return (ux + uy + uz).normalized();
	}
	// 归一化
	void norm()
	{
		ux.norm();
		uy.norm();
		uz.norm();
	}
	ucoord3 normcopy()
	{
		ucoord3 c = *this;
		c.ux.norm();
		c.uy.norm();
		c.uz.norm();
		return c;
	}
	// 本征向量（坐标系作为旋转变换时候的特征）
	inline vec3 eigenvec() const
	{
		return toquat().axis();
	}
	inline real dot(const vec3& v) const
	{
		return v.dot(ux) + v.dot(uy) + v.dot(uz);
	}
	inline real dot(const ucoord3& c) const
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
	// 坐标系到欧拉角(pit,yaw,roll)
	vec3 coord2eulers() const
	{
		real c_eps = 1e-5;

		const ucoord3& rm = *this;
		float sy = sqrt(rm.ux.x * rm.ux.x + rm.uy.x * rm.uy.x);
		bool singular = sy < c_eps;

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

	// 转化为四元数
	quat toquat() const
	{
		const ucoord3& rm = *this;
		float trace = rm.ux.x + rm.uy.y + rm.uz.z;
		float w, x, y, z;

		if (trace > 0)
		{
			float s = 0.5f / sqrt(trace + 1.0f);
			w = 0.25f / s;
			x = (rm.uy.z - rm.uz.y) * s;
			y = (rm.uz.x - rm.ux.z) * s;
			z = (rm.ux.y - rm.uy.x) * s;
		}
		else if (rm.ux.x > rm.uy.y && rm.ux.x > rm.uz.z)
		{
			float s = 2.0f * sqrt(1.0f + rm.ux.x - rm.uy.y - rm.uz.z);
			w = (rm.uy.z - rm.uz.y) / s;
			x = 0.25f * s;
			y = (rm.uy.x + rm.ux.y) / s;
			z = (rm.uz.x + rm.ux.z) / s;
		}
		else if (rm.uy.y > rm.uz.z)
		{
			float s = 2.0f * sqrt(1.0f + rm.uy.y - rm.ux.x - rm.uz.z);
			w = (rm.uz.x - rm.ux.z) / s;
			x = (rm.uy.x + rm.ux.y) / s;
			y = 0.25f * s;
			z = (rm.uz.y + rm.uy.z) / s;
		}
		else
		{
			float s = 2.0f * sqrt(1.0f + rm.uz.z - rm.ux.x - rm.uy.y);
			w = (rm.ux.y - rm.uy.x) / s;
			x = (rm.uz.x + rm.ux.z) / s;
			y = (rm.uz.y + rm.uy.z) / s;
			z = 0.25f * s;
		}
		return quat(w, x, y, z);
	}
	// 梯度坐标系 = 梯度 X 切空间
	// 相当于一阶坐标系的导数
	// C2 = UG * C1
	static inline ucoord3 ugrad(const ucoord3& c1, const ucoord3& c2)
	{
		return c1.reversed() * c2;
	}
	void dump(const std::string& name = "") const
	{
		PRINT("----" << name << "---");
		PRINTVEC3(ux);
		PRINTVEC3(uy);
		PRINTVEC3(uz);
	}
};
#ifdef PMDLL
const ucoord3 ucoord3::ZERO = {};
const ucoord3 ucoord3::ONE = {};
#endif

// ******************************************************************
//  |/_
// VC     3d Rotation & Scaling Coordinate System
// ******************************************************************
struct vcoord3 : ucoord3
{
	static const vcoord3 ZERO;
	static const vcoord3 ONE;

	vec3 s = vec3::ONE;		// 缩放

	vcoord3() {}
	vcoord3(const vcoord3& c) : ucoord3(c.ux, c.uy, c.uz), s(c.s)	{}
	vcoord3(const ucoord3& c) : ucoord3(c){}
	vcoord3(const ucoord3& c, const vec3& _s) : ucoord3(c), s(_s){}
	vcoord3( const vec3& _ux, const vec3& _uy, const vec3& _uz, const vec3& _s) : ucoord3(_ux, _uy, _uz), s(_s){ }
	vcoord3( const vec3& _ux, const vec3& _uy, const vec3& _uz) : ucoord3(_ux, _uy, _uz){}
	vcoord3(const quaternion& q) : ucoord3(q){}
	vcoord3(const quaternion& q, const vec3& _s) : ucoord3(q), s(_s) {}

	// 乘法：在坐标系下定义一个向量
	inline friend vec3 operator * (const vec3& p, const vcoord3& c)
	{
		return c.ux * (c.s.x * p.x) + c.uy * (c.s.y * p.y) + c.uz * (c.s.z * p.z);
	}
	inline friend void operator *= (vec3& p, const vcoord3& c)
	{
		p = p * c;
	}
	inline vcoord3 operator * (const vec3& v) const
	{
		return (*this) * vcoord3(vec3::UX * v.x, vec3::UY * v.y, vec3::UZ * v.z);
	}
	inline void operator *= (const vec3& v)
	{
		*this = (*this) * v;
	}
	inline friend real operator * (const real& s, const vcoord3& c)
	{
		return s * ((c.s.x + c.s.y + c.s.z) / 3.0);
	}
	vcoord3 operator * (real s) const
	{
		vcoord3 c = *this;
		{// C*S 缩放乘法
			c.s.x *= s; c.s.y *= s; c.s.z *= s;
		}
		return c;
	}
	void operator *= (real s)
	{
		*this = (*this) * s;
	}
	vcoord3 operator * (const vcoord3& c) const
	{// Cchild * Cparent * ...
		vcoord3 rc;
		rc.ux = ux.x * c.ux + ux.y * c.uy + ux.z * c.uz;
		rc.uy = uy.x * c.ux + uy.y * c.uy + uy.z * c.uz;
		rc.uz = uz.x * c.ux + uz.y * c.uy + uz.z * c.uz;

		rc.s = s * c.s;
		return rc;
	}
	inline void operator *= (const vcoord3& c)
	{
		*this = (*this) * c;
	}
	vcoord3 operator * (const quaternion& q) const
	{
		vcoord3 rc = *this;
		rc.ux = q * ux;
		rc.uy = q * uy;
		rc.uz = q * uz;
		return rc;
	}
	inline void operator *= (const quaternion& q)
	{
		*this = (*this) * q;
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
	inline friend vec3 operator / (const vec3& v, const vcoord3& c)
	{
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
	inline friend void operator /= (vec3& p, const vcoord3& c)
	{
		p = p / c;
	}
	inline vcoord3 operator / (const vec3& v) const
	{
		return (*this) / vcoord3(vec3::UX * v.x, vec3::UY * v.y, vec3::UZ * v.z);
	}
	inline void operator /= (const vec3& v)
	{
		*this = (*this) / v;
	}

	vcoord3 operator / (real s) const
	{// C/S 缩放除法
		vcoord3 c = *this;
		c.s /= s;
		return c;
	}
	inline void operator /= (real s)
	{
		*this = (*this) / s;
	}
	// oper(/) = C1 * C2^-1
	vcoord3 operator / (const vcoord3& c) const
	{
		vcoord3 rc;
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
		return rc;
	}
	void operator /= (const vcoord3& c)
	{
		*this = (*this) / c;
	}
	vcoord3 operator / (const quaternion& q) const
	{
		return (*this) * q.conjcopy();
	}
	void operator /= (const quaternion& q)
	{
		*this = (*this) / q;
	}

	// 归一化
	void norm()
	{
		ucoord3::norm();
		s.norm();
	}
	vcoord3 normcopy()
	{
		vcoord3 c = *this;
		c.norm();
		return c;
	}
};
#ifdef PMDLL
const vcoord3 vcoord3::ZERO = {};
const vcoord3 vcoord3::ONE = {};
#endif

// ******************************************************************
//  |/_
// C     3d Coordinate System
// ******************************************************************
struct coord3 : vcoord3
{
	static const coord3 ZERO;
	static const coord3 ONE;

	vec3 o;					// 原点

	coord3() {}
	coord3(const coord3& c) : vcoord3(c.ux, c.uy, c.uz, c.s), o(c.o){}
	coord3(const ucoord3& c) : vcoord3(c){}
	coord3(const vec3& _o, const vec3& _s, const vec3& _ux, const vec3& _uy, const vec3& _uz) : vcoord3(_ux, _uy, _uz, _s), o(_o){}
	coord3(const vec3& _o, const vec3& _ux, const vec3& _uy, const vec3& _uz) : vcoord3(_ux, _uy, _uz), o(_o){ }
	coord3(const vec3& _ux, const vec3& _uy, const vec3& _uz) : vcoord3(_ux, _uy, _uz){}
	coord3(const vec3& _ux, const vec3& _uy) : vcoord3(_ux, _uy, ux.cross(uy)){}
	coord3(const vec3& _p) : o(_p){}
	coord3(const ucoord3& c,const vec3& _o) : vcoord3(c), o(_o){}
	coord3(const vec3& _o, const ucoord3& c) : vcoord3(c), o(_o){}
	coord3(const ucoord3& c, const vec3& _s, const vec3& _o) : vcoord3(c, _s), o(_o){}
	coord3(const vec3& _o, const vec3& _s, const ucoord3& c) : vcoord3(c, _s), o(_o) {}
	coord3(real ang, const vec3& ax)
	{
		ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);
	}
	coord3(real x, real y, real z)
	{
		o = vec3(x,y,z);
	}	
	coord3(real x, real y, real z, real rx, real ry, real rz)
	{
		quaternion q(rx, ry, rz);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
		o = vec3(x, y, z);
	}
	coord3(const quaternion& q) : vcoord3(q){}
	coord3(const vec3& p, const quaternion& q, const vec3& _s = vec3::ONE) : vcoord3(q, _s), o(p){}

	// 移动差
	inline void fromvecsT(const vec3& v1, const vec3& v2)
	{
		o = v2 - v1;
	}
	inline operator quat () const
	{
		return toquat();
	}
	inline operator vec3 () const
	{
		return o;
	}
	inline vec3 VX() const { return ux * s.x; }
	inline vec3 VY() const { return uy * s.y; }
	inline vec3 VZ() const { return uz * s.z; }

	inline vec3 X() const { return ux * s.x + vec3::UX * o.x; }
	inline vec3 Y() const { return uy * s.y + vec3::UY * o.y; }
	inline vec3 Z() const { return uz * s.z + vec3::UZ * o.z; }

	// 旋转坐标系
	inline const ucoord3& ucoord() const
	{
		return static_cast<const ucoord3&>(*this);
	}
	inline void ucoord(const ucoord3& ucd)
	{
		ux = ucd.ux; uy = ucd.uy; uz = ucd.uz;
	}
	inline const ucoord3& ucrd() const
	{
		return static_cast<const ucoord3&>(*this);
	}
	inline void ucrd(const ucoord3& ucd)
	{
		ux = ucd.ux; uy = ucd.uy; uz = ucd.uz;
	}
	inline const ucoord3& UC() const
	{
		return static_cast<const ucoord3&>(*this);
	}
	inline ucoord3 UC(const ucoord3& ucd)
	{
		ux = ucd.ux; uy = ucd.uy; uz = ucd.uz;
	}
	inline const ucoord3& RC() const
	{
		return static_cast<const ucoord3&>(*this);
	}
	inline ucoord3 RC(const ucoord3& ucd)
	{
		ux = ucd.ux; uy = ucd.uy; uz = ucd.uz;
	}
	// 向量坐标系 = 方向 X 缩放
	inline vcoord3 vcoord() const
	{
		return static_cast<const vcoord3&>(*this);
	}
	inline coord3 vcrd() const
	{
		return static_cast<const vcoord3&>(*this);
	}
	inline coord3 VC() const
	{
		return static_cast<const vcoord3&>(*this);
	}
	// 姿态
	inline coord3 pose()
	{
		return { ucoord(), vec3::ONE, o };
	}
	// 位置
	inline vec3 pos() const
	{
		return o;
	}
	// 总向量
	inline vec3 sumvec() const
	{
		return o + VX() + VY() + VZ();
	}
	// 向量
	inline vec3 tovec() const
	{
		return ux * s.x + uy * s.y + uz * s.z;
	}
	inline quaternion Q() const
	{
		return toquat();
	}
	inline void Q(const quaternion& q)
	{
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}

	coord3 operator = (const coord3& c)
	{
		o = c.o;
		s = c.s;
		ux = c.ux; uy = c.uy; uz = c.uz;
		return (*this);
	}
	inline bool is_same_dirs(const coord3& c) const
	{
		return ux == c.ux && uy == c.uy && uz == c.uz;
	}
	inline bool operator == (const coord3& c) const
	{
		return o == c.o && s == c.s && is_same_dirs(c);
	}
	inline bool operator != (const coord3& c) const
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
	inline void operator += (const coord3& c)
	{
		*this = (*this) + c;
	}
	inline coord3 operator + (const vec3& v) const
	{
		coord3 c = (*this); c.o += v;
		return c;
	}
	inline void operator += (const vec3& v)
	{
		*this = *this + v;
	}
	inline friend vec3 operator + (const vec3& p, const coord3& c)
	{
		return p + c.o;
	}
	inline friend void operator += (vec3& p, const coord3& c)
	{
		p = p + c;
	}
	inline friend vec3 operator - (const vec3& p, const coord3& c)
	{
		return p - c.o;
	}
	inline friend void operator -= (vec3& p, const coord3& c)
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
	inline coord3 operator - (const vec3& v) const
	{
		coord3 c = (*this); c.o -= v;
		return c;
	}
	inline void operator -= (const vec3& v)
	{
		*this = *this - v;
	}

	// 乘法：在坐标系下定义一个向量
	inline friend vec3 operator * (const vec3& p, const coord3& c)
	{
		return c.ux * (c.s.x * p.x) + c.uy * (c.s.y * p.y) + c.uz * (c.s.z * p.z) + c.o;
	}
	inline friend void operator *= (vec3& p, const coord3& c)
	{
		p = p * c;
	}
	inline coord3 operator * (const vec3& v) const
	{
		return (*this) * coord3(vec3::UX * v.x, vec3::UY * v.y, vec3::UZ * v.z);
	}
	inline void operator *= (const vec3& v)
	{
		*this = (*this) * v;
	}
	inline friend real operator * (const real& s, const coord3& c)
	{
		return s * ((c.s.x + c.s.y + c.s.z) / 3.0);
	}
	coord3 operator * (real s) const
	{
		coord3 c = *this;
		{// C*S 缩放乘法
			c.s.x *= s; c.s.y *= s; c.s.z *= s;
		}
		{// C*S 移动乘法
			c.o.x *= s; c.o.y *= s; c.o.z *= s;
		}
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
		rc.o = c.o + o.x * c.s.x * c.ux + o.y * c.s.y * c.uy + o.z * c.s.z * c.uz;
		return rc;
	}
	inline void operator *= (const coord3& c)
	{
		*this = (*this) * c;
	}
	coord3 operator * (const quaternion& q) const
	{
		coord3 rc = *this;
		rc.ux = q * ux;
		rc.uy = q * uy;
		rc.uz = q * uz;
		rc.o = q * rc.o;
		return rc;
	}
	inline void operator *= (const quaternion& q)
	{
		*this = (*this) * q;
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
	inline friend vec3 operator / (const vec3& p, const coord3& c)
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
	inline friend void operator /= (vec3& p, const coord3& c)
	{
		p = p / c;
	}
	inline coord3 operator / (const vec3& v) const
	{
		return (*this) / coord3(vec3::UX * v.x, vec3::UY * v.y, vec3::UZ * v.z);
	}
	inline void operator /= (const vec3& v)
	{
		*this = (*this) / v;
	}

	coord3 operator / (real s) const
	{// C/S 缩放除法
		coord3 c = *this;
		c.s /= s;
		c.o /= s;
		return c;
	}
	inline void operator /= (real s)
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
	coord3 operator / (const quaternion& q) const
	{
		return (*this) * q.conjcopy();
	}
	void operator /= (const quaternion& q)
	{
		*this = (*this) / q;
	}
	// oper(//) = C1^-1 * C2
	coord3 operator % (const coord3& c) const
	{
		return (*this).reversed() * c;
	}
	coord3 operator ^ (const vec3& v) const
	{
		coord3 c = *this;
		c.ux = vec3::lerp(vec3::UX, c.ux, v.x); c.ux.norm();
		c.uy = vec3::lerp(vec3::UY, c.uy, v.y); c.uy.norm();
		c.uz = vec3::lerp(vec3::UZ, c.uz, v.z); c.uz.norm();

		c.s = vec3::lerp(vec3::ONE, c.s, v);
		c.o = vec3::lerp(vec3::ZERO, c.o, v);

		return c;
	}
	coord3 operator ^ (real t) const
	{
		ucoord3 uc = ucoord();
		uc ^= t;
		/*c.ux = vec3::lerp(vec3::UX, c.ux, t); c.ux.norm();
		c.uy = vec3::lerp(vec3::UY, c.uy, t); c.uy.norm();
		c.uz = vec3::lerp(vec3::UZ, c.uz, t); c.uz.norm();*/

		vec3 s = vec3::lerp(vec3::ONE, s, t);
		vec3 o = vec3::lerp(vec3::ZERO, o, t);

		return coord3(uc, s, o);
	}
	// 归一化
#define ISZERO(a) (fabs(a) < 1e-10)
	void norm(bool bscl = true)
	{
		s.x = ux.len(); if (!ISZERO(s.x)) ux /= s.x;
		s.y = uy.len(); if (!ISZERO(s.y)) uy /= s.y;
		s.z = uz.len(); if (!ISZERO(s.z)) uz /= s.z;
		if (!bscl)
			s = vec3::ONE;
	}
	coord3 normcopy(bool bscl = true) const
	{
		coord3 c = *this;
		c.norm();
		return c;
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
		return coord3(
			VX().cross(v),
			VY().cross(v),
			VZ().cross(v)
		);
	}
	// 梯度坐标系 = 梯度 X 切空间
	// 相当于一阶坐标系的导数
	// C2 = UG * C1
	// V2 - V1 = G * V1 = (UG - ONE) * V1
	// G = UG - ONE
	static coord3 grad(const coord3& c1, const coord3& c2)
	{
		return c2 / c1 - ONE;
	}
	static coord3 G(const coord3& c1, const coord3& c2)
	{
		return c2 / c1 - ONE;
	}
	std::string serialise()
	{
		vec3 eu = coord2eulers();
		return o.serialise() + "," + eu.serialise();
	}
	void dump(const std::string& name = "") const
	{
		PRINT("----" << name << "---");
		PRINTVEC3(ux);
		PRINTVEC3(uy);
		PRINTVEC3(uz);
		PRINTVEC3(s);
		PRINTVEC3(o);
	}

	/// 便捷函数 ///
	void rot(real angle, const vec3& ax)
	{
		quaternion q(angle, ax);
		ux *= q;
		uy *= q;
		uz *= q;
	}
	void rot(const quaternion& q)
	{
		ux = q * ux;
		uy = q * uy;
		uz = q * uz;
	}
	coord3 rotcopy(real ang, const vec3& ax) const
	{
		coord3 c = *this;
		c.ux.rot(ang, ax);
		c.uy.rot(ang, ax);
		c.uz.rot(ang, ax);
		return c;
	}
	coord3 rotcopy(const quaternion& q) const
	{
		coord3 c = *this;
		c.ux *= q;
		c.uy *= q;
		c.uz *= q;
		return c;
	}
	void uxto(const vec3& _ux)
	{
		*this *= quat(ux, _ux);
	}
	void uyto(const vec3& _uy)
	{
		*this *= quat(uy, _uy);
	}
	void uzto(const vec3& _uz)
	{
		*this *= quat(uz, _uz);
	}
	void moveto(const coord3& c, real alpha = 1)
	{
		o = vec3::lerp(o, c.o, alpha);
	}
};
#ifdef PMDLL
const coord3 coord3::ZERO = {};
const coord3 coord3::ONE = {};
#endif
