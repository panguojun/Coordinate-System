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
*   without changing) under V. Observing V in a curved coordinate system, V is different at points.
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

//#define	NON_UNIFORM_SCALE
// ********************************************************************************************
//  |/_
// UC     3d Rotation Coordinate System
// ********************************************************************************************
struct ucoord3
{
	static const ucoord3 ONE;

	vec3 ux = vec3::UX;		// basis 单位化基向量
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
		quaternion q(ang, ax);
		ux = q * ux;
		uy = q * uy;
		uz = q * uz;

		/*ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);*/
	}
    ucoord3(real pit, real yaw, real rol)
	{
		quaternion q(pit, yaw, rol);
		ux = q * ux;
		uy = q * uy;
		uz = q * uz;
	}
	ucoord3(const quaternion& q)
	{
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}

	// uy方向 推测ux,uz
	void fromquat(const quaternion& q)
	{
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}
	void fromuy(const vec3& _uy)
	{
		quat q; q.fromvectors(uy, _uy);
		fromquat(q);
	}
	// 引用四元数的欧拉角转化
	void frompyr(real pit, real yaw, real rol)
	{
		fromquat({ pit, yaw, rol });
	}
	void frompyr(const vec3& pyr)
	{
		fromquat(quaternion(pyr.x, pyr.y, pyr.z));
	}
	vec3 topyr() const
	{
		return Q().toeulers();
	}
	// 坐标系的欧拉角转化
	vec3 toeulers() const
	{
		return coord2eulers();
	}
	// 旋转差
	void from_vecs_R(const vec3& v1, const vec3& v2)
	{
		vec3  v = v1.cross(v2);
		float c = v1.dot(v2);
		float k = 1.0 / (1.0 + c);
		
		ux = { v.x * v.x * k + c,      v.y * v.x * k - v.z,    v.z * v.x * k + v.y };
		uy = { v.x * v.y * k + v.z,    v.y * v.y * k + c,      v.z * v.y * k - v.x };
		uz = { v.x * v.z * k - v.y,    v.y * v.z * k + v.x,    v.z * v.z * k + c   };
	}
	// 轴，向量1，2
	void from_ax_vecs(const vec3& ax, const vec3& v1, const vec3& v2)
	{
		vec3 pv1 = v1.crossdot(ax);
		vec3 pv2 = v2.crossdot(ax);
		real ang = acos(pv1.dot(pv2));
		quaternion q; q.ang_axis(ang, ax);
		fromquat(q);
	}
	bool same_dirs(const ucoord3& c) const
	{
		return ux == c.ux && uy == c.uy && uz == c.uz;
	}
	bool operator == (const ucoord3& c) const
	{
		return same_dirs(c);
	}
	bool operator != (const ucoord3& c) const
	{
		return !same_dirs(c);
	}
	vec3 operator[] (int index) const
	{
		if(index == 0)
			return ux;
		else if (index == 1)
			return uy;
		else if (index == 2)
			return uz;

		return vec3::ZERO;
	}

	// 乘法：在坐标系下定义一个向量，或者向量向父空间还原
	friend vec3 operator * (const vec3& v, const ucoord3& c)
	{
		return c.ux * (v.x) + c.uy * (v.y) + c.uz * (v.z);
	}
	ucoord3 operator * (const ucoord3& c) const
	{// C_child * C_parent * ...
		ucoord3 rc;
		rc.ux = ux.x * c.ux + ux.y * c.uy + ux.z * c.uz;
		rc.uy = uy.x * c.ux + uy.y * c.uy + uy.z * c.uz;
		rc.uz = uz.x * c.ux + uz.y * c.uy + uz.z * c.uz;

		return rc;
	}
	friend quaternion operator * (const quaternion& q, const ucoord3& c)
	{
		return  q * c.toquat();
	}
	ucoord3 operator * (const quaternion& q) const
	{
		ucoord3 rc;
		rc.ux = q * ux;
		rc.uy = q * uy;
		rc.uz = q * uz;
		return rc;
	}
	friend void operator *= (vec3& v, const ucoord3& c)
	{
		v = v * c;
	}
	void operator *= (const ucoord3& c)
	{
		*this = (*this) * c;
	}
	void operator *= (const quaternion& q)
	{
		ux = q * ux;
		uy = q * uy;
		uz = q * uz;
	}
	// 除法：向量向坐标系投影
    DEVICE_CALLABLE friend vec3 operator/(const vec3& v, const ucoord3& c)
	{
		return vec3(v.dot(c.ux), v.dot(c.uy), v.dot(c.uz));
	}
    DEVICE_CALLABLE friend void operator/=(vec3& v, const ucoord3& c)
	{
		v = v / c;
	}
	// oper(/) = C1 * C2^-1
    DEVICE_CALLABLE ucoord3 operator/(const ucoord3& c) const
	{
		ucoord3 rc;
		rc.ux = vec3(ux.dot(c.ux), ux.dot(c.uy), ux.dot(c.uz));
		rc.uy = vec3(uy.dot(c.ux), uy.dot(c.uy), uy.dot(c.uz));
		rc.uz = vec3(uz.dot(c.ux), uz.dot(c.uy), uz.dot(c.uz));
		return rc;
	}
    DEVICE_CALLABLE void operator/=(const ucoord3& c)
	{
		*this = (*this) / c;
	}
	friend quaternion operator / (const quaternion& q, const ucoord3& c)
	{
		return q * c.toquat().conjcopy();
	}
	ucoord3 operator / (const quaternion& q) const
	{
		return (*this) * q.conjcopy();
	}
	void operator /= (const quaternion& q)
	{
		*this = (*this) / q;
	}
	// oper(\) = C1^-1 * C2
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
		c.ux = vec3::lerp(vec3::UX, c.ux, v.x); c.ux.norm();
		c.uy = vec3::lerp(vec3::UY, c.uy, v.y); c.uy.norm();
		c.uz = vec3::lerp(vec3::UZ, c.uz, v.z); c.uz.norm();

		return c;
	}
	void operator ^= (const vec3& v)
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
	void operator ^= (real f)
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
		quaternion q(ang, ax);
		ux = q * ux;
		uy = q * uy;
		uz = q * uz;
	}
	vec3 dir() const
	{
		return (ux + uy + uz).normalized();
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
	ucoord3 eulers2coord(const vec3& eulers)
	{
		float x = eulers.x;
		float y = eulers.y;
		float z = eulers.z;

		float cx = cos(x);
		float sx = sin(x);
		float cy = cos(y);
		float sy = sin(y);
		float cz = cos(z);
		float sz = sin(z);

		ucoord3 result;
		result.ux.x = cy * cz;
		result.ux.y = -cy * sz;
		result.ux.z = sy;

		result.uy.x = cx * sz + sx * sy * cz;
		result.uy.y = cx * cz - sx * sy * sz;
		result.uy.z = -sx * cy;

		result.uz.x = sx * sz - cx * sy * cz;
		result.uz.y = sx * cz + cx * sy * sz;
		result.uz.z = cx * cy;

		return result;
	}

	// 转化为四元数, 注意四元数的乘法顺序
	quat toquat() const
	{
		return quat(coord2eulers());
	}
	quat Q() const
	{
		return quat(coord2eulers());
	}
	// 梯度坐标系 = 梯度 X 切空间
	// 相当于一阶坐标系的导数
	// C2 = UG * C1
	static ucoord3 ugrad(const ucoord3& c1, const ucoord3& c2)
	{
		return c1.reversed() * c2;
	}
	static ucoord3 R(const ucoord3& c1, const ucoord3& c2)
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

	// 方便函数, 注意 angle(u,_u) != +/-PI
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
	ucoord3 uxtoed(const vec3& _ux) const
	{
		return (*this) * quat(ux, _ux);
	}
	ucoord3 uytoed(const vec3& _uy) const
	{
		return (*this) * quat(uy, _uy);
	}
	ucoord3 uztoed(const vec3& _uz) const
	{
		return (*this) * quat(uz, _uz);
	}
};
#if defined(PMDLL) || !defined(PM_IMPLEMENTED)
const ucoord3 ucoord3::ONE = {};
#endif

// ******************************************************************
//  |/_
// VC     3d Rotation & Scaling Coordinate System
// ******************************************************************
struct vcoord3 : ucoord3
{
	static const vcoord3 ONE;

	vec3 s = vec3::ONE;		// 缩放

	vcoord3() {}
	vcoord3(const ucoord3& c) : ucoord3(c){}
	vcoord3(const ucoord3& c, const vec3& _s) : ucoord3(c), s(_s){}
	vcoord3( const vec3& _ux, const vec3& _uy, const vec3& _uz, const vec3& _s) : ucoord3(_ux, _uy, _uz), s(_s){ }
	vcoord3( const vec3& _ux, const vec3& _uy, const vec3& _uz) : ucoord3(_ux, _uy, _uz){}
	vcoord3(const quaternion& q) : ucoord3(q){}
	vcoord3(const quaternion& q, const vec3& _s) : ucoord3(q), s(_s) {}
	vcoord3(const vec3& _s) : s(_s) {}

	vec3 VX() const { return ux * s.x; }
	vec3 VY() const { return uy * s.y; }
	vec3 VZ() const { return uz * s.z; }

	const ucoord3& R() const
	{
		return static_cast<const ucoord3&>(*this);
	}
	const ucoord3& UC() const
	{
		return static_cast<const ucoord3&>(*this);
	}
	void UC(const ucoord3& ucd)
	{
		ux = ucd.ux; uy = ucd.uy; uz = ucd.uz;
	}

	// 乘法：在坐标系下定义一个向量
	friend vec3 operator * (const vec3& p, const vcoord3& c)
	{
		return c.ux * (c.s.x * p.x) + c.uy * (c.s.y * p.y) + c.uz * (c.s.z * p.z);
	}
	friend void operator *= (vec3& p, const vcoord3& c)
	{
		p = p * c;
	}
	vcoord3 operator * (const vec3& v) const
	{
		return (*this) * vcoord3(vec3::UX * v.x, vec3::UY * v.y, vec3::UZ * v.z);
	}
	void operator *= (const vec3& v)
	{
		*this = (*this) * v;
	}
	friend real operator * (const real& s, const vcoord3& c)
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
#ifdef	NON_UNIFORM_SCALE
		rc.ux = (ux.x * s.x) * (c.ux * c.s.x) + (ux.y * s.x) * (c.uy * c.s.y) + (ux.z * s.x) * (c.uz * c.s.z);
		rc.uy = (uy.x * s.y) * (c.ux * c.s.x) + (uy.y * s.y) * (c.uy * c.s.y) + (uy.z * s.y) * (c.uz * c.s.z);
		rc.uz = (uz.x * s.z) * (c.ux * c.s.x) + (uz.y * s.z) * (c.uy * c.s.y) + (uz.z * s.z) * (c.uz * c.s.z);
		rc.norm();
#else
		rc = ucoord3::operator*(c);
		rc.s = s * c.s;
#endif
		return rc;
	}
	void operator *= (const vcoord3& c)
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
	void operator *= (const quaternion& q)
	{
		*this = (*this) * q;
	}
	
	// 除法：向量向坐标系投影 注意：要保证ux,uy,uz是单位向量！
	friend vec3 operator / (const vec3& v, const vcoord3& c)
	{
		return vec3(v.dot(c.ux) / c.s.x, v.dot(c.uy) / c.s.y, v.dot(c.uz) / c.s.z);
	}
	friend void operator /= (vec3& p, const vcoord3& c)
	{
		p = p / c;
	}
	vcoord3 operator / (const vec3& v) const
	{
		return (*this) / vcoord3(vec3::UX * v.x, vec3::UY * v.y, vec3::UZ * v.z);
	}
	void operator /= (const vec3& v)
	{
		*this = (*this) / v;
	}

	vcoord3 operator / (real s) const
	{// C/S 缩放除法
		vcoord3 c = *this;
		c.s /= s;
		return c;
	}
	void operator /= (real s)
	{
		*this = (*this) / s;
	}
	// oper(/) = C1 * C2^-1
	vcoord3 operator / (const vcoord3& c) const
	{
		vcoord3 rc;
#ifdef	NON_UNIFORM_SCALE
		vec3 vx = VX();
		vec3 vy = VY();
		vec3 vz = VZ();

		vec3 cvx = c.ux / c.s.x;
		vec3 cvy = c.uy / c.s.y;
		vec3 cvz = c.uz / c.s.z;

		rc.ux = vec3(vx.dot(cvx), vx.dot(cvy), vx.dot(cvz));
		rc.uy = vec3(vy.dot(cvx), vy.dot(cvy), vy.dot(cvz));
		rc.uz = vec3(vz.dot(cvx), vz.dot(cvy), vz.dot(cvz));

		rc.norm();
#else
		rc = ucoord3::operator/(c);
		rc.s = s / c.s;
#endif
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
	void norm(bool bscl = true)
	{
		s.x = ux.len(); if (!ISZERO(s.x)) ux /= s.x;
		s.y = uy.len(); if (!ISZERO(s.y)) uy /= s.y;
		s.z = uz.len(); if (!ISZERO(s.z)) uz /= s.z;
		if (!bscl)
			s = vec3::ONE;
	}
	vcoord3 normcopy(bool bscl = true) const
	{
		vcoord3 c = *this;
		c.norm(bscl);
		return c;
	}

	// 倒置
	void reverse()
	{
		(*this) = ONE / (*this);
	}
	vcoord3 reversed() const
	{
		return ONE / (*this);
	}
	// DOT
	real dot(const vec3& v) const
	{
		return v.dot(ux) * s.x + v.dot(uy) * s.y + v.dot(uz) * s.z;
	}
	real dot(const vcoord3& c) const
	{
		return c.VX().dot(VX()) + c.VY().dot(VY()) + c.VZ().dot(VZ());
	}

	void dump(const std::string& name = "") const
	{
		PRINT("----" << name << "---");
		PRINTVEC3(ux);
		PRINTVEC3(uy);
		PRINTVEC3(uz);
		PRINTVEC3(s);
	}
};
#if  defined(PMDLL) || !defined(PM_IMPLEMENTED)
const vcoord3 vcoord3::ONE	= { };
#endif

// ******************************************************************
//  |/_
// C     3d Coordinate System
// ******************************************************************
struct coord3 : vcoord3
{
	static const coord3 ZERO;
	static const coord3 ONE;

	vec3 o;		// 原点

	DEVICE_CALLABLE coord3() {}
	DEVICE_CALLABLE coord3(const ucoord3& uc) : vcoord3(uc){}
	DEVICE_CALLABLE coord3(const vcoord3& vc) : vcoord3(vc) {}
	DEVICE_CALLABLE coord3(const vec3& _o,  const vec3& _s, const vec3& _ux, const vec3& _uy, const vec3& _uz) : vcoord3(_ux, _uy, _uz, _s), o(_o){}
	DEVICE_CALLABLE coord3(const vec3& _o,  const vec3& _ux, const vec3& _uy, const vec3& _uz) : vcoord3(_ux, _uy, _uz), o(_o){ }
	DEVICE_CALLABLE coord3(const vec3& _ux, const vec3& _uy, const vec3& _uz) : vcoord3(_ux, _uy, _uz){}
	DEVICE_CALLABLE coord3(const vec3& _ux, const vec3& _uy) : vcoord3(_ux, _uy, ux.cross(uy)){}
	DEVICE_CALLABLE coord3(const vec3& _p) : o(_p){}
	DEVICE_CALLABLE coord3(const ucoord3& c,const vec3& _o) : vcoord3(c), o(_o){}
	DEVICE_CALLABLE coord3(const vec3& _o,  const ucoord3& c) : vcoord3(c), o(_o){}
	DEVICE_CALLABLE coord3(const ucoord3& c,const vec3& _s, const vec3& _o) : vcoord3(c, _s), o(_o){}
	DEVICE_CALLABLE coord3(const vec3& _o,  const vec3& _s, const ucoord3& c) : vcoord3(c, _s), o(_o) {}
	DEVICE_CALLABLE coord3(real ang, const vec3& ax) : vcoord3(quaternion(ang, ax)) {}
    DEVICE_CALLABLE coord3(real x, real y, real z) : o(x, y, z) {}
    DEVICE_CALLABLE coord3(const quaternion& q) : vcoord3(q) {}
    DEVICE_CALLABLE coord3(const vec3& p, const quaternion& q, const vec3& _s = vec3::ONE) : vcoord3(q, _s), o(p) {}
	DEVICE_CALLABLE coord3(real x, real y, real z, real rx, real ry, real rz)
	{
		quaternion q(rx, ry, rz);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
		o = vec3(x, y, z);
	}
	DEVICE_CALLABLE operator quat() const
	{
		return toquat();
	}
    DEVICE_CALLABLE operator vec3() const
	{
		return o;
	}
    DEVICE_CALLABLE vec3 VX() const { return ux * s.x; }
	DEVICE_CALLABLE vec3 VY() const { return uy * s.y; }
	DEVICE_CALLABLE vec3 VZ() const { return uz * s.z; }

	DEVICE_CALLABLE vec3 X() const { return ux * s.x + vec3::UX * o.x; }
	DEVICE_CALLABLE vec3 Y() const { return uy * s.y + vec3::UY * o.y; }
	DEVICE_CALLABLE vec3 Z() const { return uz * s.z + vec3::UZ * o.z; }

	// 旋转坐标系
    DEVICE_CALLABLE const ucoord3& ucoord() const
	{
		return static_cast<const ucoord3&>(*this);
	}
    DEVICE_CALLABLE void ucoord(const ucoord3& ucd)
	{
		ux = ucd.ux; uy = ucd.uy; uz = ucd.uz;
	}
    DEVICE_CALLABLE void ucoord(vec3 _ux, vec3 _uy, vec3 _uz)
	{
		ux = _ux; uy = _uy; uz = _uz;
	}
	DEVICE_CALLABLE const ucoord3& R() const
	{
		return static_cast<const ucoord3&>(*this);
	}
    DEVICE_CALLABLE const ucoord3& UC() const
	{
		return static_cast<const ucoord3&>(*this);
	}
    DEVICE_CALLABLE void UC(const ucoord3& ucd)
	{
		ux = ucd.ux; uy = ucd.uy; uz = ucd.uz;
	}
    DEVICE_CALLABLE void UC(vec3 _ux, vec3 _uy, vec3 _uz)
	{
		ux = _ux; uy = _uy; uz = _uz;
	}
	// 向量坐标系 = 方向 X 缩放
    DEVICE_CALLABLE const vcoord3& vcoord() const
	{
		return static_cast<const vcoord3&>(*this);
	}
    DEVICE_CALLABLE const vcoord3& VC() const
	{
		return static_cast<const vcoord3&>(*this);
	}
	// 姿态
    DEVICE_CALLABLE coord3 pose()
	{
		return { ucoord(), vec3::ONE, o };
	}
	// 位置
    DEVICE_CALLABLE vec3 pos() const
	{
		return o;
	}
	// 总向量
    DEVICE_CALLABLE vec3 sumvec() const
	{
		return o + VX() + VY() + VZ();
	}
	// 向量
    DEVICE_CALLABLE vec3 tovec() const
	{
		return ux * s.x + uy * s.y + uz * s.z;
	}
    DEVICE_CALLABLE quaternion Q() const
	{
		return toquat();
	}
    DEVICE_CALLABLE void Q(const quaternion& q)
	{
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}
    DEVICE_CALLABLE coord3 operator=(const coord3& c)
	{
		o = c.o;
		s = c.s;
		ux = c.ux; uy = c.uy; uz = c.uz;
		return (*this);
	}
    DEVICE_CALLABLE bool equal_dirs(const coord3& c) const
	{
		return ux == c.ux && uy == c.uy && uz == c.uz;
	}
    DEVICE_CALLABLE bool operator==(const coord3& c) const
	{
		return o == c.o && s == c.s && equal_dirs(c);
	}
    DEVICE_CALLABLE bool operator!=(const coord3& c) const
	{
		return o != c.o || s != c.s || !equal_dirs(c);
	}
	// +/- 运算
    DEVICE_CALLABLE coord3 operator+(const coord3& c) const
	{
		coord3 rc;
		vec3 _ux = VX() + c.VX();
		vec3 _uy = VY() + c.VY();
		vec3 _uz = VZ() + c.VZ();

		rc.s.x = _ux.len();
		if (!ISZERO(rc.s.x))
		{
			_ux /= rc.s.x;
			rc.ux = _ux;
		}
		rc.s.y = _uy.len();
		if (!ISZERO(rc.s.y))
		{
			_uy /= rc.s.y;
			rc.uy = _uy;
		}
		rc.s.z = _uz.len();
		if (!ISZERO(rc.s.z))
		{
			_uz /= rc.s.z;
			rc.uz = _uz;
		}

		rc.o = o + c.o;
		return rc;
	}
    DEVICE_CALLABLE void operator+=(const coord3& c)
	{
		*this = (*this) + c;
	}
    DEVICE_CALLABLE coord3 operator+(const vec3& v) const
	{
		coord3 c = (*this); c.o += v;
		return c;
	}
    DEVICE_CALLABLE void operator+=(const vec3& v)
	{
		*this = *this + v;
	}
    DEVICE_CALLABLE friend vec3 operator+(const vec3& p, const coord3& c)
	{
		return p + c.o;
	}
    DEVICE_CALLABLE friend void operator+=(vec3& p, const coord3& c)
	{
		p = p + c;
	}
    DEVICE_CALLABLE friend vec3 operator-(const vec3& p, const coord3& c)
	{
		return p - c.o;
	}
    DEVICE_CALLABLE friend void operator-=(vec3& p, const coord3& c)
	{
		p = p - c;
	}
    DEVICE_CALLABLE coord3 operator-(const coord3& c) const
	{
		coord3 rc;
		vec3 _ux = VX() - c.VX();
		vec3 _uy = VY() - c.VY();
		vec3 _uz = VZ() - c.VZ();

		rc.s.x = _ux.len();
		if (!ISZERO(rc.s.x))
		{
			_ux /= rc.s.x;
			rc.ux = _ux;
		}
		rc.s.y = _uy.len();
		if (!ISZERO(rc.s.y))
		{
			_uy /= rc.s.y;
			rc.uy = _uy;
		}
		rc.s.z = _uz.len();
		if (!ISZERO(rc.s.z))
		{
			_uz /= rc.s.z;
			rc.uz = _uz;
		}

		rc.o = o - c.o;
		return rc;
	}
    DEVICE_CALLABLE coord3 operator-(const vec3& v) const
	{
		coord3 c = (*this); c.o -= v;
		return c;
	}
    DEVICE_CALLABLE void operator-=(const vec3& v)
	{
		*this = *this - v;
	}

	// 乘法：在坐标系下定义一个向量
    DEVICE_CALLABLE friend vec3 operator*(const vec3& p, const coord3& c)
	{
		return c.ux * (c.s.x * p.x) + c.uy * (c.s.y * p.y) + c.uz * (c.s.z * p.z) + c.o;
	}
    DEVICE_CALLABLE friend void operator*=(vec3& p, const coord3& c)
	{
		p = p * c;
	}
    DEVICE_CALLABLE coord3 operator*(const vec3& v) const
	{
		return (*this) * coord3(vec3::UX * v.x, vec3::UY * v.y, vec3::UZ * v.z);
	}
    DEVICE_CALLABLE void operator*=(const vec3& v)
	{
		*this = (*this) * v;
	}
    DEVICE_CALLABLE friend real operator*(const real& s, const coord3& c)
	{
		return s * ((c.s.x + c.s.y + c.s.z) / 3.0);
	}
    DEVICE_CALLABLE coord3 operator*(real s) const
	{
		coord3 c = *this;
		{// C*S 缩放乘法
			c.s.x *= s; c.s.y *= s; c.s.z *= s;
		}
		return c;
	}
    DEVICE_CALLABLE void operator*=(real s)
	{
		*this = (*this) * s;
	}
    DEVICE_CALLABLE coord3 operator*(const coord3& c) const
	{// Cchild * Cparent * ...
		coord3 rc = vcoord3::operator*(c);
		rc.o = c.o + (o.x * c.s.x) * c.ux + (o.y * c.s.y) * c.uy + (o.z * c.s.z) * c.uz;
		return rc;
	}
    DEVICE_CALLABLE void operator*=(const coord3& c)
	{
		*this = (*this) * c;
	}
    DEVICE_CALLABLE coord3 operator*(const quaternion& q) const
	{
		coord3 rc = *this;
		rc.ux = q * ux;
		rc.uy = q * uy;
		rc.uz = q * uz;
		rc.o = q * rc.o;
		return rc;
	}
    DEVICE_CALLABLE void operator*=(const quaternion& q)
	{
		*this = (*this) * q;
	}

	// 除法：向量向坐标系投影 注意：要保证ux,uy,uz是单位向量！
    DEVICE_CALLABLE friend vec3 operator/(const vec3& p, const coord3& c)
	{
		vec3 v = p - c.o;
		v = v / c.s;
		return vec3(v.dot(c.ux), v.dot(c.uy), v.dot(c.uz));
	}
    DEVICE_CALLABLE friend void operator/=(vec3& p, const coord3& c)
	{
		p = p / c;
	}
    DEVICE_CALLABLE coord3 operator/(const vec3& v) const
	{
		return (*this) / coord3(vec3::UX * v.x, vec3::UY * v.y, vec3::UZ * v.z);
	}
    DEVICE_CALLABLE void operator/=(const vec3& v)
	{
		*this = (*this) / v;
	}

	DEVICE_CALLABLE coord3 operator/(real s) const
	{// C/S 缩放除法
		coord3 c = *this;
		c.s /= s;
		c.o /= s;
		return c;
	}
    DEVICE_CALLABLE void operator/=(real s)
	{
		*this = (*this) / s;
	}
	// oper(/) = C1 * C2^ - 1
    DEVICE_CALLABLE coord3 operator/(const coord3& c) const
	{
		coord3 rc = vcoord3::operator/(c);
		rc.o = o - c.o;
		rc.o = vec3(rc.o.dot(c.ux) / c.s.x, rc.o.dot(c.uy) / c.s.y, rc.o.dot(c.uz) / c.s.z);
		return rc;
	}
    DEVICE_CALLABLE void operator/=(const coord3& c)
	{
		*this = (*this) / c;
	}
    DEVICE_CALLABLE coord3 operator/(const quaternion& q) const
	{
		return (*this) * q.conjcopy();
	}
    DEVICE_CALLABLE void operator/=(const quaternion& q)
	{
		*this = (*this) / q;
	}
	// oper(\) = C1^-1 * C2
    DEVICE_CALLABLE coord3 operator%(const coord3& c) const
	{
		return (*this).reversed() * c;
	}
    DEVICE_CALLABLE coord3 operator^(const vec3& v) const
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
	// 倒置
	void reverse()
	{
		(*this) = ONE / (*this);
	}
	coord3 reversed() const
	{
		return ONE / (*this);
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
	std::string serialise() const
	{
		vec3 eu = coord2eulers();
		return o.serialise() + "," + eu.serialise();
	}
	void dump(const std::string& name = "") const
	{
		PRINT("|/_ : " << name);
		PRINTVEC3(ux);
		PRINTVEC3(uy);
		PRINTVEC3(uz);
		PRINTVEC3(s);
		PRINTVEC3(o);
		PRINT("");
	}
	// GUID
	std::size_t hash() const
	{
		std::size_t hash = 0;
		if (sizeof(real) == sizeof(float))
		{
			// 哈希原点坐标
			hash_combine(hash, o.x);
			hash_combine(hash, o.y);
			hash_combine(hash, o.z);

			// 哈希缩放因子
			hash_combine(hash, s.x);
			hash_combine(hash, s.y);
			hash_combine(hash, s.z);

			// 哈希基向量
			hash_combine(hash, ux.x);
			hash_combine(hash, ux.y);
			hash_combine(hash, ux.z);

			hash_combine(hash, uy.x);
			hash_combine(hash, uy.y);
			hash_combine(hash, uy.z);

			hash_combine(hash, uz.x);
			hash_combine(hash, uz.y);
			hash_combine(hash, uz.z);
		}
		return hash;
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
	coord3 roted(real ang, const vec3& ax) const
	{
		coord3 c = *this;

		quaternion q(ang, ax);
		c.ux = q * ux;
		c.uy = q * uy;
		c.uz = q * uz;

		/*c.ux.rot(ang, ax);
		c.uy.rot(ang, ax);
		c.uz.rot(ang, ax);*/
		return c;
	}
	coord3 roted(const quaternion& q) const
	{
		coord3 c = *this;
		c.ux *= q;
		c.uy *= q;
		c.uz *= q;
		return c;
	}
	void moveto(const coord3& c, real alpha = 1)
	{
		o = vec3::lerp(o, c.o, alpha);
	}
};
#if defined(PMDLL) || !defined(PM_IMPLEMENTED)
const coord3 coord3::ZERO = {ucoord3::ONE, vec3::ZERO, vec3::ZERO };
const coord3 coord3::ONE = {};
#endif
