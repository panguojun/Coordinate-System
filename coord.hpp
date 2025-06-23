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
*   Take vectors V1 and V2 at adjacent points (1) and (2) respectively,
*   corresponding to coordinate systems C1 and C2. Then:
*           		V  = V1 / C1 = V2 / C2 =>
*           		V2 = V1 * C2 / C1, let R12 = C2 / C1 =>
*           		V2 = V1 * R12
*
*   The coordinate system can be used to calculate spatial curvature. In the u,v coordinate system,
*   the Curvature tensor is:
*           	Ruv  = 	Gu*Gv - Gv*Gu - G[u,v]
*           	where:  Gu = C2 / C1 - I
*                   	Connection vector: W = [U, V] (Lie bracket operation)
*                   	G[u,v] = Gu*Wu + Gv*Wv
*/

//#define	NON_UNIFORM_SCALE	// For Differential Geometry
// ********************************************************************************************
//  |/_
// UC     3d Rotation Coordinate System（Base Coordinate System）
// ********************************************************************************************
struct ucoord3
{
	static const ucoord3 ONE;
	union {
		struct {
			// basis
			vec3 ux;
			vec3 uy;
			vec3 uz;
		};
		real m[9]; // matrix
	};
	ucoord3()
	{
		ux = vec3::UX; uy = vec3::UY; uz = vec3::UZ;
	}
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
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}
    	ucoord3(real pit, real yaw, real rol)
	{
		quaternion q(pit, yaw, rol);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
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
	vec3 toeulers() const
	{
		return coord2eulers();
	}
	void from_vecs_R(const vec3& v1, const vec3& v2)
	{
		vec3  v = v1.cross(v2);
		real c = v1.dot(v2);
		real k = 1.0 / (1.0 + c);
		
		ux = { v.x * v.x * k + c,      v.y * v.x * k - v.z,    v.z * v.x * k + v.y };
		uy = { v.x * v.y * k + v.z,    v.y * v.y * k + c,      v.z * v.y * k - v.x };
		uz = { v.x * v.z * k - v.y,    v.y * v.z * k + v.x,    v.z * v.z * k + c   };
	}
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

	// Multiplication: Define a vector in the coordinate system or restore the vector to the parent space
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
	// Division: Project the vector onto the coordinate system (for non - orthogonal coordinate systems, further expansion is recommended).
    	friend vec3 operator/(const vec3& v, const ucoord3& c)
	{
		return vec3(v.dot(c.ux), v.dot(c.uy), v.dot(c.uz));
	}
    	friend void operator/=(vec3& v, const ucoord3& c)
	{
		v = v / c;
	}
	// oper(/) = C1 * C2^-1
    	ucoord3 operator/(const ucoord3& c) const
	{
		ucoord3 rc;
		rc.ux = vec3(ux.dot(c.ux), ux.dot(c.uy), ux.dot(c.uz));
		rc.uy = vec3(uy.dot(c.ux), uy.dot(c.uy), uy.dot(c.uz));
		rc.uz = vec3(uz.dot(c.ux), uz.dot(c.uy), uz.dot(c.uz));
		return rc;
	}
    	void operator/=(const ucoord3& c)
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
	// Operation (^)
	// Multiplication operation in phase space, C^v = C·e^(th·v)
	// For example, if C represents the rotation of a vector A between two points,
	// when fusing with a scalar 0 < v < 1, the result is c = C^v; when v=0, c=ONE; when v=1, c=C
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
		return ucoord3((*this).toquat() ^ f);
	}
	void operator ^= (real f)
	{
		(*this) = (*this) ^ f;
	}
	// Transpose (swap coordinate axes)
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
	// Inversion
	void reverse()
	{
		(*this) = ONE / (*this);
	}
	ucoord3 reversed() const
	{
		return ONE / (*this);
	}
	// Flip
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
	void rot(real angle, const vec3& ax)
	{
		quaternion q(angle, ax);
		ux = q * ux;
		uy = q * uy;
		uz = q * uz;
	}
	void rot(const quaternion& q)
	{
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
	// Cross product derived from electromagnetic field calculations
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
	// Coordinate system to Euler angles
	vec3 coord2eulers() const
	{
		real c_eps = 1e-5;

		const ucoord3& rm = *this;
		real sy = sqrt(rm.ux.x * rm.ux.x + rm.uy.x * rm.uy.x);
		bool singular = sy < c_eps;

		real x, y, z;
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
		real x = eulers.x;
		real y = eulers.y;
		real z = eulers.z;

		real cx = cos(x);
		real sx = sin(x);
		real cy = cos(y);
		real sy = sin(y);
		real cz = cos(z);
		real sz = sin(z);

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
	// Convert to quaternion (note quaternion multiplication order)
	quaternion toquat() const
	{
		return quat(coord2eulers());
	}
	quaternion Q() const
	{
		return quat(coord2eulers());
	}
	void Q(const quaternion& q)
	{
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}
	void Q(real qw, real qx, real qy, real qz)
	{
		quaternion q(qw, qx, qy, qz);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}

	// Gradient coordinate system = Gradient X Tangent space
	// Equivalent to the derivative of the first-order coordinate system
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

	// Convenience function (note: angle(u, _u) != +/-PI)
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
#if !defined(PM_IMPLEMENTED)
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
	vcoord3(real x, real y, real z) : s(x, y, z) {}
	vcoord3(const ucoord3& c) : ucoord3(c){}
	vcoord3(const ucoord3& c, const vec3& _s) : ucoord3(c), s(_s){}
	vcoord3( const vec3& _ux, const vec3& _uy, const vec3& _uz, const vec3& _s) : ucoord3(_ux, _uy, _uz), s(_s){ }
	vcoord3( const vec3& _ux, const vec3& _uy, const vec3& _uz) : ucoord3(_ux, _uy, _uz){}
	vcoord3(const quaternion& q) : ucoord3(q){}
	vcoord3(const quaternion& q, const vec3& _s) : ucoord3(q), s(_s) {}
	vcoord3(const vec3& _s) : s(_s) {}
	vcoord3(const real& _s) : s(_s, _s, _s) {}

	vec3 VX() const { return ux * s.x; }
	vec3 VY() const { return uy * s.y; }
	vec3 VZ() const { return uz * s.z; }

	const ucoord3& base() const
	{
		return static_cast<const ucoord3&>(*this);
	}
	void base(const ucoord3& ucd)
	{
		ux = ucd.ux; uy = ucd.uy; uz = ucd.uz;
	}
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
	vcoord3 operator * (real _s) const
	{
		vcoord3 c = *this;
		{// C*S 缩放乘法
			c.s *= _s;
		}
		return c;
	}
	void operator *= (real _s)
	{
		*this = (*this) * _s;
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
	
	// 除法：向量向坐标系投影（对于非正交坐标系，建议再扩展）
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

	vcoord3 operator / (real _s) const
	{// C/S 缩放除法
		vcoord3 c = *this;
		c.s /= _s;
		return c;
	}
	void operator /= (real _s)
	{
		*this = (*this) / _s;
	}
	// oper(/) = C1 * C2^-1
	vcoord3 operator / (const vcoord3& c) const
	{
		vcoord3 rc;
#ifdef	NON_UNIFORM_SCALE
		vec3 vx = VX();
		vec3 vy = VY();
		vec3 vz = VZ();

		vec3 cvx = c.ux.normcopy() / c.s.x;
		vec3 cvy = c.uy.normcopy() / c.s.y;
		vec3 cvz = c.uz.normcopy() / c.s.z;

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

	// Cross Product 由电磁场计算引出的叉乘
	vcoord3 cross(const vcoord3& c) const
	{
		vec3 vx = VX();
		vec3 vy = VY();
		vec3 vz = VZ();

		vec3 cvx = c.VX();
		vec3 cvy = c.VY();
		vec3 cvz = c.VZ();

		return vcoord3(
			vec3::UX * (vy.dot(cvz) - vz.dot(cvy)),
			vec3::UY * (vz.dot(cvx) - vx.dot(cvz)),
			vec3::UZ * (vx.dot(cvy) - vy.dot(cvx))
		);
	}
	// v1 x v2 = v1 * (C x v2)
	vcoord3 cross(const vec3& v) const
	{
		return vcoord3(
			VX().cross(v),
			VY().cross(v),
			VZ().cross(v)
		);
	}

	// Dot Product 
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
#if  !defined(PM_IMPLEMENTED)
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

	union {
		vec3 o = vec3::ZERO;			// 原点
		struct {
			real x, y, z;
		};
	};

	coord3() : o(0, 0, 0) {}
	coord3(const coord3& other) : o(other.o), vcoord3(other.ux, other.uy, other.uz, other.s) {}
	coord3(real x, real y, real z) : o(x, y, z) {}

	coord3(const ucoord3& uc) : vcoord3(uc){}
	coord3(const vcoord3& vc) : vcoord3(vc) {}

	coord3(const vec3& _o) : o(_o) {}
	coord3(const vec3& _o,  const vec3& _s,  const vec3& _ux, const vec3& _uy, const vec3& _uz) : vcoord3(_ux, _uy, _uz, _s), o(_o){}
	coord3(const vec3& _o,  const vec3& _ux, const vec3& _uy, const vec3& _uz) : vcoord3(_ux, _uy, _uz), o(_o){ }
	coord3(const vec3& _o,  const ucoord3& c) : vcoord3(c), o(_o){}
	coord3(const vec3& _o,  const vec3& _s,  const ucoord3& c) : vcoord3(c, _s), o(_o) {}

	coord3(const vec3& _ux, const vec3& _uy, const vec3& _uz) : vcoord3(_ux, _uy, _uz) {}
	coord3(const vec3& _ux, const vec3& _uy) : vcoord3(_ux, _uy, ux.cross(uy)) {}

	coord3(const ucoord3& c,const vec3& _s, const vec3& _o) : vcoord3(c, _s), o(_o){}
	coord3(const ucoord3& c,const vec3& _o) : vcoord3(c), o(_o) {}

	coord3(real ang, const vec3& ax) : vcoord3(quaternion(ang, ax)) {}
    	coord3(const quaternion& q) : vcoord3(q) {}
    	coord3(const vec3& p, const quaternion& q, const vec3& _s = vec3::ONE) : vcoord3(q, _s), o(p) {}
	coord3(real x, real y, real z, real qw, real qx, real qy, real qz, real sx, real sy, real sz) : vcoord3(quaternion(qw, qx, qy, qz), vec3(sx, sy, sz)), o(x, y, z) {}
	coord3(real x, real y, real z, real qw, real qx, real qy, real qz) : vcoord3(quaternion(qw, qx, qy, qz), s), o(x, y, z) {}
	coord3(real x, real y, real z, real rx, real ry, real rz)
	{
		real ang2rad = PI / 180.0;
		quaternion q(rx * ang2rad, ry * ang2rad, rz * ang2rad);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
		o = vec3(x, y, z);
	}

	static coord3 from_axes(const vec3& ux, const vec3& uy, const vec3& uz) {
		return coord3(vec3::ZERO, vec3::ONE, ux, uy, uz);
	}
	static coord3 from_angle(real angle, const vec3& axis) {
		quaternion q(angle, axis);
		return coord3(vec3::ZERO, q);
	}

	operator quaternion() const
	{
		return toquat();
	}
    	operator vec3() const
	{
		return o;
	}

    	vec3 VX()	const { return ux * s.x; }
	vec3 VY()	const { return uy * s.y; }
	vec3 VZ()	const { return uz * s.z; }

	void VX(const vec3& vx)	{ real r = vx.len(); ux = vx / r; s.x = r; }
	void VY(const vec3& vy)	{ real r = vy.len(); uy = vy / r; s.y = r; }
	void VZ(const vec3& vz)	{ real r = vz.len(); uz = vz / r; s.z = r; }

	vec3 X()	const { return ux * s.x + vec3::UX * o.x; }
	vec3 Y()	const { return uy * s.y + vec3::UY * o.y; }
	vec3 Z()	const { return uz * s.z + vec3::UZ * o.z; }

	// 旋转坐标系
    	const ucoord3& ucoord() const
	{
		return static_cast<const ucoord3&>(*this);
	}
    	void ucoord(const ucoord3& ucd)
	{
		ux = ucd.ux; uy = ucd.uy; uz = ucd.uz;
	}
    	void ucoord(vec3 _ux, vec3 _uy, vec3 _uz)
	{
		ux = _ux; uy = _uy; uz = _uz;
	}
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
    	void UC(vec3 _ux, vec3 _uy, vec3 _uz)
	{
		ux = _ux; uy = _uy; uz = _uz;
	}
	// 向量坐标系 = 方向 X 缩放
    	const vcoord3& vcoord() const
	{
		return static_cast<const vcoord3&>(*this);
	}
    	const vcoord3& VC() const
	{
		return static_cast<const vcoord3&>(*this);
	}
	// 姿态
    	coord3 pose()
	{
		return { ucoord(), vec3::ONE, o };
	}
	// 位置
    	vec3 pos() const
	{
		return o;
	}
	// 向量
    	vec3 tovec() const
	{
		return ux * s.x + uy * s.y + uz * s.z;
	}
    	coord3 operator=(const coord3& c)
	{
		o = c.o;
		s = c.s;
		ux = c.ux; uy = c.uy; uz = c.uz;
		return (*this);
	}
    	bool equal_dirs(const coord3& c) const
	{
		return ux == c.ux && uy == c.uy && uz == c.uz;
	}
    	bool operator==(const coord3& c) const
	{
		return o == c.o && s == c.s && equal_dirs(c);
	}
    	bool operator!=(const coord3& c) const
	{
		return o != c.o || s != c.s || !equal_dirs(c);
	}
	// +/- 运算
    	coord3 operator+(const coord3& c) const
	{
		coord3 rc;
		vec3 _ux = VX() + c.VX();
		vec3 _uy = VY() + c.VY();
		vec3 _uz = VZ() + c.VZ();
		/*
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
		}*/

		rc.o = o + c.o;
		return rc;
	}
    	coord3 operator+=(const coord3& c)
	{
		*this = (*this) + c;
		return *this;
	}
    	coord3 operator+(const vec3& v) const
	{
		coord3 c = (*this); c.o += v;
		return c;
	}
    	coord3 operator+=(const vec3& v)
	{
		*this = *this + v;
		return *this;
	}
    	friend vec3 operator+(const vec3& p, const coord3& c)
	{
		return p + c.o;
	}
    	friend void operator+=(vec3& p, const coord3& c)
	{
		p = p + c;
	}
    	friend vec3 operator-(const vec3& p, const coord3& c)
	{
		return p - c.o;
	}
    	friend void operator-=(vec3& p, const coord3& c)
	{
		p = p - c;
	}
    	coord3 operator-(const coord3& c) const
	{
		coord3 rc;
		vec3 _ux = VX() - c.VX();
		vec3 _uy = VY() - c.VY();
		vec3 _uz = VZ() - c.VZ();
		/*
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
		}*/

		rc.o = o - c.o;
		return rc;
	}
	coord3 operator-() const
	{
		coord3 c = (*this);
		c.o = -c.o;
		return c;
	}
    	coord3 operator-(const vec3& v) const
	{
		coord3 c = (*this); c.o -= v;
		return c;
	}
    	coord3 operator-=(const vec3& v)
	{
		*this = *this - v;
		return *this;
	}

	// 乘法：在坐标系下定义一个向量
    	friend vec3 operator*(const vec3& p, const coord3& c)
	{
		return c.ux * (c.s.x * p.x) + c.uy * (c.s.y * p.y) + c.uz * (c.s.z * p.z) + c.o;
	}
    	friend void operator*=(vec3& p, const coord3& c)
	{
		p = p * c;
	}
    	coord3 operator*(const vec3& v) const
	{
		return (*this) * coord3(vec3::UX * v.x, vec3::UY * v.y, vec3::UZ * v.z);
	}
    	void operator*=(const vec3& v)
	{
		*this = (*this) * v;
	}
    	coord3 operator*(real _s) const
	{
		coord3 c = *this;
		{// C*S 缩放乘法
			c.s *= _s;
			c.o *= _s;
		}
		return c;
	}
    	void operator*=(real _s)
	{
		*this = (*this) * _s;
	}
    	coord3 operator*(const coord3& c) const
	{// Cchild * Cparent * ...
		coord3 rc = vcoord3::operator*(c);
		rc.o = c.o + (o.x * c.s.x) * c.ux + (o.y * c.s.y) * c.uy + (o.z * c.s.z) * c.uz;
		return rc;
	}
    	coord3 operator*=(const coord3& c)
	{
		*this = (*this) * c;
		return *this;
	}
	coord3 operator*(const vcoord3& c) const
	{// Cchild * Cparent * ...
		coord3 rc = vcoord3::operator*(c);
		rc.o = (o.x * c.s.x) * c.ux + (o.y * c.s.y) * c.uy + (o.z * c.s.z) * c.uz;
		return rc;
	}
	coord3 operator*=(const vcoord3& c)
	{
		*this = (*this) * c;
		return *this;
	}
	coord3 operator*(const ucoord3& c) const
	{// Cchild * Cparent * ...
		coord3 rc = ucoord3::operator*(c);
		rc.o = (o.x) * c.ux + (o.y) * c.uy + (o.z) * c.uz;
		return rc;
	}
	coord3 operator*=(const ucoord3& c)
	{
		*this = (*this) * c;
		return *this;
	}
    	coord3 operator*(const quaternion& q) const
	{
		coord3 rc = *this;
		rc.ux = q * ux;
		rc.uy = q * uy;
		rc.uz = q * uz;
		rc.o = q * rc.o;
		return rc;
	}
    	coord3 operator*=(const quaternion& q)
	{
		*this = (*this) * q;
		return *this;
	}

	// 除法：向量向坐标系投影（对于非正交坐标系，建议再扩展）
    	friend vec3 operator/(const vec3& p, const coord3& c)
	{
		vec3 v = p - c.o;
		v = v / c.s;
		return vec3(v.dot(c.ux), v.dot(c.uy), v.dot(c.uz));
	}
    	friend void operator/=(vec3& p, const coord3& c)
	{
		p = p / c;
	}
    	coord3 operator/(const vec3& v) const
	{
		return (*this) / coord3(vec3::UX * v.x, vec3::UY * v.y, vec3::UZ * v.z);
	}
    	void operator/=(const vec3& v)
	{
		*this = (*this) / v;
	}

	coord3 operator/(real _s) const
	{// C/S 缩放除法
		coord3 c = *this;
		c.s /= _s;
		c.o /= _s;
		return c;
	}
    	void operator/=(real _s)
	{
		*this = (*this) / _s;
	}
	// oper(/) = C1 * C2^ - 1
    	coord3 operator/(const coord3& c) const
	{
		coord3 rc = vcoord3::operator/(c);
		rc.o = o - c.o;
		rc.o = vec3(rc.o.dot(c.ux) / c.s.x, rc.o.dot(c.uy) / c.s.y, rc.o.dot(c.uz) / c.s.z);
		return rc;
	}
    	coord3 operator/=(const coord3& c)
	{
		*this = (*this) / c;
		return *this;
	}
	coord3 operator/(const vcoord3& c) const
	{
		coord3 rc = vcoord3::operator/(c);
		rc.o = o;
		rc.o = vec3(rc.o.dot(c.ux) / c.s.x, rc.o.dot(c.uy) / c.s.y, rc.o.dot(c.uz) / c.s.z);
		return rc;
	}
	coord3 operator/=(const vcoord3& c)
	{
		*this = (*this) / c;
		return *this;
	}
	coord3 operator/(const ucoord3& c) const
	{
		coord3 rc = ucoord3::operator/(c);
		rc.o = o;
		rc.o = vec3(rc.o.dot(c.ux), rc.o.dot(c.uy), rc.o.dot(c.uz));
		return rc;
	}
	coord3 operator/=(const ucoord3& c)
	{
		*this = (*this) / c;
		return *this;
	}
    	coord3 operator/(const quaternion& q) const
	{
		return (*this) * q.conjcopy();
	}
    	void operator/=(const quaternion& q)
	{
		*this = (*this) / q;
	}
	// oper(\) = C1^-1 * C2
    	coord3 operator%(const coord3& c) const
	{
		return (*this).reversed() * c;
	}
    	coord3 operator^(const vec3& v) const
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

		vec3 _s = vec3::lerp(vec3::ONE, s, t);
		vec3 _o = vec3::lerp(vec3::ZERO, o, t);

		return coord3(uc, _s, _o);
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
	// 梯度坐标系
	// V2 = V1 * (C2 / C1 - I)
	// G = C2 / C1 - I
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
	}
};
#if !defined(PM_IMPLEMENTED)
const coord3 coord3::ZERO = {ucoord3::ONE, vec3::ZERO, vec3::ZERO };
const coord3 coord3::ONE = {};
#endif
