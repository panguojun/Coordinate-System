/*********************************************************************
*							����ϵ
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
* 	����ϵ�����ҵ�����װ�����ڼ�����任������������㷨���ܽ��һЩ
* 	����ϵ�任��ص����⡣
* 	����ϵ���������Ⱥ�����ơ�
*	����ϵ������������ɣ�C = M(λ�ã� + S�����ţ� * R����ת��
*  *  *  *  *  *  *  *  *  *  ���  *  *  *  *  *  *  *  *  *  *  *  *
*	����ϵ������� C������ϵ֮��ı任����д��G = C1//C2,GRAD�ݶȵ���˼
*			oper(/)  = C1 * C2^-1
*			oper(//) = C1^-1 * C2, oper(//) = gradcoord()
*	����ϵ��������: [C1,C2] = C1*C2 - C2*C1
*	������˵��
*	����һ����������ϵ(��������ƽֱ�ռ䣬�������������ƶ�������)��V,����
*	������ϵ�¹۲�V����ͬ����V�ǲ�ͬ�ģ��ʶ�����ϵ��λ���йأ�ȡ��������
*	��1),(2)�㴦������V1,V2����Ӧ����ϵC1,C2����ô��
*			V = V1 * C1 = V2 * C2 =>
*			V2 = V1 * C1 / C2, �� G12 = C1 / C2 =>
*			V2 = V1 * G12
*	
*	����������ϵ����������ϵx,y���ƽ����ͶӰ�õ���u,v������G12�ֱ�������
*	�����϶�ӦGu,Gv, ��(u1,v1)��(u2,v2) ��������·���Ĳ���ټ����������
*	�����ʹ�ʽΪ��
*			Ruv = Gu*Gv - Gv*Gu*Gu*Wu*Gv*Wv
*			W = (U + V*Gu) - (V + U*Gv)
*/

//#define	Parallel_Projection		 // ����������ϵ��ƽ��ͶӰ

// *******************************************************************
//  |_
// C     2d Coordinate System
// *******************************************************************
struct coord2
{
	static const coord2 ONE;

	vec2 ux = vec2::UX;		// ����
	vec2 uy = vec2::UY;

	vec2 s = vec2::ONE;		// ����
	vec2 o;				// ԭ��

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
	coord2(real ang)
	{
		ux.rot(ang);
		uy.rot(ang);
	}
	vec2 VX() const { return ux * s.x; }
	vec2 VY() const { return uy * s.y; }

	void rot(real ang)
	{
		ux.rot(ang);
		uy.rot(ang);
	}
	bool is_same_dirs(const coord2& c) const
	{
		return ux == c.ux && uy == c.uy;
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
		coord2 rc = (*this);
		rc.o = o + v;
		return rc;
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
		coord2 rc = (*this);
		rc.o = o - v;
		return rc;
	}
	void operator -= (const vec2& v)
	{
		*this = *this - v;
	}
	// ������ϵ�¶���һ������
	friend vec2 operator * (crvec2 p, const coord2& c)
	{
		return c.ux * (c.s.x * p.x) + c.uy * (c.s.y * p.y) + c.o;
	}
	coord2 operator * (crvec2 p) const
	{
		coord2 c = *this;
		c.ux = lerp(vec2::UX, c.VX(), p.x);
		c.uy = lerp(vec2::UY, c.VY(), p.y);
		c.norm();
		c.o.x *= p.x; c.o.y *= p.y;
		return c;
	}
	coord2 operator * (const coord2& c) const
	{
		coord2 rc;
		rc.ux = ux.x * c.ux + ux.y * c.uy;;
		rc.uy = uy.x * c.ux + uy.y * c.uy;;
		rc.s = s * c.s;
		rc.o = o + c.ux * o.x + c.uy * o.y;
		return rc;
	}
#ifdef Parallel_Projection
	// ����������ϵ��ƽ��ͶӰ Parallel projection
	static real pl_dot(crvec2 v, crvec2 ax1, crvec2 ax2)
	{
		real co = ax1.dot(ax2);
		real si = sqrt(1 - co * co);
		real sc = (co / si);
		return (v.dot(ax1) - v.cross(ax1) * sc);
	}
#endif
	// ����������ϵͶӰ ע�⣺Ҫ��֤ux,uy�ǵ�λ������
	friend vec2 operator / (crvec2 p, const coord2& c)
	{
		vec2 v = p - c.o;
#ifdef Parallel_Projection
		{// ���ڷ��������
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
		{// ���ڷ��������
			rc.ux = vec2(pl_dot(ux, c.ux, c.uy) / c.s.x, pl_dot(ux, c.uy, c.ux) / c.s.y);
			rc.uy = vec2(pl_dot(uy, c.ux, c.uy) / c.s.x, pl_dot(uy, c.uy, c.ux) / c.s.y);
	}
#else
		rc.ux = vec2(ux.dot(c.ux) / c.s.x, ux.dot(c.uy) / c.s.y);
		rc.uy = vec2(uy.dot(c.ux) / c.s.x, uy.dot(c.uy) / c.s.y);
#endif
		rc.o -= c.o;
		return rc;
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
	void reverse()
	{
		(*this) = ONE / (*this);
	}
	coord2 reversed() const
	{
		return ONE / (*this);
	}
	// �ݶ�����ϵ
	static coord2 gradcoord(const coord2& c1, const coord2& c2)
	{
		return c1.reversed() * c2;
	}
	real dot(crvec2 v) const
	{
		return v.dot(ux) * s.x + v.dot(uy) * s.y;
	}
	void dump() const
	{
		//PRINT("-------");
		PRINT("ux: " << ux.x << "," << ux.y);
		PRINT("uy: " << uy.x << "," << uy.y);
		PRINTVEC2(s);
		//PRINT("uz: " << uz.x << "," << uz.y << "," << uz.z);
		//PRINT("o: " << o.x << "," << o.y << "," << o.z);
	}
};
const coord2 coord2::ONE = coord2();

// ******************************************************************
//  |/_
// C     3d Coordinate System
// ******************************************************************
struct coord3
{
	static const coord3 ZERO;
	static const coord3 ONE;

	vec3 ux = vec3::UX;		// ����
	vec3 uy = vec3::UY;
	vec3 uz = vec3::UZ;

	vec3 s = vec3::ONE;		// ����
	vec3 o;					// ԭ��

	coord3() {}
	coord3(const coord3& c)
	{
		ux = c.ux; uy = c.uy; uz = c.uz;
		s = c.s;
		o = c.o;
	}
	coord3(crvec _ux, crvec _uy, crvec _uz)
	{
		ux = _ux; uy = _uy; uz = _uz;
	}
	coord3(crvec _ux, crvec _uy)
	{
		ux = _ux; uy = _uy; uz = ux.cross(uy);
	}
	explicit coord3(crvec _p)
	{
		o = _p;
	}
	coord3(real ang, crvec ax)
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
	// �ƶ���
	void fromvectorsT(crvec v1, crvec v2)
	{
		quaternion q;
		q.fromvectors(v1, v2);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
	}
	// ��ת��
	void fromvectorsR(crvec v1, crvec v2)
	{
		o = v2 - v1;
	}
	void fromaxvecs(crvec ax, crvec v1, crvec v2)
	{
		vec3 pv1 = v1.crossdot(ax);
		vec3 pv2 = v2.crossdot(ax);
		real ang = acos(pv1.dot(pv2));
		quaternion q;
		q.fromangleaxis(ang, ax);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;;
	}

	vec3 VX() const { return ux * s.x; }
	vec3 VY() const { return uy * s.y; }
	vec3 VZ() const { return uz * s.z; }

	vec3 X() const { return ux * s.x + vec3::UX * o.x; }
	vec3 Y() const { return uy * s.y + vec3::UX * o.y; }
	vec3 Z() const { return uz * s.z + vec3::UX * o.z; }

	// ��һ������������ϵ
	coord3 ucoord() const
	{
		coord3 c = *this;
		c.norm(false);
		c.o = vec3::ZERO;
		return c;
	}
	// λ��
	inline vec3 pos()
	{
		return o;
	}
	// ���� X ����
	inline coord3 vcoord()
	{
		coord3 c = *this;
		c.o = vec3::ZERO;
		return c;
	}
	quaternion toquat() const
	{
		coord3 c = ucoord();
		vec3 pyr = c.coord2eulers();
		quaternion q;
		q.fromeuler(pyr.x, pyr.y, pyr.z);
		return q;
	}
	bool is_same_dirs(const coord3& c) const
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
	
	// +/- ����
	coord3 operator + (const coord3& c) const
	{
		coord3 rc;
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
	{// C+V �ƶ�
		coord3 rc = (*this);
		rc.o = o + v;
		return rc;
	}
	void operator += (const vec3& v)
	{
		*this = *this + v;
	}
	coord3 operator - (const coord3& c) const
	{
		coord3 rc;
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

	// �˷���������ϵ�¶���һ������
	friend vec3 operator * (crvec p, const coord3& c)
	{
		return c.ux * (c.s.x * p.x) + c.uy * (c.s.y * p.y) + c.uz * (c.s.z * p.z) + c.o;
	}
	friend void operator *= (rvec p, const coord3& c)
	{
		p = p * c;
	}
	coord3 operator * (crvec p) const
	{// C*V ���ų˷�
		coord3 c = *this;
		c.s.x *= p.x;c.s.y *= p.y;c.s.z *= p.z;
		return c;
	}
	coord3 operator * (const coord3& c) const
	{// Cchild * Cparent * ...
		coord3 rc;
		rc.ux = ux.x * c.ux + ux.y * c.uy + ux.z * c.uz;
		rc.uy = uy.x * c.ux + uy.y * c.uy + uy.z * c.uz;
		rc.uz = uz.x * c.ux + uz.y * c.uy + uz.z * c.uz;

		rc.s = s * c.s;
		rc.o = c.o + c.ux * o.x + c.uy * o.y + c.uz * o.z;
		return rc;
	}
	void operator *= (const coord3& c)
	{
		*this = (*this) * c;
	}
	coord3 operator * (const quaternion& q) const
	{
		coord3 rc = *this;
		rc.ux = q * ux;
		rc.uy = q * uy;
		rc.uz = q * uz;
		return rc;
	}
	void operator *= (const quaternion& q)
	{
		*this = (*this) * q;
	}

	// ����������������ϵͶӰ ע�⣺Ҫ��֤ux,uy,uz�ǵ�λ������
#ifdef Parallel_Projection
	// ����������ϵ��ƽ��ͶӰ Parallel projection
	static real pl_prj(crvec v, crvec ax1, crvec ax2)
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
	friend vec3 operator / (crvec p, const coord3& c)
	{
		vec3 v = p - c.o;
#ifdef Parallel_Projection
		{// ���ڷ��������
			return vec3(
				pl_prj(v - c.uz * v.dot(c.uz), c.ux, c.uy) / c.s.x,
				pl_prj(v - c.ux * v.dot(c.ux), c.uy, c.uz) / c.s.y,
				pl_prj(v - c.uy * v.dot(c.uy), c.uz, c.ux) / c.s.z);
		}
#endif
		return vec3(v.dot(c.ux) / c.s.x, v.dot(c.uy) / c.s.y, v.dot(c.uz) / c.s.z);
	}
	friend void operator /= (vec p, const coord3& c)
	{
		p = p / c;
	}
	// oper(/) = C1 * C2^-1
	coord3 operator / (const coord3& c) const
	{
		coord3 rc;
#ifdef Parallel_Projection
		{// ���ڷ��������
			rc.ux = PL_PRJ3(ux);
			rc.uy = PL_PRJ3(uy);
			rc.uz = PL_PRJ3(uz);
		}
#else
		rc.ux = vec3(ux.dot(c.ux) / c.s.x, ux.dot(c.uy) / c.s.y, ux.dot(c.uz) / c.s.z);
		rc.uy = vec3(uy.dot(c.ux) / c.s.x, uy.dot(c.uy) / c.s.y, uy.dot(c.uz) / c.s.z);
		rc.uz = vec3(uz.dot(c.ux) / c.s.x, uz.dot(c.uy) / c.s.y, uz.dot(c.uz) / c.s.z);
#endif
		rc.o -= c.o;
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
	// oper(^)
	// ��ռ�ĳ˷�����,Ce^(th*v)
	// ��C��ʾĳ����A����������ת��
	// �ں�����0<v<1,c=C^v; v=0ʱc=ONE,v=1ʱc=C
	coord3 operator ^ (crvec v) const
	{
		coord3 c = *this;
		c.ux = lerp(vec3::UX, c.VX(), v.x);
		c.uy = lerp(vec3::UY, c.VY(), v.y);
		c.uz = lerp(vec3::UZ, c.VZ(), v.z);
		c.norm();
		return c;
	}
	
	// ��һ��
	void norm(bool bscl = true)
	{
#define ISZERO(a) (fabs(a) < 1e-10)
		s.x = ux.len(); if (!ISZERO(s.x)) ux /= s.x;
		s.y = uy.len(); if (!ISZERO(s.y)) uy /= s.y;
		s.z = uz.len(); if (!ISZERO(s.z)) uz /= s.z;
		if (!bscl)
			s = vec3::ONE;
	}
	// ת��
	void transpose()
	{
		vec3 ux = vec3(ux.x, uy.x, uz.x);
		vec3 uy = vec3(ux.y, uy.y, uz.y);
		vec3 uz = vec3(ux.z, uy.z, uz.z);
		(*this).ux = ux;
		(*this).uy = uy;
		(*this).uz = uz;
	}
	coord3 transposed()
	{
		coord3 c = (*this);
		c.ux = vec3(ux.x, uy.x, uz.x);
		c.uy = vec3(ux.y, uy.y, uz.y);
		c.uz = vec3(ux.z, uy.z, uz.z);
		return c;
	}
	// ����
	void reverse()
	{
		(*this) = ONE / (*this);
	}
	coord3 reversed() const
	{
		return ONE / (*this);
	}
	// ��ת
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
	void rot(real ang, crvec ax)
	{
		ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);
	}
	vec3 sumvec() const
	{
		return ux * s.x + uy * s.y + uz * s.z;
	}
	// ��������������ϵ��Ϊ��ת�任ʱ���������
	vec3 eigenvec() const
	{
		return toquat().axis();
	}
	real dot(crvec v) const
	{
		return v.dot(ux) * s.x + v.dot(uy) * s.y + v.dot(uz) * s.z;
	}
	real dot(const coord3& c) const
	{
		return c.VX().dot(VX()) + c.VY().dot(VY()) + c.VZ().dot(VZ());
	}
	// ������������Ĳ�ˣ����ӷ���Ⱥ��
	coord3 lie_cross(const coord3& c) const
	{
		return (*this) * c - c * (*this);
	}
	// �ɵ�ų����������Ĳ��
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
	// ����ϵ��ŷ���ǣ�Ҫ��֤�ǹ�һ������������ϵ
	vec3 coord2eulers() const
	{
		const coord3& rm = *this;
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
	// �ݶ�����ϵ = �ݶ� X �пռ�
	// �൱��һ������ϵ�ĵ���
	static coord3 gradcoord(const coord3& c1, const coord3& c2)
	{
		return c1.reversed() * c2;
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