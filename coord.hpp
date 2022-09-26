/********************************************************
*			坐标系
* *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
* 坐标系类是我单独封装，用于简化坐标变换，衍生出许多算法，能解决一些
* 坐标系变换相关的问题。
* 
* 坐标系是物理规范的简化，物理学因为规范而变得数学上很困难，如果能用
* 算法简化规范运算,可以还原一个初中生眼中的物理世界！
* 物理学范式: 测量 = 规范 * 本征 * 量纲
* 补充了2D 非正交情况下的平行投影
* 
*/

const real deta_d = 0.0001f;	// 空间微分精度
const real deta_t = 0.0001f;	// 时间微分精度
// ***********************************************
// coord2
// ***********************************************
struct coord2
{
	vec2 ux = vec2::UX;		// 方向
	vec2 uy = vec2::UY;

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
	coord2(real ang)
	{
		ux.rot(ang);
		uy.rot(ang);
	}
	vec2 VX() const { return ux * s.x; }
	vec2 VY() const { return uy * s.y; }

	void rot(real ang)
	{
		//	o.rot(ang, ax);
		ux.rot(ang);
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
	coord2 operator - (const coord2& c) const
	{
		coord2 rc;
		{
			rc.ux = VX() - c.VX();
			rc.uy = VY() - c.VY();
			//rc.norm();
		}
		rc.o = o - c.o;
		return rc;
	}
	// 在坐标系下定义一个向量
	vec2 operator * (crvec2 p) const
	{
		return ux * (s.x * p.x) + uy * (s.y * p.y) + o;
	}
	friend vec2 operator * (crvec2 p, const coord2& c)
	{
		return c.ux * (c.s.x * p.x) + c.uy * (c.s.y * p.y) + c.o;
	}
	coord2 operator * (const coord2& c) const
	{
		coord2 rc;
		rc.ux = ux * c.ux.x + uy * c.ux.y;
		rc.uy = ux * c.uy.x + uy * c.uy.y;
		rc.s = s * c.s;
		rc.o = o + ux * c.o.x + uy * c.o.y;
		return rc;
	}
	// Parallel projection
	static real pl_prj(crvec2 v, crvec2 ax1, crvec2 ax2)
	{
		real co = ax1.dot(ax2);
		real si = sqrt(1 - co * co);
		real sc = (co / si);
		return (v.dot(ax1)- v.cross(ax1) * sc);
	}
	// 向量向坐标系投影 注意：要保证ux,uy是单位向量！
	friend vec2 operator / (crvec2 p, const coord2& c)
	{
		vec2 v = p - c.o;
		//{// 对于非正交情况
		//	return vec2(pl_prj(v, c.ux, c.uy) / c.s.x, pl_prj(v, c.uy, c.ux) / c.s.y);
		//}
		return vec2(v.dot(c.ux) / c.s.x, v.dot(c.uy) / c.s.y);
	}
	coord2 operator / (const coord2& c) const
	{
		coord2 rc;
		//{// 对于非正交情况
		//	rc.ux = vec2(pl_prj(ux, c.ux, c.uy) / c.s.x, pl_prj(ux, c.uy, c.ux) / c.s.y);
		//	rc.uy = vec2(pl_prj(uy, c.ux, c.uy) / c.s.x, pl_prj(uy, c.uy, c.ux) / c.s.y);
		//}
		rc.ux = vec2(ux.dot(c.ux) / c.s.x, ux.dot(c.uy) / c.s.y);
		rc.uy = vec2(uy.dot(c.ux) / c.s.x, uy.dot(c.uy) / c.s.y);
		rc.o -= c.o;
		return rc;
	}
	void norm(bool bscl = true)
	{
#define ISZERO(a) (fabs(a) < 1e-10)
		s.x = ux.len(); if (!ISZERO(s.x)) ux /= s.x;
		s.y = uy.len(); if (!ISZERO(s.y)) uy /= s.y;

		if (!bscl)
			s = vec2::ONE;
	}
	// 本征向量
	vec2 eigenvec() const
	{
		return (ux + uy) * s;
	}
	real dot(crvec2 v) const
	{
		return v.dot(eigenvec());
	}
	void dump() const
	{
		//PRINT("-------");
		PRINT("ux: " << ux.x << "," << ux.y);
		PRINT("uy: " << uy.x << "," << uy.y);
		PRINTVEC2(s);
	}
};

// ***********************************************
// coord3
// ***********************************************
struct coord3
{
	vec3 ux = vec3::UX;		// 方向
	vec3 uy = vec3::UY;
	vec3 uz = vec3::UZ;

	vec3 s = vec3::ONE;		// 缩放
	vec3 o;				// 原点

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
		uz = ux.cross(uy);
	}
	coord3(real ang, crvec ax)
	{
		ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);
	}
	coord3(real pit, real yaw, real rol)
	{
		coord3 cx(pit, vec3::UX);
		coord3 cy(yaw, vec3::UY);
		coord3 cz(rol, vec3::UZ);
		*this = cx * cy * cz;
	}
	vec3 VX() const { return ux * s.x; }
	vec3 VY() const { return uy * s.y; }
	vec3 VZ() const { return uz * s.z; }

	void rot(real ang, crvec ax)
	{
		//o.rot(ang, ax);
		ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);
	}
	bool is_same_dirs(const coord3& c) const
	{
		return ux == c.ux && uy == c.uy && uz == c.uz;
	}
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
	coord3 operator - (const coord3& c) const
	{
		coord3 rc;
		{
			rc.ux = VX() - c.VX();
			rc.uy = VY() - c.VY();
			rc.uz = VZ() - c.VZ();
			//rc.norm();
		}
		rc.o = o - c.o;
		return rc;
	}
	// 在坐标系下定义一个向量
	vec3 operator * (crvec p) const
	{
		return ux * (s.x * p.x) + uy * (s.y * p.y) + uz * (s.z * p.z) + o;
	}
	friend vec3 operator * (crvec p, const coord3& c)
	{
		return c.ux * (c.s.x * p.x) + c.uy * (c.s.y * p.y) + c.uz * (c.s.z * p.z) + c.o;
	}
	coord3 operator * (const coord3& c) const
	{
		coord3 rc;
		rc.ux = ux.x * c.ux + uy.x * c.uy + uz.x * c.uz;
		rc.uy = ux.y * c.ux + uy.y * c.uy + uz.y * c.uz;
		rc.uz = ux.z * c.ux + uy.z * c.uy + uz.z * c.uz;
		rc.s = s * c.s;
		rc.o = o + ux * c.o.x + uy * c.o.y + uz * c.o.z;
		return rc;
	}
	coord3 operator * (const quaternion& q) const
	{
		coord3 rc = *this;
		rc.ux = q * ux;
		rc.uy = q * uy;
		rc.uz = q * uz;
		return rc;
	}
	// Parallel projection
	static real pl_prj(crvec v, crvec ax1, crvec ax2)
	{
		vec3 ax = ax1.cross(ax2);ax.norm();
		real co = ax1.dot(ax2);
		real si = sqrt(1 - co * co);
		real sc = (co / si);
		return (v.dot(ax1) - v.cross(ax1).dot(ax) * sc);
	} 
	// 向量向坐标系投影 注意：要保证ux,uy,uz是单位向量！
	friend vec3 operator / (crvec p, const coord3& c)
	{
		vec3 v = p - c.o;
		//{// 对于非正交情况
		//	return vec3(pl_prj(v-c.uz*v.dot(c.uz), c.ux, c.uy) / c.s.x, 
		//		    pl_prj(v-c.ux*v.dot(c.ux), c.uy, c.uz) / c.s.y, 
		//		    pl_prj(v-c.uy*v.dot(c.uy), c.uz, c.ux) / c.s.z);
		//}
		return vec3(v.dot(c.ux) / c.s.x, v.dot(c.uy) / c.s.y, v.dot(c.uz) / c.s.z);
	}
#define PL_PRJ3(v) vec3(pl_prj(v-c.uz*v.dot(c.uz), c.ux, c.uy) / c.s.x, \
			pl_prj(v-c.ux*v.dot(c.ux), c.uy, c.uz) / c.s.y, \
			pl_prj(v-c.uy*v.dot(c.uy), c.uz, c.ux) / c.s.z)
	coord3 operator / (const coord3& c) const
	{
		coord3 rc;
		/*{// 对于非正交情况
			rc.ux = PL_PRJ3(ux);
			rc.uy = PL_PRJ3(uy);
			rc.uz = PL_PRJ3(uz);
		}*/
		rc.ux = vec3(ux.dot(c.ux) / c.s.x, ux.dot(c.uy) / c.s.y, ux.dot(c.uz) / c.s.z);
		rc.uy = vec3(uy.dot(c.ux) / c.s.x, uy.dot(c.uy) / c.s.y, uy.dot(c.uz) / c.s.z);
		rc.uz = vec3(uz.dot(c.ux) / c.s.x, uz.dot(c.uy) / c.s.y, uz.dot(c.uz) / c.s.z);
		rc.o -= c.o;
		return rc;
	}
	void norm(bool bscl = true)
	{
#define ISZERO(a) (fabs(a) < 1e-10)
		s.x = ux.len(); if (!ISZERO(s.x)) ux /= s.x;
		s.y = uy.len(); if (!ISZERO(s.y)) uy /= s.y;
		s.z = uz.len(); if (!ISZERO(s.z)) uz /= s.z;

		if (!bscl)
			s = vec3::ONE;
	}
	vec3 sumvec() const
	{
		return (ux + uy + uz) * s;
	}
	// 本征值
	void eigenvalue() const
	{
		vec3 sv = sumvec();
		return ux.dot(sv) + uy.dot(sv) + uz.dot(sv);
	}
	real dot(crvec v) const
	{
		return v.dot(sumvec());
	}
	vec3 cross(const coord3& c) const
	{
		vec3 vx = VX();
		vec3 vy = VY();
		vec3 vz = VZ();

		vec3 cvx = c.VX();
		vec3 cvy = c.VY();
		vec3 cvz = c.VZ();

		return vec3(
			vy.dot(cvz) - vz.dot(cvy),
			vz.dot(cvx) - vx.dot(cvz),
			vx.dot(cvy) - vy.dot(cvx)
		);
	}
	// 曲率
	coord3 curvature(std::function<void(coord3& c, vec3 q)> coord_at, crvec q, crvec v)
	{
		const real delta_x = 0.01f;
		const real delta_y = 0.01f;
		vec3 q11 = q + vec3(delta_x, 0, 0);
		vec3 q21 = q + vec3(0, delta_y, 0);
		vec3 q12 = q + vec3(delta_x, delta_y, 0);

		coord3 c11;
		coord_at(c11, q11);
		vec3 p11 = q11 * c11;
		c11.norm(false);

		coord3 c21;
		coord_at(c21, q21);
		vec3 p21 = q21 * c21;
		c21.norm(false);

		coord3 c12;
		coord_at(c12, q12);
		vec3 p12 = q12 * c12;
		c12.norm(false);

		coord3 grad1 = c21 / c11; grad1.norm(false);
		coord3 grad2 = c12 / c11; grad2.norm(false);

		PRINT("--- g1 ---");
		coord3 g1 = grad1 * grad2; g1.norm(false); g1.dump();
		PRINT("--- g2 ---");
		coord3 g2 = grad2 * grad2; g2.norm(false); g2.dump();
		coord3 R = (g1 - g2);
		PRINT("--- R ---");
		R.dump();
		vec3 deta = v * R;
		deta.norm();
		PRINTVEC3(deta);

		return R;
	}
	void dump() const
	{
		//PRINT("-------");
		PRINT("ux: " << ux.x << "," << ux.y << "," << ux.z);
		PRINT("uy: " << uy.x << "," << uy.y << "," << uy.z);
		PRINT("uz: " << uz.x << "," << uz.y << "," << uz.z);
		PRINTVEC3(s);
		PRINTVEC3(o);
	}
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
		PRINT("rx: " << x << ", ry: " << y  << ", rz: " << z);
		return vec3(x, y, z);
	}
};

// **********************************************************************
// 梯度 / 时间变化率
// **********************************************************************
#define GRAD_V3(Fai, p, t) \
        vec3((Fai(p + vec3(deta_d,0.0,0.0), t) - Fai(p, t)) / deta_d,\
			 (Fai(p + vec3(0.0,deta_d,0.0), t) - Fai(p, t)) / deta_d, \
			 (Fai(p + vec3(0.0,0.0,deta_d), t) - Fai(p, t)) / deta_d)

#define GRAD_C3(A, p, t) \
    coord3( \
		(A(p + vec3(1.0,0.0,0.0) * deta_d, t) - A(p, t)) / deta_d, \
		(A(p + vec3(0.0,1.0,0.0) * deta_d, t) - A(p, t)) / deta_d, \
		(A(p + vec3(0.0,0.0,1.0) * deta_d, t) - A(p, t)) / deta_d)

#define DT_A(A, p, t) (A(p,t + deta_t) - A(p, t)) / deta_t
