/********************************************************************
*				坐标系
* *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
* 坐标系类是我单独封装，用于简化坐标变换，衍生出许多算法，能解决一些
* 坐标系变换相关的问题。
* 
* 坐标系是物理规范的简化，物理学因为规范而变得数学上很困难，如果能用
* 算法简化规范运算,可以还原一个初中生眼中的物理世界！
* 物理学范式: 测量 = 规范 * 本征 * 量纲
* 
*/

const real deta_d = 0.0001f;	// 空间微分精度
const real deta_t = 0.0001f;	// 时间微分精度
extern void edgeax(const VECLIST& e, vec& ux, vec& uy, vec& uz);
struct coord3
{
	vec3 ux = vec3::UX;		// 方向
	vec3 uy = vec3::UY;
	vec3 uz = vec3::UZ;

	vec3 scl = vec3::ONE;	// 缩放

	vec3 o;					// 原点

	coord3() {}
	coord3(const coord3& c)
	{
		ux = c.ux; uy = c.uy; uz = c.uz;
		scl = c.scl;
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
		coord3 cy(pit, vec3::UY);
		coord3 cz(pit, vec3::UZ);
		*this = cx * cy * cz;
	}
	coord3(const VECLIST& e)
	{
		edgeax(e, ux, uy, uz);
	}
	vec3 VX() const { return ux * scl.x; }
	vec3 VY() const { return uy * scl.y; }
	vec3 VZ() const { return uz * scl.z; }

	void rot(real ang, crvec ax)
	{
		//	o.rot(ang, ax);
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
		if (is_same_dirs(c))
		{
			rc.ux = ux; rc.uy = uy; rc.uz = uz;
		}
		else
		{
			rc.ux = VX() - c.VX();
			rc.uy = VY() - c.VY();
			rc.uz = VZ() - c.VZ();
			rc.norm();
		}
		rc.o = o - c.o;
		return rc;
	}
	vec3 operator * (crvec p) const
	{
		return ux * (scl.x * p.x) + uy * (scl.y * p.y) + uz * (scl.z * p.z) + o;
	}
	friend vec3 operator * (crvec p, const coord3& c)
	{
		return c.ux * (c.scl.x * p.x) + c.uy * (c.scl.y * p.y) + c.uz * (c.scl.z * p.z) + c.o;
	}
	coord3 operator * (const coord3& c) const
	{
		coord3 rc;
		rc.ux = ux * c.ux.x + uy * c.ux.y + uz * c.ux.z;
		rc.uy = ux * c.uy.x + uy * c.uy.y + uz * c.uy.z;
		rc.uz = ux * c.uz.x + uy * c.uz.y + uz * c.uz.z;
		rc.scl = scl * c.scl;
		rc.o = o + ux * c.o.x + uy * c.o.y + uz * c.o.z;
		return rc;
	}
	friend vec3 operator / (crvec p, const coord3& c)
	{
		vec3 v = p - c.o;
		return vec3(v.dot(c.ux) / c.scl.x, v.dot(c.uy) / c.scl.y, v.dot(c.uz) / c.scl.z);
	}
	coord3 operator / (const coord3& c) const
	{
		coord3 rc;
		rc.ux = vec3(ux.dot(c.ux) / c.scl.x, ux.dot(c.uy) / c.scl.y, ux.dot(c.uz) / c.scl.z);
		rc.uy = vec3(uy.dot(c.ux) / c.scl.x, uy.dot(c.uy) / c.scl.y, uy.dot(c.uz) / c.scl.z);
		rc.uz = vec3(uz.dot(c.ux) / c.scl.x, uz.dot(c.uy) / c.scl.y, uz.dot(c.uz) / c.scl.z);
		rc.o -= c.o;
		return rc;
	}
	void norm()
	{
#define ISZERO(a) (fabs(a) < 1e-6)
		scl.x = ux.len(); if (!ISZERO(scl.x)) ux /= scl.x;
		scl.y = uy.len(); if (!ISZERO(scl.y)) uy /= scl.y;
		scl.z = uz.len(); if (!ISZERO(scl.z)) uz /= scl.z;
	}
	// 本征向量
	vec3 eigenvec() const
	{
		return (ux + uy + uz) * scl;
	}
	real dot(crvec v) const
	{
		return v.dot(eigenvec());
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
			vy.dot(cvz) - vz.dot(cvy) +
			vz.dot(cvx) - vx.dot(cvz) +
			vx.dot(cvy) - vy.dot(cvx)
		);
	}
	coord3 curvature(std::function<void(coord3& c, vec3 q)> coord_at, crvec q, crvec v)
	{
		vec3 q11 = q;
		vec3 q21 = q + vec3(deta_d, 0, 0);
		vec3 q12 = q + vec3(0, deta_d, 0);
		vec3 q22 = q + vec3(deta_d, deta_d, 0);

		coord3 c11;
		coord_at(c11, q11);
		vec3 p11 = q11 * c11;

		coord3 c21;
		coord_at(c21, q21);
		vec3 p21 = q21 * c21;

		coord3 c12;
		coord_at(c12, q12);
		vec3 p12 = q12 * c12;

		coord3 grad1 = c21 / c11;
		coord3 grad2 = c12 / c11;

		vec3 deta = v * grad1 * grad2 - v * grad2 * grad1; // 非阿贝尔群
		deta /= deta_d;

		PRINTVEC3(v * grad1); PRINTVEC3(v * grad2);
		PRINTVEC3(deta);

		coord3 ret = grad1 * grad2 - grad2 * grad1;
		return ret;
	}
	void dump() const
	{
		PRINT("-------");
		PRINT("ux: " << ux.x << "," << ux.y << "," << ux.z);
		PRINT("uy: " << uy.x << "," << uy.y << "," << uy.z);
		PRINT("uz: " << uz.x << "," << uz.y << "," << uz.z);
		PRINT("o: " << o.x << "," << o.y << "," << o.z);
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
		PRINT("rx: " << x * 180 / PI << ", ry: " << y * 180 / PI << ", rz: " << z * 180 / PI);
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
