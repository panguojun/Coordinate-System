/*********************************************************************
*				坐标系
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
* 	坐标系类是我单独封装，用于简化坐标变换，衍生出许多算法，能解决一些
* 	坐标系变换相关的问题。
* 	坐标系的运算跟李群很相似。
* 
*  *  *  *  *  *  *  *  *  *  详解  *  *  *  *  *  *  *  *  *  *  *  *  
* 	坐标系本体符号 C，坐标系之间的变换可以写成G（G=C1/C2,GRAD梯度的意思）
* 	坐标系的李括号: [C1,C2] = C1*C2 - C2*C1

* 	具体来说：
* 	定义一个内禀坐标系(假设它是平直空间，向量可以随意移动而不变)下V,在弯
* 	曲坐标系下观察V，不同点上V是不同的，故而坐标系跟位置有关，取相邻两点
* 	（1),(2)点处有向量V1,V2，对应坐标系C1,C2，那么：
*		V = V1 * C1 = V2 * C2 => 
*		V2 = V1 * C1 / C2, 令 G12 = C1 / C2 =>
*		V2 = V1 * G12
* 
* 	在弯曲坐标系下内禀坐标系x,y轴的平行线投影得到的u,v曲线上G12分别在两个
* 	方向上对应Gu,Gv, 从(u1,v1)到(u2,v2) 计算两个路径的差别再加上修正项可
* 	得曲率公式为：
*		Ruv = Gu*Gv - Gv*Gu*Gu*Wu*Gv*Wv
*		W = (U + V*Gu) - (V + U*Gv) 
*/

//#define	Parallel_Projection	  // 非正交坐标系下平行投影
// *******************************************************************
//  |_
// C     2d Coordinate System
// *******************************************************************
struct coord2
{
	static const coord2 ONE;

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
	coord2 operator - (const coord2& c) const
	{
		coord2 rc;
		rc.ux = VX() - c.VX();
		rc.uy = VY() - c.VY();
		rc.norm();
		rc.o = o - c.o;
		return rc;
	}
	// 在坐标系下定义一个向量
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
		rc.ux = ux.x * c.ux + uy.x * c.ux;
		rc.uy = ux.y * c.uy + uy.y * c.uy;
		rc.s = s * c.s;
		rc.o = o + ux * c.o.x + uy * c.o.y;
		return rc;
	}
#ifdef Parallel_Projection
	// 非正交坐标系下平行投影 Parallel projection
	static real pl_dot(crvec2 v, crvec2 ax1, crvec2 ax2)
	{
		real co = ax1.dot(ax2);
		real si = sqrt(1 - co * co);
		real sc = (co / si);
		return (v.dot(ax1)- v.cross(ax1) * sc);
	}
#endif
	// 向量向坐标系投影 注意：要保证ux,uy是单位向量！
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
	void revert()
	{
		(*this) = ONE / (*this);
	}
	coord2 revertcopy() const
	{
		return ONE / (*this);
	}
	// 梯度坐标系
	static coord2 gradcoord(const coord2& c1, const coord2& c2)
	{
		return c1.revertcopy() * c2;
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
		ux = _ux; uy = _uy; uz = ux.cross(uy);
	}
	coord3(crvec _p)
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
	void fromvectors(crvec v1, crvec v2)
	{
		quaternion q;
		q.fromvectors(v1,v2);
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
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

	// 归一化的正交坐标系
	coord3 ucoord() const
	{
		coord3 c = *this;
		c.norm(false);
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
		rc.ux = VX() - c.VX();
		rc.uy = VY() - c.VY();
		rc.uz = VZ() - c.VZ();
		rc.norm();
		rc.o = o - c.o;
		return rc;
	}
	// 乘法：在坐标系下定义一个向量
	friend vec3 operator * (crvec p, const coord3& c)
	{
		return c.ux * (c.s.x * p.x) + c.uy * (c.s.y * p.y) + c.uz * (c.s.z * p.z) + c.o;
	}
	coord3 operator * (crvec p) const
	{
		coord3 c = *this;
		c.ux = lerp(vec3::UX, c.VX(), p.x);
		c.uy = lerp(vec3::UY, c.VY(), p.y);
		c.uz = lerp(vec3::UZ, c.VZ(), p.z);
		c.norm();
		c.o.x *= p.x; c.o.y *= p.y; c.o.z *= p.z;
		return c;
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

	// 除法：向量向坐标系投影 注意：要保证ux,uy,uz是单位向量！
#ifdef Parallel_Projection
	// 非正交坐标系下平行投影 Parallel projection
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
		{// 对于非正交情况
			return vec3(
					pl_prj(v-c.uz*v.dot(c.uz), c.ux, c.uy) / c.s.x, 
				    pl_prj(v-c.ux*v.dot(c.ux), c.uy, c.uz) / c.s.y, 
				    pl_prj(v-c.uy*v.dot(c.uy), c.uz, c.ux) / c.s.z);
		}
#endif
		return vec3(v.dot(c.ux) / c.s.x, v.dot(c.uy) / c.s.y, v.dot(c.uz) / c.s.z);
	}
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
		rc.ux = vec3(ux.dot(c.ux) / c.s.x, ux.dot(c.uy) / c.s.y, ux.dot(c.uz) / c.s.z);
		rc.uy = vec3(uy.dot(c.ux) / c.s.x, uy.dot(c.uy) / c.s.y, uy.dot(c.uz) / c.s.z);
		rc.uz = vec3(uz.dot(c.ux) / c.s.x, uz.dot(c.uy) / c.s.y, uz.dot(c.uz) / c.s.z);
#endif
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
	// 转置
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
	void reverse()
	{
		(*this) = ONE / (*this);
	}
	coord3 reversed() const
	{
		return ONE / (*this);
	}
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
	// 本征向量（坐标系作为旋转变换时候的特征）
	vec3 eigenvec() const
	{
		return toquat().axis();
	}
	real dot(crvec v) const
	{
		return v.dot(ux) * s.x + v.dot(uy) * s.y + v.dot(uz) * s.z;
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
	// 坐标系到欧拉角，要保证是归一化的正交坐标系
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
	// 梯度坐标系 = 梯度 X 切空间
	// 相当于一阶坐标系的导数
	static coord3 gradcoord(const coord3& c1, const coord3& c2)
	{
		return c1.reversed() * c2;
	}
	// 曲率测试算法
	coord3 curvature(std::function<void(coord3& c, vec3 q)> coord_at, crvec q)
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

		coord3 grad1 = coord3::gradcoord(c11, c21); grad1.norm(false);
		coord3 grad2 = coord3::gradcoord(c11, c12); grad2.norm(false);

		PRINT("--- g1 ---");
		coord3 g1 = grad1 * grad2; g1.norm(false); g1.dump();
		PRINT("--- g2 ---");
		coord3 g2 = grad2 * grad2; g2.norm(false); g2.dump();
		coord3 R = (g1 - g2); // 这个是（二阶）李括号，就是在坐标系映射下的李括号，在[X,Y]=0时可以对应曲率
		//R.dump();
	/*	vec3 deta = v * R;
		deta.norm();
		PRINTVEC3(deta);*/

		return R;
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
const coord3 coord3::ZERO = {0};
const coord3 coord3::ONE = coord3();
