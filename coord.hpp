/**
	Coordinate structure define 
*/

// **********************************************************************
// coord_t
// **********************************************************************
struct coord_t
{
	vec3 ux = vec3::UX;		// 方向
	vec3 uy = vec3::UY;
	vec3 uz = vec3::UZ;

	vec3 scl = vec3::ONE;	// 缩放

	vec3 o;					// 空间
	real t;					// 时间
	vec3 vel;				// 运动参考系
	
	coord_t() {}
	coord_t(crvec _ux, crvec _uy, crvec _uz)
	{
		ux = _ux; uy = _uy; uz = _uz;
	}
	void rot(real ang, crvec ax)
	{
	//	o.rot(ang, ax);
		ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);
	}
	// C1*C2*C3* ... *Cloc * Vloc （transfrom)
	vec3 operator * (crvec v) const
	{
		return ux * (scl.x * v.x) + uy * (scl.y * v.y) + uz * (scl.z * v.z) + o;
	}
	// V * C1 * C2 ...
	friend vec3 operator * (crvec v, const coord_t& c)
	{
		return c.ux * (c.scl.x * v.x) + c.uy * (c.scl.y * v.y) + c.uz * (c.scl.z * v.z) + c.o;
	}
	coord_t operator * (coord_t& c) const
	{
		coord_t rc;
		rc.ux = ux * c.ux.x + uy * c.ux.y + uz * c.ux.z;
		rc.uy = ux * c.uy.x + uy * c.uy.y + uz * c.uy.z;
		rc.ux = ux * c.uz.x + uy * c.uz.y + uz * c.uz.z;
		rc.scl = scl * c.scl;
		rc.o += ux * c.o.x + uy * c.o.y + uz * c.o.z;
		return rc;
	}
	// Vworld/C1/C2/C3/ ... /Cloc（projection)
	friend vec3 operator / (crvec v, const coord_t& c)
	{
		vec3 dv = v - c.o;
		return vec3(dv.dot(c.ux)/ c.scl.x, dv.dot(c.uy) / c.scl.y, dv.dot(c.uz) / c.scl.z);
	}
	coord_t operator / (const coord_t& c)
	{
		coord_t rc;
		rc.ux = vec3(ux.dot(c.ux) / c.scl.x, ux.dot(c.uy) / c.scl.y, ux.dot(c.uz) / c.scl.z);
		rc.uy = vec3(uy.dot(c.ux) / c.scl.x, uy.dot(c.uy) / c.scl.y, uy.dot(c.uz) / c.scl.z);
		rc.uz = vec3(uz.dot(c.ux) / c.scl.x, uz.dot(c.uy) / c.scl.y, uz.dot(c.uz) / c.scl.z);
		rc.o -= c.o;
		return rc;
	}
	void norm()
	{
		scl.x = ux.len(); ux /= scl.x;
		scl.y = uy.len(); uy /= scl.y;
		scl.z = uz.len(); uz /= scl.z;
	}
	real dot(crvec v)
	{
		vec3 dv = v - o;
		return (dv.dot(ux) * scl.x + dv.dot(uy) * scl.y + dv.dot(uz) * scl.z;
	}
	vec3 cross(crvec v)
	{
		vec3 dv = v - o;
		return dv.cross(ux) * scl.x + dv.cross(uy) * scl.y + dv.cross(uz) * scl.z;
	}
	void dump()
	{
		PRINT("-------");
		PRINT("ux: " << ux.x << "," << ux.y << "," << ux.z);
		PRINT("uy: " << uy.x << "," << uy.y << "," << uy.z);
		PRINT("uz: " << uz.x << "," << uz.y << "," << uz.z);
	}
	static vec3 coord2eulers(coord_t& rm)
	{
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
		return vec3(x, y, z);
	}
};
