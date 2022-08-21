/**
	Coordinate structure define 
*/

// **********************************************************************
// coord_t
// **********************************************************************
struct coord_t
{
	vec3 ux = vec3::UX;			// 方向
	vec3 uy = vec3::UY;
	vec3 uz = vec3::UZ;

	vec3 scl = vec3::ONE;			// 缩放

	vec3 o;					// 空间位置
	
	coord_t() {}
	coord_t(crvec _ux, crvec _uy, crvec _uz)
	{
		ux = _ux; uy = _uy; uz = _uz;
	}
	void rot(real ang, crvec ax)
	{
		ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);
	}
	vec3 operator * (crvec p) const
	{
		return ux * (scl.x * p.x) + uy * (scl.y * p.y) + uz * (scl.z * p.z) + o;
	}
	friend vec3 operator * (crvec p, const coord_t& c)
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
	friend vec3 operator / (crvec p, const coord_t& c)
	{
		vec3 v = p - c.o;
		return vec3(v.dot(c.ux)/ c.scl.x, v.dot(c.uy) / c.scl.y, v.dot(c.uz) / c.scl.z);
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
		return (v.dot(ux) * scl.x + v.dot(uy) * scl.y + v.dot(uz) * scl.z;
	}
	real dot(const coord_t& c)
	{
		return 	(ux.dot(c.ux) * c.scl.x + ux.dot(c.uy) * c.scl.y + ux.dot(c.uz) * c.scl.z) +
			(uy.dot(c.ux) * c.scl.x + uy.dot(c.uy) * c.scl.y + uy.dot(c.uz) * c.scl.z) +
			(uz.dot(c.ux) * c.scl.x + uz.dot(c.uy) * c.scl.y + uz.dot(c.uz) * c.scl.z);
	}
	vec3 cross(crvec v)
	{
		return v.cross(ux) * scl.x + v.cross(uy) * scl.y + v.cross(uz) * scl.z;
	}
	vec3 cross(const coord_t& c)
	{
		return 	(ux.cross(c.ux) * c.scl.x + ux.cross(c.uy) * c.scl.y + ux.cross(c.uz) * c.scl.z) +
			(uy.cross(c.ux) * c.scl.x + uy.cross(c.uy) * c.scl.y + uy.cross(c.uz) * c.scl.z) +
			(uz.cross(c.ux) * c.scl.x + uz.cross(c.uy) * c.scl.y + uz.cross(c.uz) * c.scl.z);
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
