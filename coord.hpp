/**
	Coordinate structure define 
*/

struct coord
{
	vec3 ux = vec3::UX;			// 方向
	vec3 uy = vec3::UY;
	vec3 uz = vec3::UZ;

	vec3 scl = vec3::ONE;			// 缩放

	vec3 o;					// 空间位置
	
	coord() {}
	coord(crvec _ux, crvec _uy, crvec _uz)
	{
		ux = _ux; uy = _uy; uz = _uz;
	}
	
	vec3 UX0() const { return ux * scl.x; }
	vec3 UY0() const { return uy * scl.y; }
	vec3 UZ0() const { return uz * scl.z; }
	
	vec3 (funAX)();
	funAX UX = UX0;
	funAX UY = UY0;
	funAX UZ = UZ0;
	
	vec3 operator * (crvec p) const
	{
		return ux * (scl.x * p.x) + uy * (scl.y * p.y) + uz * (scl.z * p.z) + o;
	}
	friend vec3 operator * (crvec p, const coord& c)
	{
		return c.ux * (c.scl.x * v.x) + c.uy * (c.scl.y * v.y) + c.uz * (c.scl.z * v.z) + c.o;
	}
	coord operator * (coord& c) const
	{
		coord_t rc;
		rc.ux = ux * c.ux.x + uy * c.ux.y + uz * c.ux.z;
		rc.uy = ux * c.uy.x + uy * c.uy.y + uz * c.uy.z;
		rc.uz = ux * c.uz.x + uy * c.uz.y + uz * c.uz.z;
		rc.scl = scl * c.scl;
		rc.o = o + ux * c.o.x + uy * c.o.y + uz * c.o.z;
		return rc;
	}
	friend vec3 operator / (crvec p, const coord& c)
	{
		vec3 v = p - c.o;
		return vec3(v.dot(c.ux)/ c.scl.x, v.dot(c.uy) / c.scl.y, v.dot(c.uz) / c.scl.z);
	}
	coord operator / (const coord& c)
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
	vec3 eigenvec()
	{
		return (ux + uy + uz) * scl;
	}
	real dot(crvec v)
	{
		return v.dot(eigenvec());
	}
	vec3 cross(const coord3& c)
	{
		return vec3(
			uy.dot(c.UZ()) - uz.dot(c.UY()) +
			uz.dot(c.UX()) - ux.dot(c.UZ()) +
			ux.dot(c.UY()) - uy.dot(c.UX())
		);
	}
	coord3 deta_coord3(vec3 deta)
	{
		coord3 rc(UX(),UY(),UZ());
		rc.scl = deta;
	}
	void rot(real ang, crvec ax)
	{
		ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);
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
	void dump()
	{
		PRINT("-------");
		PRINT("ux: " << ux.x << "," << ux.y << "," << ux.z);
		PRINT("uy: " << uy.x << "," << uy.y << "," << uy.z);
		PRINT("uz: " << uz.x << "," << uz.y << "," << uz.z);
	}
};
