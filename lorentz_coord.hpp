/**
 * Lorentz Coordinate System
 *
 * The Lorentz Coordinate System is a coordinate system designed to describe space-time transformations and handle the time dimension.
 * It is an extension of the complex coordinate system, where the rotation of the coordinate system can be represented by a quaternion composed of a vector and a complex angle.
 * The Lorentz Coordinate System consists of three vectors, ux, uy, and uz, which represent the x, y, and z-axis directions respectively.
 * Additionally, it includes a complex number power, which represents the power in the time dimension.
 * The Lorentz Coordinate System allows for operations such as multiplication, division, cross product, transpose, reverse, flip, rotation, and conversion to Euler angles.
 * These operations enable the transformation of vectors between different coordinate systems and the calculation of relativistic effects.
 * The design philosophy of the Lorentz Coordinate System is to provide a flexible and intuitive way to handle space-time transformations and relativistic effects.
 * By incorporating complex numbers and quaternions, it allows for a more comprehensive representation of rotations and transformations in 3D space.
 */

struct lorentz_coord
{
	static const lorentz_coord ZERO;
	static const lorentz_coord ONE;
	vec3 ux = vec3::UX; // x-axis direction
	vec3 uy = vec3::UY; // y-axis direction
	vec3 uz = vec3::UZ; // z-axis direction
	complex power; // power in the time dimension

	lorentz_coord() {}

	/**
	 * Constructor
	 * @param c - lorentz_coord object
	 */
	lorentz_coord(const lorentz_coord& c)
	{
		ux = c.ux; uy = c.uy; uz = c.uz;
		power = c.power;
	}

	/**
	 * Constructor
	 * @param _ux - x-axis direction
	 * @param _uy - y-axis direction
	 * @param _uz - z-axis direction
	 * @param _power - power in the time dimension
	 */
	lorentz_coord(const vec3& _ux, const vec3& _uy, const vec3& _uz, const complex& _power)
	{
		ux = _ux; uy = _uy; uz = _uz;
		power = _power;
	}

	/**
	 * Constructor
	 * @param _ux - x-axis direction
	 * @param _uy - y-axis direction
	 * @param _power - power in the time dimension
	 */
	lorentz_coord(const vec3& _ux, const vec3& _uy, const complex& _power)
	{
		ux = _ux; uy = _uy; uz = ux.cross(uy);
		power = _power;
	}

	/**
	 * Constructor
	 * @param ang - angle
	 * @param ax - axis
	 * @param _power - power in the time dimension
	 */
	lorentz_coord(real ang, const vec3& ax, const complex& _power)
	{
		ux.rot(ang, ax);
		uy.rot(ang, ax);
		uz.rot(ang, ax);
		power = _power;
	}

	/**
	 * Constructor
	 * @param pit - pitch
	 * @param yaw - yaw
	 * @param rol - roll
	 * @param _power - power in the time dimension
	 */
	lorentz_coord(real pit, real yaw, real rol, const complex& _power)
	{
		ux.rot(pit, vec3::UX);
		uy.rot(yaw, vec3::UY);
		uz.rot(rol, vec3::UZ);
		power = _power;
	}

	/**
	 * Constructor
	 * @param q - quaternion
	 * @param _power - power in the time dimension
	 */
	lorentz_coord(const quaternion& q, const complex& _power)
	{
		ux = q * vec3::UX;
		uy = q * vec3::UY;
		uz = q * vec3::UZ;
		power = _power;
	}

	/**
	 * Multiplication: define a vector in the coordinate system or restore a vector to the parent space
	 * @param p - vector
	 * @param c - lorentz_coord object
	 * @return result of the multiplication
	 */
	friend vec3 operator * (const vec3& p, const lorentz_coord& c)
	{  
		return (c.ux * p.x + c.uy * p.y + c.uz * p.z) * cos(c.power);
	}

	/**
	 * Multiplication: Cchild * Cparent * ...
	 * @param c - lorentz_coord object
	 * @return result of the multiplication
	 */
	lorentz_coord operator * (const lorentz_coord& c) const
	{
		lorentz_coord rc;
		rc.ux = ux.x * c.ux + ux.y * c.uy + ux.z * c.uz;
		rc.uy = uy.x * c.ux + uy.y * c.uy + uy.z * c.uz;
		rc.uz = uz.x * c.ux + uz.y * c.uy + uz.z * c.uz;
		rc.power = power + c.power;
		return rc;
	}

	/**
	 * Multiplication: q * C
	 * @param q - quaternion
	 * @param c - lorentz_coord object
	 * @return result of the multiplication
	 */
	friend quaternion operator * (const quaternion& q, const lorentz_coord& c)
	{
		return q * c.to_quaternion();
	}

	/**
	 * Multiplication: C * q
	 * @param q - quaternion
	 * @return result of the multiplication
	 */
	lorentz_coord operator * (const quaternion& q) const
	{
		lorentz_coord rc;
		rc.ux = q * ux;
		rc.uy = q * uy;
		rc.uz = q * uz;
		rc.power = power + q.angle();
		return rc;
	}

	/**
	 * Division: vector projection onto the coordinate system
	 * @param v - vector
	 * @param c - lorentz_coord object
	 * @return result of the division
	 */
	friend vec3 operator / (const vec3& v, const lorentz_coord& c)
	{
		return vec3(v.dot(c.ux), v.dot(c.uy), v.dot(c.uz)) * cos(-c.power);
	}

	/**
	 * Division: C1 * C2^-1
	 * @param c - lorentz_coord object
	 * @return result of the division
	 */
	lorentz_coord operator / (const lorentz_coord& c) const
	{
		lorentz_coord rc;
		rc.ux = vec3(ux.dot(c.ux), ux.dot(c.uy), ux.dot(c.uz));
		rc.uy = vec3(uy.dot(c.ux), uy.dot(c.uy), uy.dot(c.uz));
		rc.uz = vec3(uz.dot(c.ux), uz.dot(c.uy), uz.dot(c.uz));
		rc.power = power - c.power;
		return rc;
	}

	/**
	 * Division: q / C
	 * @param q - quaternion
	 * @param c - lorentz_coord object
	 * @return result of the division
	 */
	friend quaternion operator / (const quaternion& q, const lorentz_coord& c)
	{
		return q * c.to_quaternion().conjugate();
	}

	/**
	 * Division: C / q
	 * @param q - quaternion
	 * @return result of the division
	 */
	lorentz_coord operator / (const quaternion& q) const
	{
		return (*this) * q.conjugate();
	}

	/**
	 * Cross product: C1 x C2
	 * @param c - lorentz_coord object
	 * @return result of the cross product
	 */
	lorentz_coord cross(const lorentz_coord& c) const
	{
		return lorentz_coord(
			ux.cross(c.uz) - ux.cross(c.uy),
			uy.cross(c.ux) - uy.cross(c.uz),
			uz.cross(c.uy) - uz.cross(c.ux),
			power + c.power
		);
	}

	/**
	 * Cross product: C x v
	 * @param v - vector
	 * @return result of the cross product
	 */
	lorentz_coord cross(const vec3& v) const
	{
		return lorentz_coord(
			ux.cross(v),
			uy.cross(v),
			uz.cross(v),
			power
		);
	}

	/**
	 * Transpose (axis exchange)
	 */
	void transpose()
	{
		vec3 _ux = vec3(ux.x, uy.x, uz.x);
		vec3 _uy = vec3(ux.y, uy.y, uz.y);
		vec3 _uz = vec3(ux.z, uy.z, uz.z);
		ux = _ux; uy = _uy; uz = _uz;
	}

	/**
	 * Transpose (axis exchange)
	 * @return transposed lorentz_coord object
	 */
	lorentz_coord transposed() const
	{
		lorentz_coord c = (*this);
		c.ux = vec3(ux.x, uy.x, uz.x);
		c.uy = vec3(ux.y, uy.y, uz.y);
		c.uz = vec3(ux.z, uy.z, uz.z);
		return c;
	}

	/**
	 * Reverse
	 */
	void reverse()
	{
		(*this) = ONE / (*this);
	}

	/**
	 * Reverse
	 * @return reversed lorentz_coord object
	 */
	lorentz_coord reversed() const
	{
		return ONE / (*this);
	}

	/**
	 * Flip X-axis
	 */
	void flipX()
	{
		ux = -ux;
	}

	/**
	 * Flip Y-axis
	 */
	void flipY()
	{
		uy = -uy;
	}

	/**
	 * Flip Z-axis
	 */
	void flipZ()
	{
		uz = -uz;
	}

	/**
	 * Rotation
	 * @param ang - angle
	 * @param ax - axis
	 */
	void rotate(real ang, const vec3& ax)
	{
		ux.rotate(ang, ax);
		uy.rotate(ang, ax);
		uz.rotate(ang, ax);
	}

	/**
	 * Convert coordinate system to Euler angles
	 * @return Euler angles
	 */
	vec3 coord_to_eulers() const
	{
		const lorentz_coord& rm = *this;
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

	/**
	 * Dot product with a vector
	 * @param v - vector
	 * @return dot product
	 */
	real dot(const vec3& v) const
	{
		return v.dot(ux) + v.dot(uy) + v.dot(uz);
	}

	/**
	 * Dot product with another coordinate system
	 * @param c - lorentz_coord object
	 * @return dot product
	 */
	real dot(const lorentz_coord& c) const
	{
		return c.ux.dot(ux) + c.uy.dot(uy) + c.uz.dot(uz);
	}

	/**
	 * Gradient coordinate system = Gradient X Tangent space
	 * @param c1 - initial lorentz_coord object
	 * @param c2 - transformed lorentz_coord object
	 * @return gradient lorentz_coord object
	 */
	lorentz_coord gradient(const lorentz_coord& c1, const lorentz_coord& c2)
	{
		return c1.reversed() * c2;
	}

	/**
	 * Convert to quaternion
	 * @return quaternion
	 */
	quaternion to_quaternion() const
	{
		vec3 pyr = coord_to_eulers();
		quaternion q;
		q.from_euler(pyr.x, pyr.y, pyr.z);
		return q;
	}

	/**
	 * Lerp operator
	 * @param v - vector
	 * @return lerped lorentz_coord object
	 */
	lorentz_coord operator ^ (const vec3& v) const
	{
		lorentz_coord c = *this;
		c.ux = vec3::lerp(vec3::UX, c.ux, v.x); c.ux.normalize();
		c.uy = vec3::lerp(vec3::UY, c.uy, v.y); c.uy.normalize();
		c.uz = vec3::lerp(vec3::UZ, c.uz, v.z); c.uz.normalize();
		return c;
	}

	/**
	 * Lerp operator
	 * @param f - factor
	 * @return lerped lorentz_coord object
	 */
	lorentz_coord operator ^ (real f) const
	{
		return lorentz_coord((*this).to_quaternion() ^ f);
	}
};

const lorentz_coord lorentz_coord::ZERO = lorentz_coord(vec3::ZERO, vec3::ZERO, vec3::ZERO, complex::ZERO);
const lorentz_coord lorentz_coord::ONE = lorentz_coord(vec3::UX, vec3::UY, vec3::UZ, complex::ONE);
