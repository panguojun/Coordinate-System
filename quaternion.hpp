/**************************************************************************
                         Quaternion

	Quaternions are extensions of complex numbers, with unit quaternions used for rotation operations,
	specifically designed for rotation transformations in coordinate system theory.

	Quaternions support normalization, conjugation, inverse, and arithmetic operations,
	and are an important component of frame field composite operator theory.

**************************************************************************/

// Forward declarations for real and vector3
#ifndef REAL_DEFINED
typedef float real;
#define REAL_DEFINED
#endif

// Forward declaration
struct vector3;

struct quaternion
{
	static const quaternion ONE;
	static const quaternion UX;
	static const quaternion UY;
	static const quaternion UZ;

	real w = 1, x = 0, y = 0, z = 0;

	//-----------------------------------------------------------------------
	// Constructors
	quaternion() { }
	quaternion(real fW, real fX, real fY, real fZ)
	{
		w = fW;
		x = fX;
		y = fY;
		z = fZ;
	}
	quaternion(real pitch, real yaw, real roll)
	{
		from_eulers(pitch, yaw, roll);
	}
	quaternion(const vector3& pyr)
	{
		from_eulers(pyr.x, pyr.y, pyr.z);
	}
	quaternion(const quaternion& rkQ)
	{
		w = rkQ.w;
		x = rkQ.x;
		y = rkQ.y;
		z = rkQ.z;
	}
	quaternion(real rfAngle, const vector3& rkAxis)
	{
		ang_axis(rfAngle, rkAxis);
	}
	quaternion(const vector3& v1, const vector3& v2)
	{
		from_vectors(v1, v2);
	}

	//-----------------------------------------------------------------------
	// Arithmetic operations
	quaternion operator+(const quaternion& rkQ) const
	{
		return quaternion(w + rkQ.w, x + rkQ.x, y + rkQ.y, z + rkQ.z);
	}

	quaternion operator-(const quaternion& rkQ) const
	{
		return quaternion(w - rkQ.w, x - rkQ.x, y - rkQ.y, z - rkQ.z);
	}
	quaternion operator-() const
	{
		return quaternion(-w, -x, -y, -z);
	}

	//-----------------------------------------------------------------------
	// Vector rotation (core application of quaternions)
	vector3 operator*(const vector3& v) const
	{
		// nVidia SDK implementation
		vector3 uv, uuv;
		vector3 qvec(x, y, z);
		uv = qvec.cross(v);
		uuv = qvec.cross(uv);
		uv = uv * (2.0f * w);
		uuv = uuv * 2.0f;

		return v + uv + uuv;
	}
	vector3 friend operator*(const vector3& v, const quaternion& q)
	{
		return q * v;
	}
	void friend operator*=(vector3& v, const quaternion& q)
	{
		v = q * v;
	}

	// Quaternion multiplication (rotation composition)
	quaternion operator*(const quaternion& rkQ) const
	{
		return quaternion
		(
			w * rkQ.w - x * rkQ.x - y * rkQ.y - z * rkQ.z,
			w * rkQ.x + x * rkQ.w + y * rkQ.z - z * rkQ.y,
			w * rkQ.y + y * rkQ.w + z * rkQ.x - x * rkQ.z,
			w * rkQ.z + z * rkQ.w + x * rkQ.y - y * rkQ.x
		);
	}
	void operator*=(const quaternion& rkQ)
	{
		(*this) = (*this) * rkQ;
	}
	quaternion operator*(real fScalar) const
	{
		return quaternion(fScalar * w, fScalar * x, fScalar * y, fScalar * z);
	}
	quaternion friend operator*(real fScalar, const quaternion& rkQ)
	{
		return quaternion(fScalar * rkQ.w, fScalar * rkQ.x, fScalar * rkQ.y,
			fScalar * rkQ.z);
	}

	//-----------------------------------------------------------------------
	// Division operations
	quaternion operator/(real fScalar) const
	{
		return quaternion(w / fScalar, x / fScalar, y / fScalar, z / fScalar);
	}
	void operator/=(real fScalar)
	{
		w /= fScalar;
		x /= fScalar;
		y /= fScalar;
		z /= fScalar;
	}
	quaternion operator/(const quaternion& q) const
	{
		return (*this) * q.conjcopy();
	}
	vector3 friend operator/(const vector3& v, const quaternion& q)
	{
		return q.conjcopy() * v;
	}

	//-----------------------------------------------------------------------
	// Comparison operations
	bool operator==(const quaternion& rkQ) const
	{
		const real eps = 1e-5f;
		return (fabs(w - rkQ.w) < eps) &&
			(fabs(x - rkQ.x) < eps) &&
			(fabs(y - rkQ.y) < eps) &&
			(fabs(z - rkQ.z) < eps);
	}
	bool operator!=(const quaternion& rkQ) const
	{
		return !(*this == rkQ);
	}

	//-----------------------------------------------------------------------
	// Basic properties
	real dot(const quaternion& rkQ) const
	{
		return w * rkQ.w + x * rkQ.x + y * rkQ.y + z * rkQ.z;
	}
	real length() const
	{
		return sqrt(w * w + x * x + y * y + z * z);
	}
	vector3 xyz() const
	{
		return vector3(x, y, z);
	}
	vector3 axis() const
	{
		return vector3(x, y, z).normalized();
	}
	real angle() const
	{
		if (w >= 1.0f)
			return 0.0f;
		if (w <= -1.0f)
			return 3.14159265359f;
		real ang = acos(w) * 2;
		if (ang > 3.14159265359f)
			return ang - 3.14159265359f * 2;
		if (ang < -3.14159265359f)
			return ang + 3.14159265359f * 2;
		return ang;
	}

	//-----------------------------------------------------------------------
	// Normalization
	quaternion normalize(void)
	{
		real len = length();
		if (len != 0)
		{
			real factor = 1.0f / len;
			*this = *this * factor;
		}
		return *this;
	}
	quaternion normalized(void) const
	{
		real len = length();
		if (len == 0)
		{
			return quaternion::ONE;
		}
		return (*this) / len;
	}

	//-----------------------------------------------------------------------
	// Conjugation (inverse operation for coordinate system transformations)
	void conj()
	{
		this->x = -x; this->y = -y; this->z = -z;
	}
	quaternion conjcopy() const
	{
		return quaternion(w, -x, -y, -z);
	}

	//-----------------------------------------------------------------------
	// Inverse
	quaternion inverse() const {
		real lenSquared = w * w + x * x + y * y + z * z;
		if (lenSquared != 0) {
			real factor = 1.0f / lenSquared;
			return conjcopy() * factor;
		}
		return quaternion::ONE;
	}

	//-----------------------------------------------------------------------
	// Construct quaternion from two vectors (foundation of frame field transformations)
	quaternion from_vectors(const vector3& v1, const vector3& v2)
	{
		const real eps = 1e-5f;
		if (fabs(v1.x - v2.x) <= eps && fabs(v1.y - v2.y) <= eps && fabs(v1.z - v2.z) <= eps)
		{
			(*this) = quaternion::ONE;
		}
		else
		{
			real dot = v1.dot(v2);
			if (fabs(dot + 1.0f) < eps)	// Handle 180 degree case
			{
				vector3 ax;
				vector3 uz = vector3::UZ();
				if (fabs(v1.x) < eps && fabs(v1.y) < eps)
					uz = -vector3::UX();
				ax = uz.cross(v1).normalized();
				ang_axis(3.14159265359f, ax.normalized());
			}
			else if (dot > -1.0f + eps && dot <= 1.0f - eps) // Handle general case
			{
				vector3 axis = v1.cross(v2).normalized();
				real angle = acos(dot);
				ang_axis(angle, axis);
			}
		}
		return (*this);
	}

	//-----------------------------------------------------------------------
	// Angle-axis definition (standard representation for frame field rotations)
	quaternion ang_axis(real rfAngle, const vector3& rkAxis)
	{
		real fHalfAngle = 0.5f * rfAngle;
		real fSin = sin(fHalfAngle);
		w = cos(fHalfAngle);
		x = fSin * rkAxis.x;
		y = fSin * rkAxis.y;
		z = fSin * rkAxis.z;

		return (*this);
	}

	//-----------------------------------------------------------------------
	// Euler angle conversion
	void from_eulers(real roll, real pitch, real yaw)
	{
		real t0 = cos(yaw * 0.5f);
		real t1 = sin(yaw * 0.5f);
		real t2 = cos(roll * 0.5f);
		real t3 = sin(roll * 0.5f);
		real t4 = cos(pitch * 0.5f);
		real t5 = sin(pitch * 0.5f);

		w = t2 * t4 * t0 + t3 * t5 * t1;
		x = t3 * t4 * t0 - t2 * t5 * t1;
		y = t2 * t5 * t0 + t3 * t4 * t1;
		z = t2 * t4 * t1 - t3 * t5 * t0;
	}
	vector3 to_eulers() const
	{
		vector3 v;
		real epsilon = 0.00001f;
		real halfpi = 0.5f * 3.14159265359f;

		real temp = 2 * (y * z - x * w);
		if (temp >= 1 - epsilon)
		{
			v.x = halfpi;
			v.y = -atan2(y, w);
			v.z = -atan2(z, w);
		}
		else if (-temp >= 1 - epsilon)
		{
			v.x = -halfpi;
			v.y = -atan2(y, w);
			v.z = -atan2(z, w);
		}
		else
		{
			v.x = asin(temp);
			v.y = -atan2(x * z + y * w, 0.5f - x * x - y * y);
			v.z = -atan2(x * y + z * w, 0.5f - x * x - z * z);
		}
		return v;
	}

	//-----------------------------------------------------------------------
	// Exponential operation (note operator precedence)
	quaternion operator^(real t)
	{
		return slerp(quaternion::ONE, *this, t);
	}

	//-----------------------------------------------------------------------
	// Spherical linear interpolation (smooth frame field transitions)
	static quaternion slerp(const quaternion& qa, const quaternion& qb, real t) {
		quaternion qm;
		real cosHalfTheta = qa.w * qb.w + qa.x * qb.x + qa.y * qb.y + qa.z * qb.z;

		if (fabs(cosHalfTheta) >= 1.0f) {
			qm.w = qa.w; qm.x = qa.x; qm.y = qa.y; qm.z = qa.z;
			return qm;
		}

		real halfTheta = acos(cosHalfTheta);
		real sinHalfTheta = sqrt(1.0f - cosHalfTheta * cosHalfTheta);

		if (fabs(sinHalfTheta) < 0.001f) {
			qm.w = (qa.w * 0.5f + qb.w * 0.5f);
			qm.x = (qa.x * 0.5f + qb.x * 0.5f);
			qm.y = (qa.y * 0.5f + qb.y * 0.5f);
			qm.z = (qa.z * 0.5f + qb.z * 0.5f);
			return qm;
		}

		real ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
		real ratioB = sin(t * halfTheta) / sinHalfTheta;

		qm.w = (qa.w * ratioA + qb.w * ratioB);
		qm.x = (qa.x * ratioA + qb.x * ratioB);
		qm.y = (qa.y * ratioA + qb.y * ratioB);
		qm.z = (qa.z * ratioA + qb.z * ratioB);
		return qm;
	}
};

// Constant definition - declared first, avoid circular dependency initialization issues
inline const quaternion quaternion::ONE = quaternion(1, 0, 0, 0);