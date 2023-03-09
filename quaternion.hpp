/****************************************************************************
							四元数
					(本文件基于OGRE引擎的数学库)
	四元数是在复数基础上的扩展,单位四元数用于旋转操作，向量是源自于四元数，
	不过二者有差别。目前四元数跟向量之间的关系以及应用存在着争议。

	*  *  *  *  *  *  *  *  *  详解  *  *  *  *  *  *  *  *  *  *  *  *  * 
	类似于复数，四元数也拥有指数形式：e^q，结果也是一个四元数： q = e^q,
	在底层物理里应了规范变换，跟坐标系变换有一些不同，规范变换更加
	侧重相位拥有时间属性，坐标系变换偏向于空间的变换以及曲率等特征提取。

	四元数存在归一化(normalise)，共轭(conj)，求逆(invert)，乘除法等操作，
	还规定了单位一（ONE).

****************************************************************************/
struct  quaternion
{
	real w, x, y, z;
	static const quaternion ONE;
	//-----------------------------------------------------------------------
	quaternion(
		real fW = 1.0,
		real fX = 0.0, real fY = 0.0, real fZ = 0.0)
	{
		w = fW;
		x = fX;
		y = fY;
		z = fZ;
	}
	quaternion(const quaternion& rkQ)
	{
		w = rkQ.w;
		x = rkQ.x;
		y = rkQ.y;
		z = rkQ.z;
	}
	quaternion(const real& rfAngle, const vector3& rkAxis)
	{
		this->fromangleaxis(rfAngle, rkAxis);
	}
	//-----------------------------------------------------------------------
	quaternion operator+ (const quaternion& rkQ) const
	{
		return quaternion(w + rkQ.w, x + rkQ.x, y + rkQ.y, z + rkQ.z);
	}
	//-----------------------------------------------------------------------
	quaternion operator- (const quaternion& rkQ) const
	{
		return quaternion(w - rkQ.w, x - rkQ.x, y - rkQ.y, z - rkQ.z);
	}
	quaternion operator - () const
	{
		quaternion q;
		q.x = -x;
		q.y = -y;
		q.z = -z;
		q.w = -w;
		return q;
	}
	//-----------------------------------------------------------------------
	vector3 operator* (const vector3& v) const
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
	vector3 friend operator* (const vector3& v, const quaternion& q)
	{
		return q * v;
	}
	void friend  operator*= (vector3& v, const quaternion& q)
	{
		v = q * v;
	}
	quaternion operator* (const quaternion& rkQ) const
	{
		// NOTE:  Multiplication is not generally commutative, so in most
		// cases p*q != q*p.

		return quaternion
		(
			w * rkQ.w - x * rkQ.x - y * rkQ.y - z * rkQ.z,
			w * rkQ.x + x * rkQ.w + y * rkQ.z - z * rkQ.y,
			w * rkQ.y + y * rkQ.w + z * rkQ.x - x * rkQ.z,
			w * rkQ.z + z * rkQ.w + x * rkQ.y - y * rkQ.x
		);
	}
	void operator*= (const quaternion& rkQ)
	{
		(*this) = (*this) * rkQ;
	}
	quaternion operator* (real fScalar) const
	{
		return quaternion(fScalar * w, fScalar * x, fScalar * y, fScalar * z);
	}
	//-----------------------------------------------------------------------
	quaternion friend operator* (real fScalar, const quaternion& rkQ)
	{
		return quaternion(fScalar * rkQ.w, fScalar * rkQ.x, fScalar * rkQ.y,
			fScalar * rkQ.z);
	}
	//-----------------------------------------------------------------------
	quaternion operator / (real fScalar) const
	{
		return quaternion(w / fScalar, x / fScalar, y / fScalar, z / fScalar);
	}
	void operator /= (real fScalar)
	{
		w /= fScalar;
		x /= fScalar;
		y /= fScalar;
		z /= fScalar;
	}
	quaternion operator / (const quaternion& q) const
	{
		return (*this) * q.conjcopy();
	}
	vector3 friend operator / (const vector3& v, const quaternion& q)
	{
		return q.conjcopy() * v;
	}
	//-----------------------------------------------------------------------
	bool operator == (const quaternion& rkQ) const
	{
		return w == rkQ.w && x == rkQ.x && y == rkQ.y && z == rkQ.z;
	}
	//-----------------------------------------------------------------------
	real dot(const quaternion& rkQ) const
	{
		return w * rkQ.w + x * rkQ.x + y * rkQ.y + z * rkQ.z;
	}
	//-----------------------------------------------------------------------
	real length() const
	{
		return sqrt(w * w + x * x + y * y + z * z);
	}
	//-----------------------------------------------------------------------
	vec3 xyz() const
	{
		return vec3(x, y, z);
	}
	vec3 axis() const
	{
		return vec3(x, y, z).normcopy();
	}
	real angle() const
	{
		return acos(w) * 2;
	}
	//-----------------------------------------------------------------------
	void normalize(void)
	{
		real len = length();
		if (len != 0)
		{
			real factor = 1.0f / (len);
			*this = *this * factor;
		}
	}
	quaternion normalized(void)
	{
		real len = length();
		if (len != 0)
		{
			return (*this) / len;
		}
	}
	//-----------------------------------------------------------------------
	void conj() // 共轭
	{
		this->w = w;
		this->x = -x; this->y = -y; this->z = -z;
	}
	quaternion conjcopy() const
	{
		quaternion q;
		q.w = w;
		q.x = -x; q.y = -y; q.z = -z;
		return q;
	}
	//-----------------------------------------------------------------------
	// 指数上运算
	quaternion exp() const
	{
		real r = length();
		vec3 v(x, y, z);
		real th = v.len();
		vec3 n = v.normcopy();
		vec3 qv = n * (r * sin(th));
		return quaternion(r * cos(th), qv.x, qv.y, qv.z);
	}
	friend quaternion exp(const quaternion& q)
	{
		real r = q.length();
		vec3 v(q.x, q.y, q.z);
		real th = v.len();
		vec3 n = v.normcopy();
		vec3 qv = n * (r * sin(th));
		return quaternion(r * cos(th), qv.x, qv.y, qv.z);
	}
	quaternion operator ^ (int n) const
	{
		quaternion ret = *this;
		for (int i = 1; i < n; i++)
		{
			ret = ret * (*this);
		}
		return ret;
	}
	quaternion operator ^ (real t)
	{
		return slerp(t, quaternion::ONE, *this, false);
	}
	//-----------------------------------------------------------------------
	// v1, v2 是单位向量
	void fromvectors(crvec v1, crvec v2)
	{
		fromangleaxis(acos(v1.dot(v2)), v1.cross(v2).normlized());
	}
	//-----------------------------------------------------------------------
	void fromangleaxis(real rfAngle, const vector3& rkAxis)
	{
		// assert:  axis[] is unit length
		//
		// The quaternion representing the rotation is
		//   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

		real fHalfAngle(0.5 * rfAngle);
		real fSin = sin(fHalfAngle);
		w = cos(fHalfAngle);
		x = fSin * rkAxis.x;
		y = fSin * rkAxis.y;
		z = fSin * rkAxis.z;
	}
	//-----------------------------------------------------------------------
	void fromeuler(real roll, real pitch, real yaw)
	{
		real t0 = cos(yaw * 0.5);
		real t1 = sin(yaw * 0.5);
		real t2 = cos(roll * 0.5);
		real t3 = sin(roll * 0.5);
		real t4 = cos(pitch * 0.5);
		real t5 = sin(pitch * 0.5);

		w = t2 * t4 * t0 + t3 * t5 * t1;
		x = t3 * t4 * t0 - t2 * t5 * t1;
		y = t2 * t5 * t0 + t3 * t4 * t1;
		z = t2 * t4 * t1 - t3 * t5 * t0;
	}
	//-----------------------------------------------------------------------
	vec3 toeuler() const
	{
		vec3 v;

		real epsilon = 0.00001f;
		real halfpi = 0.5 * PI;

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
			v.y = -atan2(x * z + y * w, 0.5 - x * x - y * y);
			v.z = -atan2(x * y + z * w, 0.5 - x * x - z * z);
		}
		return v;
	}
	//-----------------------------------------------------------------------
	// 求导运算
	quaternion derivative_angle(const quaternion& a, const quaternion& b) const
	{
		return quaternion(1, (b / a).axis()); //角度为单位1
	}
	quaternion diff_angle(real angle) const
	{
		return quaternion(angle, axis());
	}
	quaternion diff_angle(const quaternion& deri, real angle) const
	{
		return quaternion(angle, deri.axis());
	}
	//-----------------------------------------------------------------------
	static quaternion slerp(const quaternion& qa, const quaternion& qb, double t) {
		// quaternion to return
		quaternion qm;
		// Calculate angle between them.
		double cosHalfTheta = qa.w * qb.w + qa.x * qb.x + qa.y * qb.y + qa.z * qb.z;
		// if qa=qb or qa=-qb then theta = 0 and we can return qa
		if (abs(cosHalfTheta) >= 1.0) {
			qm.w = qa.w; qm.x = qa.x; qm.y = qa.y; qm.z = qa.z;
			return qm;
		}
		// Calculate temporary values.
		double halfTheta = acos(cosHalfTheta);
		double sinHalfTheta = sqrt(1.0 - cosHalfTheta * cosHalfTheta);
		// if theta = 180 degrees then result is not fully defined
		// we could rotate around any axis normal to qa or qb
		if (fabs(sinHalfTheta) < 0.001) { // fabs is floating point absolute
			qm.w = (qa.w * 0.5 + qb.w * 0.5);
			qm.x = (qa.x * 0.5 + qb.x * 0.5);
			qm.y = (qa.y * 0.5 + qb.y * 0.5);
			qm.z = (qa.z * 0.5 + qb.z * 0.5);
			return qm;
		}
		double ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
		double ratioB = sin(t * halfTheta) / sinHalfTheta;
		//calculate Quaternion.
		qm.w = (qa.w * ratioA + qb.w * ratioB);
		qm.x = (qa.x * ratioA + qb.x * ratioB);
		qm.y = (qa.y * ratioA + qb.y * ratioB);
		qm.z = (qa.z * ratioA + qb.z * ratioB);
		return qm;
	}
	//-----------------------------------------------------------------------
	static quaternion slerp(real fT, const quaternion& rkP,
		const quaternion& rkQ, bool shortestPath)
	{
		const real msEpsilon = 1e-03;
		real fCos = rkP.dot(rkQ);
		quaternion rkT;

		// Do we need to invert rotation?
		if (fCos < 0.0f && shortestPath)
		{
			fCos = -fCos;
			rkT = -rkQ;
		}
		else
		{
			rkT = rkQ;
		}

		if (fabs(fCos) < 1 - msEpsilon)
		{
			// Standard case (slerp)
			real fSin = sqrt(1 - (fCos * fCos));
			real fAngle = atan2(fSin, fCos);
			real fInvSin = 1.0f / fSin;
			real fCoeff0 = sin((1.0f - fT) * fAngle) * fInvSin;
			real fCoeff1 = sin(fT * fAngle) * fInvSin;
			return fCoeff0 * rkP + fCoeff1 * rkT;
		}
		else
		{
			// There are two situations:
			// 1. "rkP" and "rkQ" are very close (fCos ~= +1), so we can do a linear
			//    interpolation safely.
			// 2. "rkP" and "rkQ" are almost inverse of each other (fCos ~= -1), there
			//    are an infinite number of possibilities interpolation. but we haven't
			//    have method to fix this case, so just use linear interpolation here.
			quaternion t = (1.0f - fT) * rkP + fT * rkT;
			// taking the complement requires renormalisation
			t.normalize();
			return t;
		}
	}
	//-----------------------------------------------------------------------
	quaternion nlerp(real fT, const quaternion& rkP,
		const quaternion& rkQ, bool shortestPath)
	{
		quaternion result;
		real fCos = rkP.dot(rkQ);
		if (fCos < 0.0f && shortestPath)
		{
			result = rkP + fT * ((-rkQ) - rkP);
		}
		else
		{
			result = rkP + fT * (rkQ - rkP);
		}
		result.normalize();
		return result;
	}
};
// **********************************************************************
const quaternion quaternion::ONE = quaternion(1, 0, 0, 0);