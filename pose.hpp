/**
	A pose describes a position (see Position ) and an orientation (see Orientation ) in 3D space.
*/
struct pose3
{
	vec3 p;		// position
	quat q;		// orientation

	pose3 operator + (const pose3& _p) const
	{
		pose3 ret;
		ret.p = p + _p.p;
		ret.q = q; //  暂时不考虑_p的旋转
		return ret;
	}
	void operator += (const pose3& _p)
	{
		*this = (*this) + _p;
	}
	pose3 operator - (const pose3& _p) const
	{
		pose3 ret;
		ret.p = p - _p.p;
		ret.q = q; //  暂时不考虑_p的旋转
		return ret;
	}
	void operator -= (const pose3& _p)
	{
		*this = (*this) + _p;
	}
	pose3 operator * (const coord3& c) const
	{
		pose3 ret;
		ret.p = (*this).p * c;
		ret.q = (*this).q * c;
		return ret;
	}
	void operator *= (const coord3& c)
	{
		*this = (*this) * c;
	}
	pose3 operator * (const quaternion& q) const
	{
		pose3 ret;
		ret.p = (*this).p;
		ret.q = (*this).q * q;
		return ret;
	}
	void operator *= (const quaternion& q)
	{
		*this = (*this) * q;
	}
	pose3 operator / (const coord3& c) const
	{
		pose3 ret;
		ret.p = (*this).p / c;
		ret.q = (*this).q / c;
		return ret;
	}
	void operator /= (const coord3& c)
	{
		*this = (*this) / c;
	}
};
