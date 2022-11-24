/**
*					向量的定义是来自于四元数
*					向量不是完整的数，
*					向量跟空间结构有关系，
*					如果在时空中建议使用四元数
*/
// **********************************************************************
// 2D
// **********************************************************************
struct vector2 {
	union {
		real val[2];
		struct {
			real x;
			real y;
		};
	};

	static const vector2 ZERO;
	static const vector2 ONE;
	static const vector2 UX;
	static const vector2 UY;
	static const vector2 CENTER;

	real& operator [](int ind) {
		return val[ind];
	}

	vector2() {
		x = 0;
		y = 0;
	}
	explicit vector2(real v)
	{
		x = v;
		y = v;
	}
	vector2(real _x, real _y) {
		x = _x;
		y = _y;
	}
	vector2 operator + (const vector2& _p) const
	{
		vector2 fp;
		fp.x = x + _p.x;
		fp.y = y + _p.y;

		return fp;
	}
	void operator += (const vector2& _p)
	{
		x += _p.x;
		y += _p.y;
	}
	vector2 operator - (const vector2& _p) const
	{
		vector2 fp;
		fp.x = x - _p.x;
		fp.y = y - _p.y;
		return fp;
	}
	void operator -= (const vector2& _p)
	{
		x = x - _p.x;
		y = y - _p.y;
	}
	vector2 operator - () const
	{
		vector2 fp;
		fp.x = -x;
		fp.y = -y;
		return fp;
	}
	vector2 operator * (real s) const
	{
		vector2 fp;
		fp.x = s * x;
		fp.y = s * y;
		return fp;
	}
	void operator *= (real s)
	{
		x = s * x;
		y = s * y;
	}
	friend vector2 operator * (real s, const vector2& v)
	{
		vector2 fp;
		fp.x = v.x * s;
		fp.y = v.y * s;
		return fp;
	}
	vector2 operator * (const vector2& b) const
	{
		return vector2(x * b.x - y * b.y, x * b.y + y * b.x);
	}
	vector2 operator / (real s) const
	{
		vector2 fp;
		fp.x = x / s;
		fp.y = y / s;
		return fp;
	}
	void operator /= (real s)
	{
		x = x / s;
		y = y / s;
	}
	bool operator == (const vector2& rv) const
	{
		return (fabs(x - rv.x) <= 1e-5 && fabs(y - rv.y) <= 1e-5);
	}
	bool operator != (const vector2& rv) const
	{
		return (fabs(x - rv.x) > 1e-5 || fabs(y - rv.y) > 1e-5);
	}
	real len() const
	{
		return sqrt(x * x + y * y);
	}
	real sqrlen() const
	{
		return (x * x + y * y);
	}
	real angle() const
	{
		return atan2(y, x);
	}
	void norm()
	{
		real r = len();
		if (r > 0)
		{
			x /= r;
			y /= r;
		}
	}
	vector2 normcopy() const
	{
		real r = len();
		if (r > 0)
		{
			return vector2(x / r, y / r);
		}
		return vector2::ZERO;
	}
	void rot(real angle)
	{
		(*this) = (*this) * vector2::fromanglelength(angle, 1);
	}
	vector2 rotcopy(real angle) const
	{
		return (*this) * vector2::fromanglelength(angle, 1);
	}
	void rot(real angle, const vector2& o)
	{
		vector2 v = (*this) - o;
		v = v * vector2::fromanglelength(angle, 1);
		(*this) = v + o;
	}
	vector2 rotcopy(real angle, const vector2& o) const
	{
		vector2 v = (*this) - o;
		v = v * vector2::fromanglelength(angle, 1);
		return v + o;
	}
	vector2 rotcpy(real angle, const vector2& o) const
	{
		vector2 v = (*this) - o;
		v = v * vector2::fromanglelength(angle, 1);
		return v + o;
	}
	real dot(const vector2& v) const
	{
		return x * v.x + y * v.y;
	}
	real cross(const vector2& v) const
	{
		return x * v.y - y * v.x;
	}
	static vector2 fromanglelength(real _angle, real _r);
};

const vector2 vector2::ZERO = vector2(0, 0);
const vector2 vector2::ONE = vector2(1, 1);
const vector2 vector2::UX = vector2(1, 0);
const vector2 vector2::UY = vector2(0, 1);
const vector2 vector2::CENTER = vector2(0.5, 0.5);

vector2 vector2::fromanglelength(real _angle, real _r)
{
	return vector2(_r * cos(_angle), _r * sin(_angle));
}

// **********************************************************************
// 3D
// **********************************************************************
struct vector3
{
	static const vector3 ZERO;
	static const vector3 ONE;
	static const vector3 UX;
	static const vector3 UY;
	static const vector3 UZ;
	static const vector3 CENTER;
	static real sEPSINON;
	union {
		real val[3];
		struct
		{
			real x;
			real y;
			real z;
		};
		struct
		{
			real r;
			real g;
			real b;
		};
	};
	vector3()
	{
		x = 0;
		y = 0;
		z = 0;
	}
	explicit vector3(real v)
	{
		x = v;
		y = v;
		z = v;
	}
	vector3(real _x, real _y, real _z = 0.0f)
	{
		x = _x;
		y = _y;
		z = _z;
	}
	vector3(const vector3& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}
	real& operator [](int ind) {
		return val[ind];
	}
	vector2 xy() const
	{
		vector2 fp;
		fp.x = x;
		fp.y = y;
		return fp;
	}
	void xy(const vec2& _xy)
	{
		x = _xy.x;
		y = _xy.y;
	}
	vector2 xz() const
	{
		vector2 fp;
		fp.x = x;
		fp.y = z;
		return fp;
	}
	void xz(const vec2& _xz)
	{
		x = _xz.x;
		z = _xz.y;
	}
	vector2 yz() const
	{
		vector2 fp;
		fp.x = y;
		fp.y = z;
		return fp;
	}
	vector3 zxy() const
	{
		vector3 fp;
		fp.x = z;
		fp.y = x;
		fp.z = y;
		return fp;
	}
	vector3 operator + (const vector3& _p) const
	{
		vector3 fp;
		fp.x = x + _p.x;
		fp.y = y + _p.y;
		fp.z = z + _p.z;
		return fp;
	}
	void operator += (const vector3& _p)
	{
		x += _p.x;
		y += _p.y;
		z += _p.z;
	}
	vector3 operator - (const vector3& _p) const
	{
		vector3 fp;
		fp.x = x - _p.x;
		fp.y = y - _p.y;
		fp.z = z - _p.z;
		return fp;
	}
	void operator -= (const vector3& _p)
	{
		x -= _p.x;
		y -= _p.y;
		z -= _p.z;
	}
	vector3 operator - () const
	{
		vector3 fp;
		fp.x = -x;
		fp.y = -y;
		fp.z = -z;
		return fp;
	}
	vector3 operator * (real s) const
	{
		vector3 fp;
		fp.x = s * x;
		fp.y = s * y;
		fp.z = s * z;
		return fp;
	}
	vector3 operator * (const vector3& v) const
	{
		vector3 fp;
		fp.x = v.x * x;
		fp.y = v.y * y;
		fp.z = v.z * z;
		return fp;
	}
	friend vector3 operator * (real s, const vector3& v)
	{
		vector3 fp;
		fp.x = v.x * s;
		fp.y = v.y * s;
		fp.z = v.z * s;
		return fp;
	}
	void operator *= (real s)
	{
		x = s * x;
		y = s * y;
		z = s * z;
	}
	vector3 operator / (real s) const
	{
		vector3 fp;
		fp.x = x / s;
		fp.y = y / s;
		fp.z = z / s;
		return fp;
	}
	void operator /= (real s)
	{
		x = x / s;
		y = y / s;
		z = z / s;
	}
	bool operator == (const vector3& rv) const
	{
		return (fabs(x - rv.x) <= sEPSINON && fabs(y - rv.y) <= sEPSINON && fabs(z - rv.z) <= sEPSINON);
	}
	bool operator != (const vector3& rv) const
	{
		return (fabs(x - rv.x) > sEPSINON || fabs(y - rv.y) > sEPSINON || fabs(z - rv.z) > sEPSINON);
	}
	bool operator < (const vector3& rv) const
	{
		if (x < rv.x && y < rv.y && z < rv.z)
			return true;
		return false;
	}
	bool operator > (const vector3& rv) const
	{
		if (x > rv.x && y > rv.y && z > rv.z)
			return true;
		return false;
	}
	real len() const
	{
		return sqrt(x * x + y * y + z * z);
	}
	real lenxy() const
	{
		return sqrt(x * x + y * y);
	}
	real sqrlenxy() const
	{
		return (x * x + y * y);
	}
	real lenxz() const
	{
		return sqrt(x * x + z * z);
	}
	real sqrlen() const
	{
		return (x * x + y * y + z * z);
	}
	bool norm()
	{
		real r = len();
		if (r > 1e-6)
		{
			x /= r;
			y /= r;
			z /= r;
			return true;
		}
		return false;
	}
	vector3 normcopy() const
	{
		real r = len();
		if (r > 0)
		{
			return vector3(this->x / r,
				this->y / r,
				this->z / r);
		}
		return vector3::ZERO;
	}
	void fromanglelength(real _angle, real _r)
	{
		x = _r * cos(_angle);
		y = _r * sin(_angle);
	}
	inline real angle() const
	{
		return atan2(y, x);
	}
	inline real dot(const vector3& v) const
	{
		return x * v.x + y * v.y + z * v.z;
	}
	vec3 crossdot(const vector3& n) const
	{
		const vec3& v = (*this);
		return v - n * v.dot(n);
	}
	vector3 cross(const vector3& v) const
	{
		vector3 n;
		n.x = -(y * v.z - z * v.y);
		n.y = -(z * v.x - x * v.z);
		n.z = -(x * v.y - y * v.x);
		return n;
	}
	void rot(real angle, const vector3& ax);
	vector3 rotcopy(real angle, const vector3& ax) const;
	vector3 rotcpy(real angle, const vector3& ax) const;
	static vector3 rnd(real min = 0, real max = 1);

	static vector3 rndrad(real r = 1);
};

const vector3 vector3::ZERO = vector3(0, 0, 0);
const vector3 vector3::ONE = vector3(1, 1, 1);
const vector3 vector3::UX = vector3(1, 0, 0);
const vector3 vector3::UY = vector3(0, 1, 0);
const vector3 vector3::UZ = vector3(0, 0, 1);
const vector3 vector3::CENTER = vector3(0.0, 0.0, MAXZ / 2);
real vector3::sEPSINON = EPSINON;
// ----------------------------------------
vector3 vector3::rnd(real min, real max)
{
	return vector3(rrnd(min, max), rrnd(min, max), rrnd(min, max));
}
// ----------------------------------------
vector3 vector3::rndrad(real r)
{
	return rnd(-1, 1).normcopy() * r;
}
// **********************************************************************
// 4D
// **********************************************************************
struct vector4
{
	static const vector4 ZERO;
	static const vector4 UX;
	static const vector4 UY;
	static const vector4 UZ;
	static const vector4 UW;
	static const vector4 CENTER;
	union {
		real val[4];
		struct
		{
			float x;
			float y;
			float z;
			float w;
		};
	};
	vector4()
	{
		x = 0;
		y = 0;
		z = 0;
		w = 0;
	}
	vector4(float _x, float _y, float _z, float _w)
	{
		x = _x;
		y = _y;
		z = _z;
		w = _w;
	}
	vector4(float _x, crvec _v3)
	{
		x = _x;
		y = _v3.x;
		z = _v3.y;
		w = _v3.z;
	}
	vector4(float _v)
	{
		x = _v;
		y = _v;
		z = _v;
		w = _v;
	}
	real& operator [](int ind) {
		return val[ind];
	}
	vector4 operator + (const vector4& _p) const
	{
		vector4 fp;
		fp.x = x + _p.x;
		fp.y = y + _p.y;
		fp.z = z + _p.z;
		fp.w = w + _p.w;
		return fp;
	}
	void operator += (const vector4& _p)
	{
		x += _p.x;
		y += _p.y;
		z += _p.z;
		w += _p.w;
	}
	vector4 operator - (const vector4& _p) const
	{
		vector4 fp;
		fp.x = x - _p.x;
		fp.y = y - _p.y;
		fp.z = z - _p.z;
		fp.w = w - _p.w;
		return fp;
	}
	void operator -= (const vector4& _p)
	{
		x -= _p.x;
		y -= _p.y;
		z -= _p.z;
		w -= _p.w;
	}
	vector4 operator - () const
	{
		vector4 fp;
		fp.x = -x;
		fp.y = -y;
		fp.z = -z;
		fp.w = -w;
		return fp;
	}
	vector4 operator * (float s) const
	{
		vector4 fp;
		fp.x = s * x;
		fp.y = s * y;
		fp.z = s * z;
		fp.w = s * w;
		return fp;
	}
	void operator *= (float s)
	{
		x = s * x;
		y = s * y;
		z = s * z;
		w = s * w;
	}
	vector4 operator / (float s) const
	{
		vector4 fp;
		fp.x = x / s;
		fp.y = y / s;
		fp.z = z / s;
		fp.w = w / s;
		return fp;
	}
	void operator /= (float s)
	{
		x = x / s;
		y = y / s;
		z = z / s;
		w = w / s;
	}
	bool operator == (const vector4& rv) const
	{
		return (fabs(x - rv.x) <= EPSINON && fabs(y - rv.y) <= EPSINON && fabs(z - rv.z) <= EPSINON && fabs(w - rv.w) <= EPSINON);
	}
	bool operator != (const vector4& rv) const
	{
		return (fabs(x - rv.x) > EPSINON || fabs(y - rv.y) > EPSINON || fabs(z - rv.z) > EPSINON || fabs(w - rv.w) > EPSINON);
	}
	float len() const
	{
		return sqrt(x * x + y * y + z * z + w * w);
	}
	float sqrlen() const
	{
		return (x * x + y * y + z * z + w * w);
	}
	real norm()
	{
		float r = len();
		if (r > 1e-5)
		{
			x /= r;
			y /= r;
			z /= r;
			w /= r;
		}
		else
		{
			ERRORMSG("r == 0");
		}
		return r;
	}
	vector4 normcopy()
	{
		float r = len();
		if (r > 0)
		{
			return vector4(
				this->x / r,
				this->y / r,
				this->z / r,
				this->w / r
			);
		}
		else
		{
			ERRORMSG("r == 0");
		}
		return vector4(0, 0, 0, 0);
	}
	float dot(const vector4& v) const
	{
		return x * v.x + y * v.y + z * v.z + w * v.w;
	}

	// Cross4 computes the four-dimensional cross product of the three vectors
	// U, V and W, in that order. It returns the resulting four-vector.

	vector4 cross(const vector4& V, const vector4& W) const
	{
		vector4 result;
		double A, B, C, D, E, F;       // Intermediate Values

		// Calculate intermediate values.

		A = (V.x * W.y) - (V.y * W.x);
		B = (V.x * W.z) - (V.z * W.x);
		C = (V.x * W.w) - (V.w * W.x);
		D = (V.y * W.z) - (V.z * W.y);
		E = (V.y * W.w) - (V.w * W.y);
		F = (V.z * W.w) - (V.w * W.z);

		// Calculate the result-vector components.

		result.x = (y * F) - (z * E) + (w * D);
		result.y = -(x * F) + (z * C) - (w * B);
		result.z = (x * E) - (y * C) + (w * A);
		result.w = -(x * D) + (y * B) - (z * A);

		return result;
	}
};
// **********************************************************************
const vector4 vector4::ZERO = vector4(0, 0, 0, 0);
const vector4 vector4::UX = vector4(1, 0, 0, 0);
const vector4 vector4::UY = vector4(0, 1, 0, 0);
const vector4 vector4::UZ = vector4(0, 0, 1, 0);
const vector4 vector4::UW = vector4(0, 0, 0, 1);
const vector4 vector4::CENTER = vector4(0.5, 0.5, 0.5, 0.5);

// **********************************************************************
// nD
// 多维特征空间
// 来自于人类独特的思维方式（比如人鱼，半人马等）
// 把现象分解成多个维度的属性组合，然后再分别不同维度上进行融合，
// 最后可以合成出某种新的想象，如果在“本宇宙”得到验证，那么加强其权重，
// 否则削弱其权重。
// 注意，在物理学中（粒子物理学）至今没有观察到超过（3+1）自由度的高维粒子
// 所以说上层的高维的属性（比如尺寸与颜色）必然是某种非独立变量的某种权重组合。
// **********************************************************************
struct vectorn {
	std::vector<real> val;

	static const vectorn ZERO;
	static const vectorn ONE;
	static const vectorn CENTER;

	real& operator [](int ind) {
		return val[ind];
	}

	vectorn() {
	}
	explicit vectorn(real v)
	{
		for (int i = 0; i < val.size(); i++)
		{
			val[i] = v;
		}
	}
	vectorn operator + (vectorn& _p) const
	{
		vectorn fp;
		for (int i = 0; i < val.size(); i++)
		{
			fp[i] = val[i] + _p[i];
		}
		return fp;
	}
	void operator += (vectorn& _p)
	{
		for (int i = 0; i < val.size(); i++)
		{
			val[i] += _p[i];
		}
	}
	vectorn operator - (vectorn& _p) const
	{
		vectorn fp;
		for (int i = 0; i < val.size(); i++)
		{
			fp[i] = val[i] - _p[i];
		}
		return fp;
	}
	void operator -= (vectorn& _p)
	{
		for (int i = 0; i < val.size(); i++)
		{
			val[i] -= _p[i];
		}
	}
	vectorn operator - () const
	{
		vectorn fp;
		for (int i = 0; i < val.size(); i++)
		{
			fp[i] = - val[i];
		}
		return fp;
	}
	vectorn operator * (real s) const
	{
		vectorn fp;
		for (int i = 0; i < val.size(); i++)
		{
			fp[i] = val[i] * s;
		}
		return fp;
	}
	void operator *= (real s)
	{
		for (int i = 0; i < val.size(); i++)
		{
			val[i] *= s;
		}
	}
	friend vectorn operator * (real s,  vectorn& v)
	{
		vectorn fp;
		for (int i = 0; i < v.val.size(); i++)
		{
			fp[i] = v[i] * s;
		}
		return fp;
	}
	vectorn operator / (real s) const
	{
		vectorn fp;
		for (int i = 0; i < val.size(); i++)
		{
			fp[i] = val[i] / s;
		}
		return fp;
	}
	void operator /= (real s)
	{
		for (int i = 0; i < val.size(); i++)
		{
			val[i] /= s;
		}
	}
	bool operator == (const vectorn& v) const
	{
		if (val.size() != v.val.size())
			return false;
		bool ret = true;
		for (int i = 0; i < val.size(); i++)
		{
			ret &= (fabs(val[i] - v.val[i]) <= 1e-5);
			if (!ret) break;
		}
		return ret;
	}
	bool operator != (const vectorn& v) const
	{
		if (val.size() != v.val.size())
			return true;
		bool ret = false;
		for (int i = 0; i < val.size(); i++)
		{
			ret |= (fabs(val[i] - v.val[i]) <= 1e-5);
			if (ret) break;
		}
		return ret;
	}
	real len() const
	{
		real sqred = 0;
		for (int i = 0; i < val.size(); i++)
		{
			sqred += val[i] * val[i];
		}
		return sqrt(sqred);
	}
	real sqrlen() const
	{
		real sqred = 0;
		for (int i = 0; i < val.size(); i++)
		{
			sqred += val[i] * val[i];
		}
		return sqred;
	}
	void norm()
	{
		real r = len();
		if (r > 0)
		{
			(*this) /= r;
		}
	}
	vectorn normcopy() const
	{
		real r = len();
		if (r > 0)
		{
			return (*this) / r;
		}
		return vectorn::ZERO;
	}
	real dot(const vectorn& v) const
	{
		ASSERT(val.size() == v.val.size());
		real sum = 0;
		for (int i = 0; i < val.size(); i++)
		{
			sum += val[i] * v.val[i];
		}
		return sum;
	}
};
const vectorn vectorn::ZERO		= vectorn(0);
const vectorn vectorn::ONE		= vectorn(1);
const vectorn vectorn::CENTER	= vectorn(0.5);