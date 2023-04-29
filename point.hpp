// ***************************************************************************
//				【Point】
//									
//			 主要指的是整数点，在整数空间内
// ***************************************************************************
#pragma once
struct point2
{
	int x = 0, y = 0;
	point2()
	{
		x = 0;
		y = 0;
	}
	point2(int ix, int iy)
	{
		x = ix; y = iy;
	}
	int len()
	{
		return sqrt((x) * (x)+(y) * (y));
	}
	vec2 tovec() const
	{
		return vec2(x, y);
	}
	bool operator == (const point2& rp) const
	{
		return (x == rp.x && y == rp.y);
	}
	bool operator != (const point2& rp) const
	{
		return (x != rp.x || y != rp.y);
	}
	point2 operator + (const point2& rp) const
	{
		point2 p;
		p.x = x + rp.x;
		p.y = y + rp.y;
		return  p;
	}
	point2 operator - (const point2& rp) const
	{
		point2 p;
		p.x = x - rp.x;
		p.y = y - rp.y;
		return  p;
	}
	point2 operator / (real factor) const
	{
		point2 p;
		p.x = x / factor;
		p.y = y / factor;
		return  p;
	}
	static real dis(const point2& p1, const point2& p2)
	{
		return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
	}
	static real Mdis(const point2& p1, const point2& p2)
	{
		return abs(p1.x - p2.x) + abs(p1.y - p2.y);
	}
	static real dis(int x1, int y1, int x2, int y2)
	{
		return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	}
	static int idis(int x1, int y1, int x2, int y2)
	{
		return (int)sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	}
	static real angle(int x1, int x2, int y1, int y2)
	{
		vec2 v1(x1, y1), v2(x2, y2);
		real ang = v1.normcopy().dot(v2.normcopy());
		return ang;
	}
	vec2 tovec()
	{
		vec2 p;
		p.x = x;
		p.y = y;
		return p;
	}
	/*void fromvec(const vec2& v)
	{
		x = v.x;
		y = v.y;
	}*/
	static point2 fromvec(const vec2& v)
	{
		point2 p;
		p.x = v.x;
		p.y = v.y;
		return p;
	}
};

// ***********************************************
// 3D point
// ***********************************************
struct point3
{
	static const point3 ZERO;
	static const point3 ONE;
	int x = 0, y = 0, z = 0;
	point3()
	{
		x = 0;
		y = 0;
		z = 0;
	}
	point3(int ix, int iy, int iz)
	{
		x = ix; y = iy; z = iz;
	}
	int llen() const
	{
		return ((x) * (x)+(y) * (y)+(z) * (z));
	}

	bool operator == (const point3& rp) const
	{
		return (x == rp.x && y == rp.y && z == rp.z);
	}
	point3 operator + (const point3& rp) const
	{
		point3 p;
		p.x = x + rp.x;
		p.y = y + rp.y;
		p.z = z + rp.z;
		return  p;
	}
	point3 operator - (const point3& rp) const
	{
		point3 p;
		p.x = x - rp.x;
		p.y = y - rp.y;
		p.z = z - rp.z;
		return  p;
	}
	point3 operator * (real factor) const
	{
		point3 p;
		p.x = x * factor;
		p.y = y * factor;
		p.z = z * factor;
		return  p;
	}
	point3 operator / (real factor) const
	{
		point3 p;
		p.x = x / factor;
		p.y = y / factor;
		p.z = z / factor;
		return  p;
	}
	vec3 tovec() const
	{
		return vec3(x, y, z);
	}
};
const point3 point3::ZERO = point3(0,0,0);
const point3 point3::ONE = point3(1,1,1);

struct point3f
{
	real x = 0, y = 0, z = 0;
	point3f()
	{
		x = 0;
		y = 0;
		z = 0;
	}
	point3f(real _x, real _y, real _z)
	{
		x = _x; y = _y; z = _z;
	}
	int llen() const
	{
		return ((x) * (x)+(y) * (y)+(z) * (z));
	}

	bool operator == (const point3f& rp) const
	{
		return (x == rp.x && y == rp.y && z == rp.z);
	}
	point3f operator + (const point3f& rp) const
	{
		point3f p;
		p.x = x + rp.x;
		p.y = y + rp.y;
		p.z = z + rp.z;
		return  p;
	}
	point3f operator - (const point3f& rp) const
	{
		point3f p;
		p.x = x - rp.x;
		p.y = y - rp.y;
		p.z = z - rp.z;
		return  p;
	}
	point3f operator * (real factor) const
	{
		point3f p;
		p.x = x * factor;
		p.y = y * factor;
		p.z = z * factor;
		return  p;
	}
	point3f operator / (real factor) const
	{
		point3f p;
		p.x = x / factor;
		p.y = y / factor;
		p.z = z / factor;
		return  p;
	}
};
// ***********************************************
// pointn
// ***********************************************
struct pointn {
	std::vector<int> val;
	static const pointn ZERO;

	int& operator [](int index) {
		if (val.empty())
			val = { 0,0,0,0 }; // 默认四个维度

		ASSERT(val.size() > index);
		return val[index];
	}
	int& operator [](const std::string& key) {
		if (key == "x")
		{
			return (*this)[0];
		}
		else if (key == "y")
		{
			return (*this)[1];
		}
		else if (key == "z")
		{
			return (*this)[2];
		}
		else if (key == "w")
		{
			return (*this)[3];
		}
		return (*this)[0];
	}
	pointn() {
	}
	explicit pointn(int v)
	{
		val = { v };
	}
	explicit pointn(int x, int y)
	{
		val = { x,y };
	}
	// dimsz : SZX,SZY,SZZ,...
	int toindex(const std::vector<int>& dimsz) const
	{
		ASSERT(val.size() == dimsz.size());
		int sz = 1;
		int index = 0;
		for (int i = 0; i < val.size(); i++)
		{
			index += val[i] * sz;
			sz *= dimsz[i];
		}
	}
	int toindex(int dimsz) const
	{
		int sz = 1;
		int index = 0;
		for (int i = 0; i < val.size(); i++)
		{
			index += val[i] * sz;
			sz *= dimsz;
		}
		return index;
	}
	pointn operator + (const pointn& _p) const
	{
		pointn fp;
		int steps = _MIN(val.size(), _p.val.size());
		for (int i = 0; i < steps; i++)
		{
			fp[i] = val[i] + const_cast<pointn&>(_p)[i];
		}
		return fp;
	}
	void operator += (const pointn& _p)
	{
		int steps = _MIN(val.size(), _p.val.size());
		for (int i = 0; i < steps; i++)
		{
			val[i] += const_cast<pointn&>(_p)[i];
		}
	}
	pointn operator - (const pointn& _p) const
	{
		pointn fp;
		int steps = _MIN(val.size(), _p.val.size());
		for (int i = 0; i < steps; i++)
		{
			fp[i] = val[i] - const_cast<pointn&>(_p)[i];
		}
		return fp;
	}
	void operator -= (const pointn& _p)
	{
		int steps = _MIN(val.size(), _p.val.size());
		for (int i = 0; i < steps; i++)
		{
			val[i] -= const_cast<pointn&>(_p)[i];
		}
	}
	pointn operator - () const
	{
		pointn fp;
		for (int i = 0; i < val.size(); i++)
		{
			fp[i] = -val[i];
		}
		return fp;
	}
	pointn operator * (int s) const
	{
		pointn fp;
		for (int i = 0; i < val.size(); i++)
		{
			fp[i] = val[i] * s;
		}
		return fp;
	}
	void operator *= (int s)
	{
		for (int i = 0; i < val.size(); i++)
		{
			val[i] *= s;
		}
	}
	friend pointn operator * (int s, pointn& v)
	{
		pointn fp;
		for (int i = 0; i < v.val.size(); i++)
		{
			fp[i] = v[i] * s;
		}
		return fp;
	}
	pointn operator / (int s) const
	{
		pointn fp;
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
	bool operator == (const pointn& v) const
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
	bool operator != (const pointn& v) const
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

	int ilen() const
	{
		real sqred = 0;
		for (int i = 0; i < val.size(); i++)
		{
			sqred += val[i] * val[i];
		}
		return (int)sqrt(sqred);
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
	pointn normcopy() const
	{
		real r = len();
		if (r > 0)
		{
			return (*this) / r;
		}
		return pointn::ZERO;
	}
	real dot(const pointn& v) const
	{
		ASSERT(val.size() == v.val.size());
		real sum = 0;
		for (int i = 0; i < val.size(); i++)
		{
			sum += val[i] * v.val[i];
		}
		return sum;
	}
	bool check_inside(int min, int max) const
	{
		for (int i = 0; i < val.size(); i++)
		{
			if (val[i] < min || val[i] > max - 1)
				return false;
		}
		return true;
	}
};
// ***********************************************
const pointn pointn::ZERO = pointn(0);
