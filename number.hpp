/***************************************************************************
				数
			真正的数背后有数论支持
			整数，小数，复数，四元数都是数
			而向量，矩阵不是数
			尽量降低维度 使用三维运算
			解析数论里有大学问，等待挖掘利用！
***************************************************************************/
// 复数 
// 黎曼函数Z(s) = 1^(-s) + 2^(-s) + 3^(-s) + ...
// ************************************************************************
namespace number_math
{
	struct complex
	{
		real x = 0, y = 0;

		static const complex ZERO;
		static const complex ONE;
		static const complex UX;
		static const complex UY;
		static const complex CENTER;

		complex() {
			x = 0;
			y = 0;
		}
		explicit complex(real r)
		{
			x = r;
			y = 0;
		}
		complex(real _x, real _y) {
			x = _x;
			y = _y;
		}
		void from_ra(real r, real a) {
			x = r * cos(a);
			y = r * sin(a);
		}
		friend complex operator ^ (real _r, const complex& c)
		{
			complex ret;
			real r = pow(_r, c.x);
			if (c.y == 0)
				return complex(r);

			real a = c.y;
			ret.from_ra(r, a*log(_r));
			return ret;
		}
		complex operator + (const complex& _p) const
		{
			complex fp;
			fp.x = x + _p.x;
			fp.y = y + _p.y;

			return fp;
		}
		void operator += (const complex& _p)
		{
			x += _p.x;
			y += _p.y;
		}
		complex operator - (const complex& _p) const
		{
			complex fp;
			fp.x = x - _p.x;
			fp.y = y - _p.y;
			return fp;
		}
		void operator -= (const complex& _p)
		{
			x = x - _p.x;
			y = y - _p.y;
		}
		complex operator - () const
		{
			complex fp;
			fp.x = -x;
			fp.y = -y;
			return fp;
		}
		complex operator * (real s) const
		{
			complex fp;
			fp.x = s * x;
			fp.y = s * y;
			return fp;
		}
		friend complex operator * (real s, const complex& v)
		{
			complex fp;
			fp.x = v.x * s;
			fp.y = v.y * s;
			return fp;
		}
		complex operator * (const complex& b) const
		{
			return complex(x * b.x - y * b.y, x * b.y + y * b.x);
		}
		complex operator / (const complex& c) const
		{
			complex ret;
			real rr = c.x * c.x + c.y * c.y;
			ret.x = (x * c.x + y * c.y) / rr;
			ret.y = (y * c.x + x * c.y) / rr;
			return ret;
		}
		friend complex operator / (real r, const complex& b)
		{
			return complex(r) / b;
		}
		complex operator / (real s) const
		{
			complex fp;
			fp.x = x / s;
			fp.y = y / s;
			return fp;
		}
		void operator /= (real s)
		{
			x = x / s;
			y = y / s;
		}
		void operator /= (const complex& c)
		{
			real rr = c.x * c.x + c.y * c.y;
			*this = complex((x * c.x + y * c.y) / rr, (y * c.x - x * c.y) / rr);
		}
		bool operator == (const complex& rv) const
		{
			return (fabs(x - rv.x) <= 1e-5 && fabs(y - rv.y) <= 1e-5);
		}
		bool operator != (const complex& rv) const
		{
			return (fabs(x - rv.x) > 1e-5 || fabs(y - rv.y) > 1e-5);
		}
		real len() const
		{
			return sqrt(x * x + y * y);
		}
		real angle() const
		{
			return atan2(y, x);
		}
		complex conj()// 共轭
		{
			complex c;
			c.x = x;
			c.y = -y;
			return c;
		}
		friend complex exp(const complex& c)
		{
			complex ret;
			real r = ::exp(c.x);
			if (c.y == 0)
				return complex(r);

			real a = c.y;
			ret.from_ra(r, a);
			return ret;
		}
	};
	const complex complex::ZERO = complex(0, 0);
	const complex complex::ONE = complex(1, 1);
	const complex complex::UX = complex(1, 0);
	const complex complex::UY = complex(0, 1);
	const complex complex::CENTER = complex(0.5, 0.5);

	// 黎曼函数
	inline complex riemann_function(const complex& s, 
		int n = 18 // n -> nan
		)
	{
		complex ret = complex::ZERO;
		for (int i = 1; i <= n; i++)
		{
			ret += 1.0 / (real(i) ^ (s));
		}
		return ret;
	}
	void rimann_function_vis()
	{
		for (int i = -100; i < 100; i++)
		{
			int cor = RNDCOR;
			for (int j = 0; j < 10000; j++)
			{
				number_math::complex s(0.5 + i / 50.0f, (j - 5000) / 500.0f);
				number_math::complex fv = riemann_function(s);
				//::pixel(vec2(fv.x, fv.y) / 40.0f, cor);
			}
		}

		for (int i = -100; i < 100; i++)
		{
			for (int j = -10000; j < 10000; j++)
			{
				number_math::complex s(0.5 - i / 10.0f, (j - 5000.0f) / 100.0f);
				number_math::complex fv = riemann_function(s);
				float r = sqrt(fv.x * fv.x + fv.y * fv.y);
				//pixel(vec2(s.x, s.y) / 100.0f, blendcor(0xFF0000FF, 0xFFFFFF00, r / 0.5f));
			}
		}
		for (int i = -100; i < 100; i++)
		{
			for (int j = -10000; j < 10000; j++)
			{
				number_math::complex s(0.5 + i / 10.0f, (j - 5000.0f) / 100.0f);
				number_math::complex fv = riemann_function(s);
				float r = sqrt(fv.x * fv.x + fv.y * fv.y);
				//pixel(vec2(s.x, s.y) / 100.0f, blendcor(0xFF0000FF, 0xFFFF0000, r / 0.5f));
			}
		}
	}

	// L函数
	real Lfunction(real s, std::function <real(int)> xn)
	{
		int n = 100; // -> nan
		real sum = 0;
		for (int i = 1; i < n; i++)
		{
			sum += xn(i) / pow(real(i), s);
		}
		return sum;
	}
	// ************************************************************************
	// 一般数定义 
	// ************************************************************************
	struct number
	{
		inline static const double GR1 = (sqrt(5.0) - 1.0) / 2.0;		// 黄金比例1
		inline static const double GR2 = (sqrt(5.0) + 1.0) / 2.0;		// 黄金比例2
		inline static const double E = 2.718281828;						// 自然对数e

		real power = 1;
		quaternion q;	// 归一化

		real val()
		{
			return ::exp(power);
		}
		void setpower(real val, real min = 0, real max = 1)
		{
			val = (val - min) / (max - min);
			if (val == 0)
				power = 0;
			else
				power = log(power);
		}
		number operator + (const number& a) const
		{
			number ra;
			ra.power = power;
			ra.q = q + a.q;
			return ra;
		}
		number operator * (const number& a) const
		{
			number ra;
			ra.power = power + a.power;
			ra.q = q * a.q;
			return ra;
		}
	};
}
