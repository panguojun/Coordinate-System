/**
*			【向量】
* 
*		向量的定义是来自于四元数
*		向量不是完整的数，
*		向量跟空间结构有关系，
*		如果在时空中建议使用四元数
*/
// **********************************************************************
// 2D
// **********************************************************************
struct vector2
{
    static const vector2 ZERO;
    static const vector2 ONE;
    static const vector2 UX;
    static const vector2 UY;
    static const vector2 CENTER;
    static const vector2 INFINITY2;

    union {
        real val[2];
        struct
        {
            real x;
            real y;
        };
    };

    real& operator[](int ind)
    {
        return val[ind];
    }
    DEVICE_CALLABLE vector2()
    {
        x = 0;
        y = 0;
    }
    DEVICE_CALLABLE vector2(const vector2& v)
    {
        x = v.x;
        y = v.y;
    }
    DEVICE_CALLABLE explicit vector2(real v)
    {
        x = v;
        y = v;
    }
    DEVICE_CALLABLE vector2(real _x, real _y)
    {
        x = _x;
        y = _y;
    }

    DEVICE_CALLABLE static vector2 ang_len(real _angle, real _r)
    {
        return vector2(_r * cos(_angle), _r * sin(_angle));
    }

    bool isINF() const
    {
        return x == INFINITY || y == INFINITY || x == -INFINITY || y == -INFINITY;
    }

    DEVICE_CALLABLE vector2 xx() const { return vector2(x, x); }
    DEVICE_CALLABLE vector2 xy() const { return vector2(x, y); }
    DEVICE_CALLABLE vector2 yx() const { return vector2(y, x); }
    DEVICE_CALLABLE vector2 yy() const { return vector2(y, y); }

    DEVICE_CALLABLE vector2 operator+(const vector2& _p) const
    {
        vector2 fp;
        fp.x = x + _p.x;
        fp.y = y + _p.y;

        return fp;
    }
    DEVICE_CALLABLE void operator+=(const vector2& _p)
    {
        x += _p.x;
        y += _p.y;
    }
    DEVICE_CALLABLE vector2 operator-(const vector2& _p) const
    {
        vector2 fp;
        fp.x = x - _p.x;
        fp.y = y - _p.y;
        return fp;
    }
    DEVICE_CALLABLE void operator-=(const vector2& _p)
    {
        x = x - _p.x;
        y = y - _p.y;
    }
    DEVICE_CALLABLE vector2 operator-() const
    {
        vector2 fp;
        fp.x = -x;
        fp.y = -y;
        return fp;
    }
    DEVICE_CALLABLE vector2 operator*(real s) const
    {
        vector2 fp;
        fp.x = s * x;
        fp.y = s * y;
        return fp;
    }
    DEVICE_CALLABLE void operator*=(real s)
    {
        x = s * x;
        y = s * y;
    }
    DEVICE_CALLABLE friend vector2 operator*(real s, const vector2& v)
    {
        vector2 fp;
        fp.x = v.x * s;
        fp.y = v.y * s;
        return fp;
    }
    DEVICE_CALLABLE vector2 operator*(const vector2& b) const
    {
        return vector2(x * b.x, y * b.y);
    }
    DEVICE_CALLABLE void operator*=(const vector2& b)
    {
        x = x * b.x;
        y = y * b.y;
    }
    DEVICE_CALLABLE vector2 operator/(real s) const
    {
        vector2 fp;
        fp.x = x / s;
        fp.y = y / s;
        return fp;
    }
    DEVICE_CALLABLE vector2 operator/(const vector2& b) const
    {
        vector2 fp;
        fp.x = x / b.x;
        fp.y = y / b.y;
        return fp;
    }
    DEVICE_CALLABLE void operator/=(const vector2& b)
    {
        x = x / b.x;
        y = y / b.y;
    }
    DEVICE_CALLABLE void operator/=(real s)
    {
        x = x / s;
        y = y / s;
    }
    DEVICE_CALLABLE bool operator==(const vector2& rv) const
    {
        return (fabs(x - rv.x) <= 1e-5 && fabs(y - rv.y) <= 1e-5);
    }
    DEVICE_CALLABLE bool operator!=(const vector2& rv) const
    {
        return (fabs(x - rv.x) > 1e-5 || fabs(y - rv.y) > 1e-5);
    }
    DEVICE_CALLABLE real len() const
    {
        return sqrt(x * x + y * y);
    }
    DEVICE_CALLABLE real length() const
    {
        return sqrt(x * x + y * y);
    }
    DEVICE_CALLABLE real sqrlen() const
    {
        return (x * x + y * y);
    }
    DEVICE_CALLABLE real angle() const
    {
        return atan2(y, x);
    }
    DEVICE_CALLABLE vector2 angle(real ang)
    {
        return ang_len(ang, 1);
    }
    DEVICE_CALLABLE void norm()
    {
        real r = len();
        if (r > 0)
        {
            x /= r;
            y /= r;
        }
    }
    DEVICE_CALLABLE void normalise()
    {
        real r = len();
        if (r > 0)
        {
            x /= r;
            y /= r;
        }
    }
    DEVICE_CALLABLE vector2 normcopy() const
    {
        real r = len();
        if (r > 0)
        {
            return vector2(x / r, y / r);
        }
        return vector2::ZERO;
    }
    DEVICE_CALLABLE vector2 normalized() const
    {
        real r = len();
        if (r > 0)
        {
            return vector2(x / r, y / r);
        }
        return vector2::ZERO;
    }
    DEVICE_CALLABLE void rot(real angle)
    {
        (*this) = complex_mul((*this), vector2::ang_len(angle, 1));
    }
    DEVICE_CALLABLE vector2 rotcopy(real angle) const
    {
        return complex_mul((*this), vector2::ang_len(angle, 1));
    }
    DEVICE_CALLABLE vector2 rotcpy(real angle) const
    {
        return complex_mul((*this), vector2::ang_len(angle, 1));
    }
    DEVICE_CALLABLE vector2 roted(real angle) const
    {
        return complex_mul((*this), vector2::ang_len(angle, 1));
    }
    DEVICE_CALLABLE void rot(real angle, const vector2& o)
    {
        vector2 v = (*this) - o;
        v = complex_mul(v, vector2::ang_len(angle, 1));
        (*this) = v + o;
    }
    DEVICE_CALLABLE vector2 rotcopy(real angle, const vector2& o) const
    {
        vector2 v = (*this) - o;
        v = complex_mul(v, vector2::ang_len(angle, 1));
        return v + o;
    }
    DEVICE_CALLABLE vector2 rotcpy(real angle, const vector2& o) const
    {
        vector2 v = (*this) - o;
        v = complex_mul(v, vector2::ang_len(angle, 1));
        return v + o;
    }
    DEVICE_CALLABLE vector2 roted(real angle, const vector2& o) const
    {
        vector2 v = (*this) - o;
        v = complex_mul(v, vector2::ang_len(angle, 1));
        return v + o;
    }
    DEVICE_CALLABLE friend vector2 complex_mul(const vector2& a, const vector2& b)
    {
        return vector2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
    }
    DEVICE_CALLABLE real dot(const vector2& v) const
    {
        return x * v.x + y * v.y;
    }
    DEVICE_CALLABLE real cross(const vector2& v) const
    {
        return x * v.y - y * v.x;
    }
    std::string serialise() const
    {
        return std::to_string(x) + "," + std::to_string(y);
    }
};
#if defined(PMDLL) || !defined(PM_IMPLEMENTED)
const vector2 vector2::ZERO = vector2(0, 0);
const vector2 vector2::ONE = vector2(1, 1);
const vector2 vector2::UX = vector2(1, 0);
const vector2 vector2::UY = vector2(0, 1);
const vector2 vector2::CENTER = vector2(0.5, 0.5);
const vector2 vector2::INFINITY2 = vector2(INFINITY, INFINITY);
#endif

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
    static const vector3 INFINITY3;
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
    DEVICE_CALLABLE vector3()
    {
        x = 0;
        y = 0;
        z = 0;
    }
    DEVICE_CALLABLE explicit vector3(real v)
    {
        x = v;
        y = v;
        z = v;
    }
    DEVICE_CALLABLE explicit vector3(real* v)
    {
        x = v[0];
        y = v[1];
        z = v[2];
    }
    DEVICE_CALLABLE vector3(real _x, real _y, real _z = 0.0f)
    {
        x = _x;
        y = _y;
        z = _z;
    }
    DEVICE_CALLABLE vector3(const vector3& v)
    {
        x = v.x;
        y = v.y;
        z = v.z;
    }
    DEVICE_CALLABLE real& operator[](int ind)
    {
        return val[ind];
    }
    DEVICE_CALLABLE real operator[](int ind) const
    {
        return val[ind];
    }

    vector2 xx() const { return vector2(x, x); }
    vector2 xy() const { return vector2(x, y); }
    vector2 xz() const { return vector2(x, z); }
    vector2 yx() const { return vector2(y, x); }
    vector2 yy() const { return vector2(y, y); }
    vector2 yz() const { return vector2(y, z); }
    vector2 zx() const { return vector2(z, x); }
    vector2 zy() const { return vector2(z, y); }
    vector2 zz() const { return vector2(z, z); }

    void xy(const vector2& _xy)
    {
        x = _xy.x;
        y = _xy.y;
    }
    void xz(const vector2& _xz)
    {
        x = _xz.x;
        z = _xz.y;
    }
    void yz(const vector2& _yz)
    {
        y = _yz.x;
        z = _yz.y;
    }

    vector3 xxx() const { return vector3(x, x, x); }
    vector3 xxy() const { return vector3(x, x, y); }
    vector3 xxz() const { return vector3(x, x, z); }
    vector3 xyx() const { return vector3(x, y, x); }
    vector3 xyy() const { return vector3(x, y, y); }
    vector3 xyz() const { return vector3(x, y, z); }
    vector3 xzx() const { return vector3(x, z, x); }
    vector3 xzy() const { return vector3(x, z, y); }
    vector3 xzz() const { return vector3(x, z, z); }
    vector3 yxx() const { return vector3(y, x, x); }
    vector3 yxy() const { return vector3(y, x, y); }
    vector3 yxz() const { return vector3(y, x, z); }
    vector3 yyx() const { return vector3(y, y, x); }
    vector3 yyy() const { return vector3(y, y, y); }
    vector3 yyz() const { return vector3(y, y, z); }
    vector3 yzx() const { return vector3(y, z, x); }
    vector3 yzy() const { return vector3(y, z, y); }
    vector3 yzz() const { return vector3(y, z, z); }
    vector3 zxx() const { return vector3(z, x, x); }
    vector3 zxy() const { return vector3(z, x, y); }
    vector3 zxz() const { return vector3(z, x, z); }
    vector3 zyx() const { return vector3(z, y, x); }
    vector3 zyy() const { return vector3(z, y, y); }
    vector3 zyz() const { return vector3(z, y, z); }
    vector3 zzx() const { return vector3(z, z, x); }
    vector3 zzy() const { return vector3(z, z, y); }
    vector3 zzz() const { return vector3(z, z, z); }
    vector3 xyo() const { return vector3(x, y, 0); }
    vector3 xoz() const { return vector3(x, 0, z); }
    vector3 oyz() const { return vector3(0, y, z); }

    DEVICE_CALLABLE vector3 operator+(const vector3& _p) const
    {
        vector3 fp;
        fp.x = x + _p.x;
        fp.y = y + _p.y;
        fp.z = z + _p.z;
        return fp;
    }
    DEVICE_CALLABLE void operator+=(const vector3& _p)
    {
        x += _p.x;
        y += _p.y;
        z += _p.z;
    }
    DEVICE_CALLABLE vector3 operator-(const vector3& _p) const
    {
        vector3 fp;
        fp.x = x - _p.x;
        fp.y = y - _p.y;
        fp.z = z - _p.z;
        return fp;
    }
    DEVICE_CALLABLE void operator-=(const vector3& _p)
    {
        x -= _p.x;
        y -= _p.y;
        z -= _p.z;
    }
    DEVICE_CALLABLE vector3 operator-() const
    {
        vector3 fp;
        fp.x = -x;
        fp.y = -y;
        fp.z = -z;
        return fp;
    }
    DEVICE_CALLABLE vector3 operator*(real s) const
    {
        vector3 fp;
        fp.x = s * x;
        fp.y = s * y;
        fp.z = s * z;
        return fp;
    }
    DEVICE_CALLABLE vector3 operator*(const vector3& v) const
    {
        vector3 fp;
        fp.x = v.x * x;
        fp.y = v.y * y;
        fp.z = v.z * z;
        return fp;
    }
    DEVICE_CALLABLE void operator*=(const vector3& s)
    {
        x = s.x * x;
        y = s.y * y;
        z = s.z * z;
    }
    DEVICE_CALLABLE friend vector3 operator*(real s, const vector3& v)
    {
        vector3 fp;
        fp.x = v.x * s;
        fp.y = v.y * s;
        fp.z = v.z * s;
        return fp;
    }
    DEVICE_CALLABLE void operator*=(real s)
    {
        x = s * x;
        y = s * y;
        z = s * z;
    }
    DEVICE_CALLABLE vector3 operator/(real s) const
    {
        vector3 fp;
        fp.x = x / s;
        fp.y = y / s;
        fp.z = z / s;
        return fp;
    }
    DEVICE_CALLABLE void operator/=(real s)
    {
        x = x / s;
        y = y / s;
        z = z / s;
    }
    DEVICE_CALLABLE vector3 operator/(const vector3& v) const
    {
        vector3 fp;
        fp.x = x / v.x;
        fp.y = y / v.y;
        fp.z = z / v.z;
        return fp;
    }
    DEVICE_CALLABLE void operator/=(const vector3& s)
    {
        x = x / s.x;
        y = y / s.y;
        z = z / s.z;
    }
    DEVICE_CALLABLE vector3 operator%(const vector3& v) const
    {
        return vector3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    DEVICE_CALLABLE bool operator==(const vector3& rv) const
    {
        return (fabs(x - rv.x) <= sEPSINON && fabs(y - rv.y) <= sEPSINON && fabs(z - rv.z) <= sEPSINON);
    }
    DEVICE_CALLABLE bool operator!=(const vector3& rv) const
    {
        return (fabs(x - rv.x) > sEPSINON || fabs(y - rv.y) > sEPSINON || fabs(z - rv.z) > sEPSINON);
    }
    DEVICE_CALLABLE bool operator<(const vector3& rv) const
    {
        if (x < rv.x && y < rv.y && z < rv.z)
            return true;
        return false;
    }
    DEVICE_CALLABLE bool operator<=(const vector3& rv) const
    {
        if (x <= rv.x && y <= rv.y && z <= rv.z)
            return true;
        return false;
    }
    DEVICE_CALLABLE bool operator>(const vector3& rv) const
    {
        if (x > rv.x && y > rv.y && z > rv.z)
            return true;
        return false;
    }
    DEVICE_CALLABLE bool operator>=(const vector3& rv) const
    {
        if (x >= rv.x && y >= rv.y && z >= rv.z)
            return true;
        return false;
    }
    DEVICE_CALLABLE bool isINF() const
    {
        return x == INFINITY || y == INFINITY || z == INFINITY ||
               x == -INFINITY || y == -INFINITY || z == -INFINITY;
    }
    DEVICE_CALLABLE vector3 flipX() const
    {
        return vector3(-x, y, z);
    }
    DEVICE_CALLABLE vector3 flipY() const
    {
        return vector3(x, -y, z);
    }
    DEVICE_CALLABLE vector3 flipZ() const
    {
        return vector3(x, y, -z);
    }
    DEVICE_CALLABLE real len() const
    {
        return sqrt(x * x + y * y + z * z);
    }
    DEVICE_CALLABLE real length() const
    {
        return sqrt(x * x + y * y + z * z);
    }
    DEVICE_CALLABLE real lenxy() const
    {
        return sqrt(x * x + y * y);
    }
    DEVICE_CALLABLE real sqrlenxy() const
    {
        return (x * x + y * y);
    }
    DEVICE_CALLABLE real lenxz() const
    {
        return sqrt(x * x + z * z);
    }
    DEVICE_CALLABLE real sqrlenxz() const
    {
        return (x * x + z * z);
    }
    DEVICE_CALLABLE real sqrlen() const
    {
        return (x * x + y * y + z * z);
    }
    DEVICE_CALLABLE real abslen() const
    {
        return abs(x) + abs(y) + abs(z);
    }
    DEVICE_CALLABLE bool norm()
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
    DEVICE_CALLABLE bool normalize()
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
    DEVICE_CALLABLE vector3 normcopy() const
    {
        real r = len();
        if (r > 0)
        {
            return vector3(this->x / r,
                           this->y / r,
                           this->z / r);
        }
        return vector3(0, 0, 0);
    }
    DEVICE_CALLABLE vector3 normalized() const
    {
        real r = len();
        if (r > 0)
            return vector3(this->x / r,
                           this->y / r,
                           this->z / r);
        return vector3(0, 0, 0);
    }
    DEVICE_CALLABLE real dot(const vector3& v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }
    DEVICE_CALLABLE vector3 crossdot(const vector3& n) const
    {
        const vector3& v = (*this);
        return v - n * v.dot(n);
    }
    DEVICE_CALLABLE vector3 cross(const vector3& v) const
    {
        vector3 n;
        // 这里使用了左手顺序！
        n.x = -(y * v.z - z * v.y);
        n.y = -(z * v.x - x * v.z);
        n.z = -(x * v.y - y * v.x);
        return n;
    }
    DEVICE_CALLABLE vector3 cross_left(const vector3& v) const
    {
        vector3 n;
        // 这里使用了左手顺序！
        n.x = -(y * v.z - z * v.y);
        n.y = -(z * v.x - x * v.z);
        n.z = -(x * v.y - y * v.x);
        return n;
    }
    DEVICE_CALLABLE vector3 cross_right(const vector3& v) const
    {
        vector3 n;
        n.x = (y * v.z - z * v.y);
        n.y = (z * v.x - x * v.z);
        n.z = (x * v.y - y * v.x);
        return n;
    }
    std::string serialise() const
    {
        return std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z);
    }

    // GUID
    std::size_t hash() const
    {
        std::size_t hash = 0;
        if (sizeof(real) == sizeof(float))
        {
            // 哈希原点坐标
            hash_combine(hash, x);
            hash_combine(hash, y);
            hash_combine(hash, z);
        }
        return hash;
    }
    static vector3 min3(const vector3& a, const vector3& b)
    {
        return vector3(_MIN(a.x, b.x), _MIN(a.y, b.y), _MIN(a.z, b.z));
    }
    static vector3 max3(const vector3& a, const vector3& b)
    {
        return vector3(_MAX(a.x, b.x), _MAX(a.y, b.y), _MAX(a.z, b.z));
    }
    static vector3 rnd(real min = 0, real max = 1)
    {
        return vector3(rrnd(min, max), rrnd(min, max), rrnd(min, max));
    }
    static vector3 rndrad(real r = 1)
    {
        return rnd(-1, 1).normcopy() * r;
    }
    static vector3 lerp(const vector3& v1, const vector3& v2, real t)
    {
        return (v1 * (1 - t) + v2);
    }
    static vector3 lerp(const vector3& v1, const vector3& v2, const vector3& t)
    {
        return (v1 * (vector3::ONE - t) + v2);
    }
    static real angle(const vector3& a, const vector3& b)
    {
        return std::acos(a.dot(b) / (a.len() * b.len()));
    }
    static real angle(const vector3& a, const vector3& b, const vector3& axis)
    {
        real angle = std::acos(a.dot(b) / (a.len() * b.len()));
        vector3 cross_dir = a.cross(b);
        if (cross_dir.dot(axis) < 0)
        {
            angle = 2 * PI - angle;
        }
        return angle;
    }
};

#if defined(PMDLL) || !defined(PM_IMPLEMENTED)
const vector3 vector3::ZERO = vector3(0, 0, 0);
const vector3 vector3::ONE = vector3(1, 1, 1);
const vector3 vector3::UX = vector3(1, 0, 0);
const vector3 vector3::UY = vector3(0, 1, 0);
const vector3 vector3::UZ = vector3(0, 0, 1);
const vector3 vector3::CENTER = vector3(0, 0, 0);
const vector3 vector3::INFINITY3 = vector3(INFINITY, INFINITY, INFINITY);
real vector3::sEPSINON = EPSINON;
#endif

// **********************************************************************
// 4D vector
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
    DEVICE_CALLABLE vector4()
    {
        x = 0;
        y = 0;
        z = 0;
        w = 0;
    }
    DEVICE_CALLABLE vector4(float _x, float _y, float _z, float _w)
    {
        x = _x;
        y = _y;
        z = _z;
        w = _w;
    }
    DEVICE_CALLABLE vector4(float _x, crvec _v3)
    {
        x = _x;
        y = _v3.x;
        z = _v3.y;
        w = _v3.z;
    }
    DEVICE_CALLABLE vector4(crvec _v3, float _w = 0)
    {
        x = _v3.x;
        y = _v3.y;
        z = _v3.z;
        w = _w;
    }
    DEVICE_CALLABLE explicit vector4(float _v)
    {
        x = _v;
        y = _v;
        z = _v;
        w = _v;
    }
    DEVICE_CALLABLE real& operator[](int ind)
    {
        return val[ind];
    }
    DEVICE_CALLABLE vector3 xyz() const { return vector3(x, y, z); }
    DEVICE_CALLABLE vector3 yzw() const { return vector3(y, z, w); }

    DEVICE_CALLABLE vector4 operator+(const vector4& _p) const
    {
        vector4 fp;
        fp.x = x + _p.x;
        fp.y = y + _p.y;
        fp.z = z + _p.z;
        fp.w = w + _p.w;
        return fp;
    }
    DEVICE_CALLABLE void operator+=(const vector4& _p)
    {
        x += _p.x;
        y += _p.y;
        z += _p.z;
        w += _p.w;
    }
    DEVICE_CALLABLE vector4 operator-(const vector4& _p) const
    {
        vector4 fp;
        fp.x = x - _p.x;
        fp.y = y - _p.y;
        fp.z = z - _p.z;
        fp.w = w - _p.w;
        return fp;
    }
    DEVICE_CALLABLE void operator-=(const vector4& _p)
    {
        x -= _p.x;
        y -= _p.y;
        z -= _p.z;
        w -= _p.w;
    }
    DEVICE_CALLABLE vector4 operator-() const
    {
        vector4 fp;
        fp.x = -x;
        fp.y = -y;
        fp.z = -z;
        fp.w = -w;
        return fp;
    }
    DEVICE_CALLABLE vector4 operator*(float s) const
    {
        vector4 fp;
        fp.x = s * x;
        fp.y = s * y;
        fp.z = s * z;
        fp.w = s * w;
        return fp;
    }
    DEVICE_CALLABLE friend vector4 operator*(real s, const vector4& v)
    {
        vector4 fp;
        fp.x = v.x * s;
        fp.y = v.y * s;
        fp.z = v.z * s;
        fp.w = v.w * s;
        return fp;
    }
    DEVICE_CALLABLE void operator*=(float s)
    {
        x = s * x;
        y = s * y;
        z = s * z;
        w = s * w;
    }
    DEVICE_CALLABLE vector4 operator/(float s) const
    {
        vector4 fp;
        fp.x = x / s;
        fp.y = y / s;
        fp.z = z / s;
        fp.w = w / s;
        return fp;
    }
    DEVICE_CALLABLE void operator/=(float s)
    {
        x = x / s;
        y = y / s;
        z = z / s;
        w = w / s;
    }
    DEVICE_CALLABLE bool operator==(const vector4& rv) const
    {
        return (fabs(x - rv.x) <= EPSINON && fabs(y - rv.y) <= EPSINON && fabs(z - rv.z) <= EPSINON && fabs(w - rv.w) <= EPSINON);
    }
    DEVICE_CALLABLE bool operator!=(const vector4& rv) const
    {
        return (fabs(x - rv.x) > EPSINON || fabs(y - rv.y) > EPSINON || fabs(z - rv.z) > EPSINON || fabs(w - rv.w) > EPSINON);
    }
    DEVICE_CALLABLE float len() const
    {
        return sqrt(x * x + y * y + z * z + w * w);
    }
    DEVICE_CALLABLE float sqrlen() const
    {
        return (x * x + y * y + z * z + w * w);
    }
    DEVICE_CALLABLE real norm()
    {
        float r = len();
        if (r > 1e-5)
        {
            x /= r;
            y /= r;
            z /= r;
            w /= r;
        }
        return r;
    }
    DEVICE_CALLABLE vector4 normcopy()
    {
        float r = len();
        if (r > 0)
        {
            return vector4(
                this->x / r,
                this->y / r,
                this->z / r,
                this->w / r);
        }
        return vector4(0, 0, 0, 0);
    }
    DEVICE_CALLABLE vector4 normalized()
    {
        float r = len();
        if (r > 0)
        {
            return vector4(
                this->x / r,
                this->y / r,
                this->z / r,
                this->w / r);
        }
        return vector4(0, 0, 0, 0);
    }

    DEVICE_CALLABLE float dot(const vector4& v) const
    {
        return x * v.x + y * v.y + z * v.z + w * v.w;
    }

    // Cross4 computes the four-dimensional cross product of the three vectors
    // U, V and W, in that order. It returns the resulting four-vector.

    DEVICE_CALLABLE vector4 cross(const vector4& V, const vector4& W) const
    {
        vector4 result;
        double A, B, C, D, E, F; // Intermediate Values

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
#if defined(PMDLL) || !defined(PM_IMPLEMENTED)
const vector4 vector4::ZERO = vector4(0, 0, 0, 0);
const vector4 vector4::UX = vector4(1, 0, 0, 0);
const vector4 vector4::UY = vector4(0, 1, 0, 0);
const vector4 vector4::UZ = vector4(0, 0, 1, 0);
const vector4 vector4::UW = vector4(0, 0, 0, 1);
const vector4 vector4::CENTER = vector4(0.5, 0.5, 0.5, 0.5);
#endif

// **********************************************************************
//							【nD Vector】
// 多维特征空间
// 来自于人类独特的思维方式（比如人鱼，半人马等）
// 把现象分解成多个维度的属性组合，然后再分别不同维度上进行融合，
// 最后可以合成出某种新的想象，如果在“本宇宙”得到验证，那么加强其权重，
// 否则削弱其权重。
// 注意，在物理学中（粒子物理学）至今没有观察到超过（3+1）自由度的高维粒子
// 所以说上层的高维的属性（比如尺寸与颜色）必然是某种非独立变量的某种权重组合。
// **********************************************************************
struct vectorn
{
    std::vector<real> val;

    static const vectorn ZERO;
    static const vectorn ONE;
    static const vectorn CENTER;

    real& operator[](int ind)
    {
        if (ind >= val.size())
        { // 自动扩容
            val.resize(ind + 1);
        }
        return val[ind];
    }
    real operator[](int ind) const
    {
        if (ind >= val.size())
        {
            return 0;
        }
        return val[ind];
    }
    vectorn() {}
    vectorn(const std::initializer_list<real>& list)
        : val(list) {}
    explicit vectorn(real v, int size)
        : val(size, v) {}
    explicit vectorn(real v)
    {
        for (int i = 0; i < val.size(); i++)
        {
            val[i] = v;
        }
    }
    vectorn& operator<<(float v)
    {
        val.push_back(v);
        return *this;
    }
    vectorn& operator>>(float& v)
    {
        v = val.back();
        val.pop_back();
        return *this;
    }
    vectorn& operator<<(double v)
    {
        val.push_back(v);
        return *this;
    }
    vectorn& operator>>(double& v)
    {
        v = val.back();
        val.pop_back();
        return *this;
    }
    vectorn& operator<<(const vector3& v)
    {
        val.push_back(v.x);
        val.push_back(v.y);
        val.push_back(v.z);
        return *this;
    }
    vectorn& operator>>(vector3& v)
    {
        v.x = val.back();
        val.pop_back();
        v.z = val.back();
        val.pop_back();
        v.z = val.back();
        val.pop_back();
        return *this;
    }

    void operator=(const vector3& v)
    {
        val[0] = v.x;
        val[1] = v.y;
        val[2] = v.z;
    }
    vectorn operator+(const vectorn& _p) const
    {
        vectorn fp;
        for (int i = 0; i < val.size(); i++)
        {
            fp[i] = val[i] + _p[i];
        }
        return fp;
    }
    void operator+=(const vectorn& _p)
    {
        for (int i = 0; i < val.size(); i++)
        {
            val[i] += _p[i];
        }
    }
    vectorn operator-(const vectorn& _p) const
    {
        vectorn fp;
        for (int i = 0; i < val.size(); i++)
        {
            fp[i] = val[i] - _p[i];
        }
        return fp;
    }
    void operator-=(const vectorn& _p)
    {
        for (int i = 0; i < val.size(); i++)
        {
            val[i] -= _p[i];
        }
    }
    vectorn operator-() const
    {
        vectorn fp;
        for (int i = 0; i < val.size(); i++)
        {
            fp[i] = -val[i];
        }
        return fp;
    }
    vectorn operator*(const vectorn& _p) const
    {
        vectorn fp;
        for (int i = 0; i < val.size(); i++)
        {
            fp[i] = val[i] * _p[i];
        }
        return fp;
    }
    void operator*=(const vectorn& _p)
    {
        for (int i = 0; i < val.size(); i++)
        {
            val[i] *= _p[i];
        }
    }
    vectorn operator*(real s) const
    {
        vectorn fp;
        for (int i = 0; i < val.size(); i++)
        {
            fp[i] = val[i] * s;
        }
        return fp;
    }
    void operator*=(real s)
    {
        for (int i = 0; i < val.size(); i++)
        {
            val[i] *= s;
        }
    }
    friend vectorn operator*(real s, vectorn& v)
    {
        vectorn fp;
        for (int i = 0; i < v.val.size(); i++)
        {
            fp[i] = v[i] * s;
        }
        return fp;
    }
    vectorn operator/(real s) const
    {
        vectorn fp;
        for (int i = 0; i < val.size(); i++)
        {
            fp[i] = val[i] / s;
        }
        return fp;
    }
    void operator/=(real s)
    {
        for (int i = 0; i < val.size(); i++)
        {
            val[i] /= s;
        }
    }
    bool operator==(const vectorn& v) const
    {
        if (val.size() != v.val.size())
            return false;
        bool ret = true;
        for (int i = 0; i < val.size(); i++)
        {
            ret &= (fabs(val[i] - v.val[i]) <= 1e-5);
            if (!ret)
                break;
        }
        return ret;
    }
    bool operator!=(const vectorn& v) const
    {
        if (val.size() != v.val.size())
            return true;
        bool ret = false;
        for (int i = 0; i < val.size(); i++)
        {
            ret |= (fabs(val[i] - v.val[i]) <= 1e-5);
            if (ret)
                break;
        }
        return ret;
    }
    int dim() const
    {
        return val.size();
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
        // ASSERT(val.size() == v.val.size());
        real sum = 0;
        for (int i = 0; i < val.size(); i++)
        {
            sum += val[i] * v.val[i];
        }
        return sum;
    }
    size_t hash() const
    {
        size_t xorResult = 0;
        for (const auto& value : val)
        {
            xorResult ^= static_cast<size_t>(value);
        }
        return xorResult;
    }
};
// **********************************************************************
#if defined(PMDLL) || !defined(PM_IMPLEMENTED)
const vectorn vectorn::ZERO = vectorn(0);
const vectorn vectorn::ONE = vectorn(1);
const vectorn vectorn::CENTER = vectorn(0.5);
#endif
