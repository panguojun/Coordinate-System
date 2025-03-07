/**
 *                          【Vector】
 *
 *                  The definition of a vector comes from quaternions.
 *                  A vector is not a complete number,
 *                  and it is related to the structure of space.
 *                  If you suggest in spacetime, consider using quaternions.
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
    static real sEPSILON;

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
    
    vector2()
    {
        x = 0;
        y = 0;
    }
    
    vector2(const vector2& v)
    {
        x = v.x;
        y = v.y;
    }
    
    explicit vector2(real v)
    {
        x = v;
        y = v;
    }
    
    vector2(real _x, real _y)
    {
        x = _x;
        y = _y;
    }

    static vector2 ang_len(real _angle, real _r)
    {
        return vector2(_r * cos(_angle), _r * sin(_angle));
    }

    bool isINF() const
    {
        return x == INFINITY || y == INFINITY || x == -INFINITY || y == -INFINITY;
    }

    vector2 xx() const { return vector2(x, x); }
    vector2 xy() const { return vector2(x, y); }
    vector2 yx() const { return vector2(y, x); }
    vector2 yy() const { return vector2(y, y); }

    vector2 operator+(const vector2& _p) const
    {
        vector2 fp;
        fp.x = x + _p.x;
        fp.y = y + _p.y;
        return fp;
    }
    
    void operator+=(const vector2& _p)
    {
        x += _p.x;
        y += _p.y;
    }
    
    vector2 operator-(const vector2& _p) const
    {
        vector2 fp;
        fp.x = x - _p.x;
        fp.y = y - _p.y;
        return fp;
    }
    
    void operator-=(const vector2& _p)
    {
        x -= _p.x;
        y -= _p.y;
    }
    
    vector2 operator-() const
    {
        return vector2(-x, -y);
    }

    bool operator<(const vector2& rv) const
    {
        return (x < rv.x && y < rv.y);
    }
    
    bool operator<=(const vector2& rv) const
    {
        return (x <= rv.x && y <= rv.y);
    }
    
    bool operator>(const vector2& rv) const
    {
        return (x > rv.x && y > rv.y);
    }
    
    bool operator>=(const vector2& rv) const
    {
        return (x >= rv.x && y >= rv.y);
    }
    
    vector2 operator*(real s) const
    {
        return vector2(s * x, s * y);
    }
    
    void operator*=(real s)
    {
        x *= s;
        y *= s;
    }
    
    friend vector2 operator*(real s, const vector2& v)
    {
        return vector2(v.x * s, v.y * s);
    }
    
    vector2 operator*(const vector2& b) const
    {
        return vector2(x * b.x, y * b.y);
    }
    
    void operator*=(const vector2& b)
    {
        x *= b.x;
        y *= b.y;
    }
    
    vector2 operator/(real s) const
    {
        return vector2(x / s, y / s);
    }
    
    vector2 operator/(const vector2& b) const
    {
        return vector2(x / b.x, y / b.y);
    }
    
    void operator/=(const vector2& b)
    {
        x /= b.x;
        y /= b.y;
    }
    
    void operator/=(real s)
    {
        x /= s;
        y /= s;
    }
    
    bool operator==(const vector2& rv) const
    {
        return (fabs(x - rv.x) <= sEPSILON && fabs(y - rv.y) <= sEPSILON);
    }
    
    bool operator!=(const vector2& rv) const
    {
        return (fabs(x - rv.x) > sEPSILON || fabs(y - rv.y) > sEPSILON);
    }
    
    real len() const
    {
        return sqrt(x * x + y * y);
    }
    
    real length() const
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
    
    vector2 angle(real ang)
    {
        return ang_len(ang, 1);
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
    
    void normalize()
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
    
    vector2 normalized() const
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
        (*this) = complex_mul((*this), vector2::ang_len(angle, 1));
    }
    
    vector2 rotcopy(real angle) const
    {
        return complex_mul((*this), vector2::ang_len(angle, 1));
    }
    
    vector2 rotcpy(real angle) const
    {
        return complex_mul((*this), vector2::ang_len(angle, 1));
    }
    
    vector2 roted(real angle) const
    {
        return complex_mul((*this), vector2::ang_len(angle, 1));
    }
    
    void rot(real angle, const vector2& o)
    {
        vector2 v = (*this) - o;
        v = complex_mul(v, vector2::ang_len(angle, 1));
        (*this) = v + o;
    }
    
    vector2 rotcopy(real angle, const vector2& o) const
    {
        vector2 v = (*this) - o;
        v = complex_mul(v, vector2::ang_len(angle, 1));
        return v + o;
    }
    
    vector2 rotcpy(real angle, const vector2& o) const
    {
        vector2 v = (*this) - o;
        v = complex_mul(v, vector2::ang_len(angle, 1));
        return v + o;
    }
    
    vector2 roted(real angle, const vector2& o) const
    {
        vector2 v = (*this) - o;
        v = complex_mul(v, vector2::ang_len(angle, 1));
        return v + o;
    }
    
    friend vector2 complex_mul(const vector2& a, const vector2& b)
    {
        return vector2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
    }
    
    real dot(const vector2& v) const
    {
        return x * v.x + y * v.y;
    }
    
    real cross(const vector2& v) const
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
real vector2::sEPSILON = EPSILON;
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
    static real sEPSILON;

    union {
        real val[3];
        struct
        {
            real x;
            real y;
            real z;
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
    
    explicit vector3(real* v)
    {
        x = v[0];
        y = v[1];
        z = v[2];
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
    
    real& operator[](int ind)
    {
        return val[ind];
    }
    
    real operator[](int ind) const
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

    vector3 operator+(const vector3& _p) const
    {
        return vector3(x + _p.x, y + _p.y, z + _p.z);
    }
    
    vector3 operator+=(const vector3& _p)
    {
        x += _p.x;
        y += _p.y;
        z += _p.z;
        return *this;
    }
    
    vector3 operator-(const vector3& _p) const
    {
        return vector3(x - _p.x, y - _p.y, z - _p.z);
    }
    
    vector3 operator-=(const vector3& _p)
    {
        x -= _p.x;
        y -= _p.y;
        z -= _p.z;
        return *this;
    }
    
    vector3 operator-() const
    {
        return vector3(-x, -y, -z);
    }
    
    vector3 operator*(real s) const
    {
        return vector3(s * x, s * y, s * z);
    }
    
    vector3 operator*(const vector3& v) const
    {
        return vector3(v.x * x, v.y * y, v.z * z);
    }
    
    vector3 operator*=(const vector3& s)
    {
        x *= s.x;
        y *= s.y;
        z *= s.z;
        return *this;
    }
    
    friend vector3 operator*(real s, const vector3& v)
    {
        return vector3(v.x * s, v.y * s, v.z * s);
    }
    
    void operator*=(real s)
    {
        x *= s;
        y *= s;
        z *= s;
    }
    
    vector3 operator/(real s) const
    {
        return vector3(x / s, y / s, z / s);
    }
    
    vector3 operator/=(real s)
    {
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }
    
    vector3 operator/(const vector3& v) const
    {
        return vector3(x / v.x, y / v.y, z / v.z);
    }
    
    vector3 operator/=(const vector3& s)
    {
        x /= s.x;
        y /= s.y;
        z /= s.z;
        return *this;
    }

    vector3 operator%(const vector3& v) const
    {
        return vector3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    
    bool operator==(const vector3& rv) const
    {
        return (fabs(x - rv.x) <= sEPSILON && fabs(y - rv.y) <= sEPSILON && fabs(z - rv.z) <= sEPSILON);
    }
    
    bool operator!=(const vector3& rv) const
    {
        return (fabs(x - rv.x) > sEPSILON || fabs(y - rv.y) > sEPSILON || fabs(z - rv.z) > sEPSILON);
    }
    
    bool operator<(const vector3& rv) const
    {
        return (x < rv.x && y < rv.y && z < rv.z);
    }
    
    bool operator<=(const vector3& rv) const
    {
        return (x <= rv.x && y <= rv.y && z <= rv.z);
    }
    
    bool operator>(const vector3& rv) const
    {
        return (x > rv.x && y > rv.y && z > rv.z);
    }
    
    bool operator>=(const vector3& rv) const
    {
        return (x >= rv.x && y >= rv.y && z >= rv.z);
    }
    
    bool isINF() const
    {
        return x == INFINITY || y == INFINITY || z == INFINITY ||
               x == -INFINITY || y == -INFINITY || z == -INFINITY;
    }
    
    vector3 flipX() const
    {
        return vector3(-x, y, z);
    }
    
    vector3 flipY() const
    {
        return vector3(x, -y, z);
    }
    
    vector3 flipZ() const
    {
        return vector3(x, y, -z);
    }
    
    real len() const
    {
        return sqrt(x * x + y * y + z * z);
    }
    
    real length() const
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
    
    real sqrlenxz() const
    {
        return (x * x + z * z);
    }
    
    real sqrlen() const
    {
        return (x * x + y * y + z * z);
    }
    
    real abslen() const
    {
        return abs(x) + abs(y) + abs(z);
    }
    
    real volum() const
    {
        return abs(x) * abs(y) * abs(z);
    }
    
    bool norm()
    {
        real _r = len();
        if (_r > 1e-6)
        {
            x /= _r;
            y /= _r;
            z /= _r;
            return true;
        }
        return false;
    }
    
    bool normalize()
    {
        real _r = len();
        if (_r > 1e-6)
        {
            x /= _r;
            y /= _r;
            z /= _r;
            return true;
        }
        return false;
    }
    
    vector3 normcopy() const
    {
        real _r = len();
        if (_r > 0)
        {
            return vector3(this->x / _r,
                           this->y / _r,
                           this->z / _r);
        }
        return vector3(0, 0, 0);
    }
    
    vector3 normalized() const
    {
        real _r = len();
        if (_r > 0)
            return vector3(this->x / _r,
                           this->y / _r,
                           this->z / _r);
        return vector3(0, 0, 0);
    }
    
    real dot(const vector3& v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }
    
    vector3 crossdot(const vector3& n) const
    {
        const vector3& v = (*this);
        return v - n * v.dot(n);
    }
    
    vector3 cross(const vector3& v) const
    {
        vector3 n;
        // Here we use the left-hand rule!
        n.x = -(y * v.z - z * v.y);
        n.y = -(z * v.x - x * v.z);
        n.z = -(x * v.y - y * v.x);
        return n;
    }
    
    vector3 cross_left(const vector3& v) const
    {
        vector3 n;
        // Here we use the left-hand rule!
        n.x = -(y * v.z - z * v.y);
        n.y = -(z * v.x - x * v.z);
        n.z = -(x * v.y - y * v.x);
        return n;
    }
    
    vector3 cross_right(const vector3& v) const
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
    std::size_t hash(int precision_level) const
    {
        std::size_t hashValue = 0;

        if constexpr (sizeof(real) == sizeof(float))
        {
            // Keep precisionLevel decimal places
            real factor = std::pow(10.0f, precision_level);

            // Hash the coordinates
            hash_combine(hashValue, std::round(x * factor) / factor);
            hash_combine(hashValue, std::round(y * factor) / factor);
            hash_combine(hashValue, std::round(z * factor) / factor);
        }
        return hashValue;
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
// **********************************************************************
const vector3 vector3::ZERO = vector3(0, 0, 0);
const vector3 vector3::ONE = vector3(1, 1, 1);
const vector3 vector3::UX = vector3(1, 0, 0);
const vector3 vector3::UY = vector3(0, 1, 0);
const vector3 vector3::UZ = vector3(0, 0, 1);
const vector3 vector3::CENTER = vector3(0, 0, 0);
const vector3 vector3::INFINITY3 = vector3(INFINITY, INFINITY, INFINITY);
real vector3::sEPSILON = EPSILON;