/**
 * Coordinate System Computing Vector Library
 *
 * Vector type definitions optimized for coordinate system theory and computation
 * Provides basic operations for 2D and 3D vectors
 * Supports geometric calculations with frame field composite operators
 */
#pragma once

// Forward declarations
typedef float real;
#ifndef DEVICE_CALLABLE
#define DEVICE_CALLABLE
#endif
#ifndef EPSILON
#define EPSILON 1e-5f
#endif

// Required math functions
#include <cmath>
#include <string>

// **********************************************************************
// 2D Vector
// **********************************************************************
struct vector2
{
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

    DEVICE_CALLABLE vector2 operator+(const vector2& _p) const
    {
        return vector2(x + _p.x, y + _p.y);
    }
    DEVICE_CALLABLE void operator+=(const vector2& _p)
    {
        x += _p.x;
        y += _p.y;
    }
    DEVICE_CALLABLE vector2 operator-(const vector2& _p) const
    {
        return vector2(x - _p.x, y - _p.y);
    }
    DEVICE_CALLABLE void operator-=(const vector2& _p)
    {
        x -= _p.x;
        y -= _p.y;
    }
    DEVICE_CALLABLE vector2 operator-() const
    {
        return vector2(-x, -y);
    }

    DEVICE_CALLABLE vector2 operator*(real s) const
    {
        return vector2(s * x, s * y);
    }
    DEVICE_CALLABLE void operator*=(real s)
    {
        x = s * x;
        y = s * y;
    }
    DEVICE_CALLABLE friend vector2 operator*(real s, const vector2& v)
    {
        return vector2(v.x * s, v.y * s);
    }
    DEVICE_CALLABLE vector2 operator*(const vector2& b) const
    {
        return vector2(x * b.x, y * b.y);
    }
    DEVICE_CALLABLE vector2 operator/(real s) const
    {
        return vector2(x / s, y / s);
    }
    DEVICE_CALLABLE void operator/=(real s)
    {
        x = x / s;
        y = y / s;
    }
    DEVICE_CALLABLE bool operator==(const vector2& rv) const
    {
        return (fabs(x - rv.x) <= sEPSILON && fabs(y - rv.y) <= sEPSILON);
    }
    DEVICE_CALLABLE bool operator!=(const vector2& rv) const
    {
        return (fabs(x - rv.x) > sEPSILON || fabs(y - rv.y) > sEPSILON);
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
    DEVICE_CALLABLE void normalize()
    {
        real r = len();
        if (r > sEPSILON)
        {
            x /= r;
            y /= r;
        }
    }
    DEVICE_CALLABLE vector2 normalized() const
    {
        real r = len();
        if (r > sEPSILON)
        {
            return vector2(x / r, y / r);
        }
        return vector2(0, 0);
    }
    DEVICE_CALLABLE real dot(const vector2& v) const
    {
        return x * v.x + y * v.y;
    }
    DEVICE_CALLABLE real cross(const vector2& v) const
    {
        return x * v.y - y * v.x;
    }
    DEVICE_CALLABLE friend vector2 complex_mul(const vector2& a, const vector2& b)
    {
        return vector2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
    }
};

// vector2 static constant initialization
real vector2::sEPSILON = EPSILON;

// **********************************************************************
// 3D Vector
// **********************************************************************
struct vector3
{
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
    DEVICE_CALLABLE vector3(real _x, real _y, real _z = 0)
    {
        x = _x;
        y = _y;
        z = _z;
    }
    DEVICE_CALLABLE vector3(const vector2& v, real _z = 0)
    {
        x = v.x;
        y = v.y;
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

    vector2 xy() const { return vector2(x, y); }
    vector2 xz() const { return vector2(x, z); }
    vector2 yz() const { return vector2(y, z); }

    DEVICE_CALLABLE vector3 operator+(const vector3& _p) const
    {
        return vector3(x + _p.x, y + _p.y, z + _p.z);
    }
    DEVICE_CALLABLE vector3 operator+=(const vector3& _p)
    {
        x += _p.x;
        y += _p.y;
        z += _p.z;
        return *this;
    }
    DEVICE_CALLABLE vector3 operator-(const vector3& _p) const
    {
        return vector3(x - _p.x, y - _p.y, z - _p.z);
    }
    DEVICE_CALLABLE vector3 operator-=(const vector3& _p)
    {
        x -= _p.x;
        y -= _p.y;
        z -= _p.z;
        return *this;
    }
    DEVICE_CALLABLE vector3 operator-() const
    {
        return vector3(-x, -y, -z);
    }
    DEVICE_CALLABLE vector3 operator*(real s) const
    {
        return vector3(s * x, s * y, s * z);
    }
    DEVICE_CALLABLE vector3 operator*(const vector3& v) const
    {
        return vector3(v.x * x, v.y * y, v.z * z);
    }
    DEVICE_CALLABLE vector3 operator*=(const vector3& s)
    {
        x = s.x * x;
        y = s.y * y;
        z = s.z * z;
        return (*this);
    }
    DEVICE_CALLABLE friend vector3 operator*(real s, const vector3& v)
    {
        return vector3(v.x * s, v.y * s, v.z * s);
    }
    DEVICE_CALLABLE void operator*=(real s)
    {
        x = s * x;
        y = s * y;
        z = s * z;
    }
    DEVICE_CALLABLE vector3 operator/(real s) const
    {
        return vector3(x / s, y / s, z / s);
    }
    DEVICE_CALLABLE vector3 operator/=(real s)
    {
        x = x / s;
        y = y / s;
        z = z / s;
        return (*this);
    }
    DEVICE_CALLABLE vector3 operator/(const vector3& v) const
    {
        return vector3(x / v.x, y / v.y, z / v.z);
    }
    DEVICE_CALLABLE vector3 operator/=(const vector3& s)
    {
        x = x / s.x;
        y = y / s.y;
        z = z / s.z;
        return (*this);
    }
    DEVICE_CALLABLE vector3 operator%(const vector3& v) const
    {
        return vector3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    DEVICE_CALLABLE bool operator==(const vector3& rv) const
    {
        return (fabs(x - rv.x) <= sEPSILON && fabs(y - rv.y) <= sEPSILON && fabs(z - rv.z) <= sEPSILON);
    }
    DEVICE_CALLABLE bool operator!=(const vector3& rv) const
    {
        return (fabs(x - rv.x) > sEPSILON || fabs(y - rv.y) > sEPSILON || fabs(z - rv.z) > sEPSILON);
    }
    DEVICE_CALLABLE real len() const
    {
        return sqrt(x * x + y * y + z * z);
    }
    DEVICE_CALLABLE real length() const
    {
        return sqrt(x * x + y * y + z * z);
    }
    DEVICE_CALLABLE real sqrlen() const
    {
        return (x * x + y * y + z * z);
    }
    DEVICE_CALLABLE real len_squared() const
    {
        return (x * x + y * y + z * z);
    }
    // Normalization
    DEVICE_CALLABLE bool normalize()
    {
        real _r = len();
        if (_r > sEPSILON)
        {
            x /= _r;
            y /= _r;
            z /= _r;
            return true;
        }
        return false;
    }
    DEVICE_CALLABLE vector3 normcopy() const
    {
        real _r = len();
        if (_r > sEPSILON)
            return vector3(this->x / _r, this->y / _r, this->z / _r);
        return vector3(0, 0, 0);
    }
    DEVICE_CALLABLE vector3 normalized() const
    {
        real _r = len();
        if (_r > sEPSILON)
            return vector3(this->x / _r, this->y / _r, this->z / _r);
        return vector3(0, 0, 0);
    }
    DEVICE_CALLABLE real dot(const vector3& v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }
    DEVICE_CALLABLE vector3 cross(const vector3& v) const
    {
        return vector3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    DEVICE_CALLABLE vector3 crossdot(const vector3& n) const
    {
        const vector3& v = (*this);
        return v - n * v.dot(n);
    }
    DEVICE_CALLABLE vector3 project(const vector3& T) const
    {
        real dotProduct = dot(T);
        real TLengthSquared = T.dot(T);
        if (TLengthSquared > sEPSILON)
        {
            return (dotProduct / TLengthSquared) * T;
        }
        return vector3(0, 0, 0);
    }

    static vector3 lerp(const vector3& v1, const vector3& v2, real t)
    {
        return v1 * (1 - t) + v2 * t;
    }
    static real angle(const vector3& a, const vector3& b)
    {
        return acos(a.dot(b) / (a.len() * b.len()));
    }

    DEVICE_CALLABLE real mean() const
    {
        return (x + y + z) / 3.0f;
    }

    std::string serialise() const
    {
        return std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z);
    }

    // Static constant access functions
    static vector3 ZERO() { return vector3(0, 0, 0); }
    static vector3 ONE() { return vector3(1, 1, 1); }
    static vector3 UX() { return vector3(1, 0, 0); }
    static vector3 UY() { return vector3(0, 1, 0); }
    static vector3 UZ() { return vector3(0, 0, 1); }
};

// vector3 static constant initialization
real vector3::sEPSILON = EPSILON;

// Type aliases required for coordinate system calculations
typedef vector3 vec3;
typedef const vector3& crvec;
typedef vector3& rvec;
typedef vector2 vec2;
typedef const vector2& crvec2;