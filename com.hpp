/*******************************************************************************************************
*                  Common Header for Coordinate System Computing
*             Foundation definitions designed for coordinate system theory and computation
*               Contains necessary type aliases, constants and utility functions
*
*  Author: Guojun Pan
*  DOI: doi.org/10.5281/zenodo.14435613 (Concept DOI - all versions)
*
*  This implementation is based on the Computable Coordinate System theory,
*  which treats coordinate systems as first-class algebraic objects.
*
*  Priority Protection: This code is part of the theoretical framework published under
*  DOI: doi.org/10.5281/zenodo.14435613
********************************************************************************************************/
#pragma once

// Required system headers
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>

// Data type definitions
typedef float real;

// Device callable marker (CUDA compatible)
#ifndef DEVICE_CALLABLE
#define DEVICE_CALLABLE
#endif

// Mathematical constants
#ifndef PI
#define PI 3.14159265358979323846264338327950288f
#endif
#define TWOPI (2 * PI)
#define EPSILON 1e-5f
#define ISZERO(x) (fabs(x) < EPSILON)

// Debug and error handling macros
#define ASSERT(x) { if(!(x)) { std::stringstream ss; ss << "ASSERT FAILED! " << __FILE__ << "(" << __LINE__ << "): " << #x; throw std::runtime_error(ss.str()); } }
#define PRINT(msg) { std::cerr << msg << std::endl; }
#define PRINTVEC3(v) PRINT(#v << " = " << (v).x << "," << (v).y << "," << (v).z)

// Mathematical utility functions
inline real _MIN(real a, real b) { return ((a) < (b) ? (a) : (b)); }
inline real _MAX(real a, real b) { return ((a) > (b) ? (a) : (b)); }

// Angle/radian conversion
inline real radians(real deg) { return deg / 180.0f * PI; }
inline real degrees(real rad) { return rad * (180.0f / PI); }

// Numerical processing
inline bool is_zero(real x, real tor = EPSILON) { return fabs(x) < tor; }
DEVICE_CALLABLE inline bool check_equal(real a, real b, const real eps = EPSILON)
{
    return (fabs(a - b) <= eps);
}

// Forward declarations (avoid circular dependencies)
struct vector2;
struct vector3;
struct quaternion;
struct ucoord3;
struct vcoord3;
struct coord3;

// Type aliases for coordinate system calculations
typedef vector3 vec3;
typedef const vector3& crvec;
typedef vector3& rvec;

typedef vector2 vec2;
typedef const vector2& crvec2;
typedef vector2& rvec2;

typedef quaternion quat;
typedef const quaternion& crquat;

// Coordinate system aliases
typedef coord3 crd3;
typedef coord3& rcd3;
typedef const coord3& crcd3;

typedef vcoord3 vcd3;
typedef const vcoord3& crvcd3;

typedef ucoord3 ucd3;
typedef const ucoord3& crucd3;

// Frame terminology aliases
typedef coord3 frame3;

// Common vector constants (these are now functions)
#define ZERO3 vec3::ZERO()
#define UNITX vec3::UX()
#define UNITY vec3::UY()
#define UNITZ vec3::UZ()

// Frame field composite operator related function declarations
inline real lerp(real a, real b, real t) { return a * (1 - t) + b * t; }
inline real sign(real x) { return (x > 0) ? 1.0f : ((x < 0) ? -1.0f : 0.0f); }

// Hash function template (for unique identification of geometric objects)
template <typename T>
inline void hash_combine(std::size_t& seed, const T& v)
{
    seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// Specialized hash_combine for real type
inline void hash_combine(std::size_t& seed, const real& v)
{
    seed ^= std::hash<real>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// **********************************************************************
// 2D Vector
// Vector type definitions optimized for coordinate system theory
// DOI: doi.org/10.5281/zenodo.18217542
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
// Vector type definitions optimized for coordinate system theory
// DOI: doi.org/10.5281/zenodo.18217542
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

    // Constructors
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

    // Array access operators
    DEVICE_CALLABLE real& operator[](int ind)
    {
        return val[ind];
    }
    DEVICE_CALLABLE real operator[](int ind) const
    {
        return val[ind];
    }

    // Component extraction
    vector2 xy() const { return vector2(x, y); }
    vector2 xz() const { return vector2(x, z); }
    vector2 yz() const { return vector2(y, z); }

    // Arithmetic operators - Addition
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

    // Arithmetic operators - Subtraction
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

    // Arithmetic operators - Multiplication
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

    // Arithmetic operators - Division
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

    // Cross product operator (%)
    DEVICE_CALLABLE vector3 operator%(const vector3& v) const
    {
        return vector3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    // Comparison operators
    DEVICE_CALLABLE bool operator==(const vector3& rv) const
    {
        return (fabs(x - rv.x) <= sEPSILON && fabs(y - rv.y) <= sEPSILON && fabs(z - rv.z) <= sEPSILON);
    }
    DEVICE_CALLABLE bool operator!=(const vector3& rv) const
    {
        return (fabs(x - rv.x) > sEPSILON || fabs(y - rv.y) > sEPSILON || fabs(z - rv.z) > sEPSILON);
    }

    // Length and normalization methods
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

    // Dot product
    DEVICE_CALLABLE real dot(const vector3& v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }

    // Cross product
    DEVICE_CALLABLE vector3 cross(const vector3& v) const
    {
        return vector3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    // Cross-dot product (perpendicular component)
    DEVICE_CALLABLE vector3 crossdot(const vector3& n) const
    {
        const vector3& v = (*this);
        return v - n * v.dot(n);
    }

    // Project onto another vector
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

    // Linear interpolation
    static vector3 lerp(const vector3& v1, const vector3& v2, real t)
    {
        return v1 * (1 - t) + v2 * t;
    }

    // Angle between two vectors
    static real angle(const vector3& a, const vector3& b)
    {
        return acos(a.dot(b) / (a.len() * b.len()));
    }

    // Mean of components
    DEVICE_CALLABLE real mean() const
    {
        return (x + y + z) / 3.0f;
    }

    // Serialization
    std::string serialise() const
    {
        return std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z);
    }

    // Static constants
    static const vector3 ZERO;
    static const vector3 ONE;
    static const vector3 UX;
    static const vector3 UY;
    static const vector3 UZ;
};

// vector3 static constant initialization
real vector3::sEPSILON = EPSILON;
const vector3 vector3::ZERO = vector3(0, 0, 0);
const vector3 vector3::ONE = vector3(1, 1, 1);
const vector3 vector3::UX = vector3(1, 0, 0);
const vector3 vector3::UY = vector3(0, 1, 0);
const vector3 vector3::UZ = vector3(0, 0, 1);

// **********************************************************************
// Quaternion
// Quaternions are extensions of complex numbers, with unit quaternions
// used for rotation operations in coordinate system theory
// DOI: doi.org/10.5281/zenodo.18217542
// **********************************************************************
struct quaternion
{
    static const quaternion ONE;
    static const quaternion UX;
    static const quaternion UY;
    static const quaternion UZ;

    real w = 1, x = 0, y = 0, z = 0;

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
        return quaternion(fScalar * rkQ.w, fScalar * rkQ.x, fScalar * rkQ.y, fScalar * rkQ.z);
    }

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

    // Comparison operators
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
            return PI;
        real ang = acos(w) * 2;
        if (ang > PI)
            return ang - PI * 2;
        if (ang < -PI)
            return ang + PI * 2;
        return ang;
    }

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

    // Conjugation (inverse operation for coordinate system transformations)
    void conj()
    {
        this->x = -x; this->y = -y; this->z = -z;
    }
    quaternion conjcopy() const
    {
        return quaternion(w, -x, -y, -z);
    }

    // Inverse
    quaternion inverse() const {
        real lenSquared = w * w + x * x + y * y + z * z;
        if (lenSquared != 0) {
            real factor = 1.0f / lenSquared;
            return conjcopy() * factor;
        }
        return quaternion::ONE;
    }

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
            if (fabs(dot + 1.0f) < eps)    // Handle 180 degree case
            {
                vector3 ax;
                vector3 uz = vector3::UZ;
                if (fabs(v1.x) < eps && fabs(v1.y) < eps)
                    uz = -vector3::UX;
                ax = uz.cross(v1).normalized();
                ang_axis(PI, ax.normalized());
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
        real halfpi = 0.5f * PI;

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

    // Exponential operation (note operator precedence)
    quaternion operator^(real t)
    {
        return slerp(quaternion::ONE, *this, t);
    }

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

// Quaternion constant definitions
inline const quaternion quaternion::ONE = quaternion(1, 0, 0, 0);

