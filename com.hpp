/*********************************************************************
*                  Common Header for Coordinate System Computing
*             Foundation definitions designed for coordinate system theory and computation
*               Contains necessary type aliases, constants and utility functions
**********************************************************************/
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
struct coord2;

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

typedef coord2 crd2;
typedef const coord2& crcd2;

// Frame terminology aliases
typedef coord3 frame3;
typedef coord2 frame2;

// Common vector constants (these are now functions)
#define ZERO3 vec3::ZERO()
#define UNITX vec3::UX()
#define UNITY vec3::UY()
#define UNITZ vec3::UZ()

// Core mathematical components will be included elsewhere
// Avoid circular dependencies, forward declarations only

// Frame field composite operator related function declarations
inline real lerp(real a, real b, real t) { return a * (1 - t) + b * t; }
inline real sign(real x) { return (x > 0) ? 1.0f : ((x < 0) ? -1.0f : 0.0f); }

// Vector math functions (function definitions will be implemented after including vector.hpp)
vec3 floor(crvec v);

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