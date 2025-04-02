/*********************************************************************
*                       Common Header File
* 
**********************************************************************/
#if defined(WINDOWS) 
    #if defined(_MSC_VER) && (defined(PMDLL) || defined(DGE_CORE))
    #define EXPORT_API __declspec(dllexport)
    #else
    #define EXPORT_API __declspec(dllimport)
    #endif
#elif defined(__unix__) || defined(__linux__)
    #define EXPORT_API __attribute__((visibility("default")))
    #include <list>
#endif

// *********************************************************************
// Macros
// *********************************************************************
#ifdef WINDOWS

#define UNUSED(x)                    (void)(x)

#ifdef ECHO_OFF
    #define PRINT(msg)              {}
#else
#define PRINT(msg)                  {MYTRACE(msg);}
#endif
    #define COSOCOR(cor)            SetConsoleTextAttribute(hConsole, cor | FOREGROUND_INTENSITY); 
#else
    #define PRINT(msg)              { std::cerr << msg << std::endl; }
    #define COSOCOR(cor)
#endif

#define PRINTV(v)                   PRINT(#v << " = " << v)
#define PRINTVEC(vm, v)             PRINT(vm << " = " << v.x << "," << v.y << "," << v.z)
#define PRINTVEC2(v)                PRINT(#v << " = " << v.x << "," << v.y)
#define PRINTV2(v)                  PRINTVEC2(v)
#define PRINTVEC3(v)                PRINT(#v << " = " << (v).x << "," << (v).y << "," << (v).z)
#define PRINTV3(v)                  PRINTVEC3(v)
#define PRINTRECT(v)                PRINT(#v << " = " << v.x << "," << v.y << "," << v.w << "," << v.h)
#define PRINTVEC4(v)                PRINT(#v << " = " << v.x << "," << v.y << "," << v.z << "," << v.w)
#define PRINTQ(v)                   PRINT(#v << " = " << v.w << "," << v.x << "," << v.y << "," << v.z)
#define PRINTTRI(p1,p2,p3)          PRINT("tri(" << p1.x << "," << p1.y << "," << p1.z << "," << p2.x << "," << p2.y << "," << p2.z << "," << p3.x << "," << p3.y << "," << p3.z << ")")
#define PRINTCVP(p)                 PRINT("cvp(" << 50 * p.x + 500 << "," << 50 * p.y + 500 << "," << 50 * p.z + 500 << ")")
#define WARNING(msg)                { COSOCOR(FOREGROUND_RED); std::cerr << "!WARNING: " << msg << std::endl; COSOCOR(FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE); }
#define CAUTION(msg)                { COSOCOR(FOREGROUND_RED | FOREGROUND_GREEN); PRINT("CAUTION: " << msg); COSOCOR(FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE); }
#define NOTICE(msg)                 { COSOCOR(FOREGROUND_GREEN | FOREGROUND_BLUE); PRINT("NOTICE: " << msg); COSOCOR(FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE); }

#define ASSERT(x)                   { if(!(x)) { std::stringstream ss; ss << "ASSERT FAILED! " << __FILE__ << "(" << __LINE__ << "): " << #x; MSGBOX(ss.str()); throw std::runtime_error(ss.str()); } }
#define NASSERT(x)                  ASSERT(!(x))
#define ASSERT_RET(x)               { if(!(x)) { std::stringstream ss; ss << "ASSERT FAILED! " << __FILE__ << "(" << __LINE__ << "): " << #x; PRINT(ss.str()); return 0; } }
#define ERRORMSG(msg)               { COSOCOR(FOREGROUND_RED); std::cerr << msg << std::endl; COSOCOR(FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE); MSGBOX(msg) }

#ifdef WINDOWS
    #define MSGBOX(msg)             { std::stringstream __ss; __ss << msg; ::MessageBoxA(0, __ss.str().c_str(), "MSG", 0); }
    #define MYTRACE(msg)            { std::stringstream __ss; __ss << "\n" << msg; std::cout << msg << std::endl; ::OutputDebugStringA(__ss.str().c_str()); }
    #define ASSERT0(x)              { if(!(x)) { std::stringstream ss; ss << "ASSERT FAILED! " << __FILE__ << "(" << __LINE__ << ")"; ::MessageBoxA(0, ss.str().c_str(), "ASSERT", 0); throw; } }
    #define ASSERT_MSG(x,msg)       { if(!(x)) { ::MessageBoxA(0, msg, "ASSERT", 0); throw std::runtime_error(msg); }}
    #define TICK                    ::GetTickCount()
#else
    #include <chrono>
    #define MSGBOX(msg)             { std::cerr << msg << std::endl; }
    #define MYTRACE(msg)            { std::cerr << msg << std::endl; }
    #define ASSERT0(x)              { if(!(x)) { std::stringstream ss; ss << "ASSERT FAILED! " << __FILE__ << "(" << __LINE__ << ")"; PRINT(ss.str()); } }
    #define ASSERT_MSG(x,msg)       { if(!(x)) throw std::runtime_error(msg); }    
    #define TICK                    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count()
#endif

#define HASHMAP                     std::unordered_map
#define NOW                         TICK

#define COND_THROW(cond, msg)      if(cond){ std::stringstream ss; ss << msg; throw(std::runtime_error(ss.str())); }

#define MAKE_SHAREDT(T)            std::make_shared<T>()
#define MAKE_SHARED(T, X)          std::make_shared<T>(X)

#define INFINITY                    std::numeric_limits<real>::infinity()

#define MAX_FLOAT                   std::numeric_limits<float>::max()

#define MIN_FLOAT                   std::numeric_limits<float>::min()

#define MAX_INT                     2147483647

#define RANDOM_DEV                  std::mt19937  

#define SET_VECEPS(eps)            vec3::sEPSILON = eps

// Common keywords
#define anyptr                      void*
#define crstr                       const std::string&
#define real_type                   float              // Precision type
#define real                        real_type          
#define uint                        unsigned int
#define int64                       long long
#define uchar                       unsigned char
#define crauto                      const auto&
#define rauto                       auto&

#define form                        struct
#define scope_begin(name)          namespace name
#define scope_end                   
#define impfun                      static
#define shrfun                      static inline
#define luaapi                      static int

#define vec3                        vector3
#define rvec                        vector3&
#define crvec                       const vector3&

#define vec2                        vector2
#define rvec2                       vector2&
#define crvec2                      const vector2&
#define ivec2                       point2
#define rivec2                      point2&
#define crivec2                     const point2&
#define ivec3                       point3
#define rivec3                      point3&
#define crivec3                     const point3&

#define pnt                         point2
#define pnt3                        point3

#define vec4                        vector4
#define rvec4                       vector4&
#define crvec4                      const vector4&

#define vecn                        vectorn
#define crvecn                      const vectorn&
#define rvecn                       vectorn&

#define vecN                        vectorN
#define vec10                       vectorN<real, 10>
#define crvecN                      const vectorN&
#define rvecN                       vectorN&

#define dir3                        direction3
#define pos3                        position3

#define quat                        quaternion
#define crquat                      const quaternion&
#define compx                       number_math::complex

#define crd2                        coord2
#define rcd2                        coord2&
#define crcd2                       const coord2&
#define ucd2                        ucoord2
#define csys2                       coord2

#define crd3                        coord3
#define rcd3                        coord3&
#define crcd3                       const coord3&
#define crcrd3                      const coord3&

#define vcd3                        vcoord3
#define rvcd3                       vcoord3&
#define crvcd3                      const vcoord3&

#define ucd3                        ucoord3
#define rucd3                       ucoord3&
#define crucd3                      const ucoord3&

#define csys                        coord3
#define crs3                        coord3

#define AX                          axisX
#define AXZ                         axisXZ
#define AXY                         axisXY

#define mat4                        matrix_lit
#define mat3                        matrix
#define mat2                        mat2x2

#define GVAR                        
#define GFUN                        
#define INVALID                     -1

#define PI                          3.14159265358979323846264338327950288
#define RAD_DEG(rad)                (rad * (PI / 180))

#define delta_d                     0.001f      // Spatial differential precision
#define delta_t                     0.001f      // Time differential precision

#ifndef RGB
#define RGB(r, g, b)               (0xFF000000 | ((unsigned int)((((unsigned char)(r))) | (((unsigned int)(unsigned char)(g)) << 8) | (((unsigned int)(unsigned char)(b)) << 16))))
#endif
inline unsigned int _RGB(unsigned char r, unsigned char g, unsigned char b) { return RGB(r, g, b); }

#define RGBA(r,g,b,a)               ((unsigned int)(((unsigned char)(r)) | (((unsigned int)(unsigned char)(g)) << 8) | (((unsigned int)(unsigned char)(b)) << 16) | (((unsigned int)(unsigned char)(a)) << 24)))

#ifndef GetRValue
#define GetRValue(rgb)             ((unsigned char)(0x000000ff  & (rgb)))
#define GetGValue(rgb)             ((unsigned char)((0x0000ff00 & (rgb))>>8))
#define GetBValue(rgb)             ((unsigned char)((0x00ff0000 & (rgb))>>16))
#endif
#define GetAValue(rgb)             ((unsigned char)((0xff000000 & (rgb))>>24))

#define GETHIGH(v)                 ((int)((v) >> 32))
#define GETLOW(v)                  ((int)(v))
#define HIGHT_LOW(h, l)            (((long long)h << 32) | l)    
#define EPSILON                     1e-3 // Precision

#define ISZERO(a)                  (fabs(a) < 1e-8)

#define MIN(a, b)                  ((a) < (b) ? (a) : (b))
inline real _MIN(real a, real b)   { return ((a) < (b) ? (a) : (b)); }
inline int  _MINi(int a, int b)    { return ((a) < (b) ? (a) : (b)); }

#define MAX(a, b)                  ((a) > (b) ? (a) : (b))
inline real _MAX(real a, real b)   { return ((a) > (b) ? (a) : (b)); }
inline int  _MAXi(int a, int b)    { return ((a) > (b) ? (a) : (b)); }

#define SWAP(A,B)                  {auto t=A; A=B; B=t;}

// Candy
#define DEFAULTVEC3                vec3::ZERO
#define itovec2(ix, iy)            vector2((ix / FIMAGESCALE - 0.5f) * 2.0f, (iy / FIMAGESCALE - 0.5f) * 2.0f)
#define mix                        lerp

#define RNDCOR                     RGB(rrnd(0, 255), rrnd(0, 255), rrnd(0, 255))

#define POLY3                      std::vector<vec3>
#define PLIST                      std::vector<vec3>
#define VECLIST                    std::vector<vec3>
#define VERLIST                    std::vector<vertex>
#define POLY                       std::vector<vertex>
#define EPOLY                      std::vector<vertex>
#define POLYS                      std::vector<POLY>
#define EPOLYS                     std::vector<EPOLY>

#define VECMAP                     std::vector<VECLIST>

#define VECLIST2                  std::vector<vec2>
#define VERLIST2                  std::vector<vertex2>
#define POLY2                     std::vector<vec2>
#define POLY2S                    std::vector<POLY2>
#define EPOLY2                    std::vector<vertex2>
#define EPOLY2S                   std::vector<EPOLY2>

#define IPAIR                     std::pair<int, int>

#define PUSH                      push_back
#define EPUSH                     emplace_back

#define ONE3                      vec3(1,1,1)
#define ZERO3                     vec3::ZERO
#define UNITX                     vec3::UX
#define UNITY                     vec3::UY
#define UNITZ                     vec3::UZ
#define UNIT_X                    vec3::UX
#define UNIT_Y                    vec3::UY
#define UNIT_Z                    vec3::UZ

#define VRED                      vec3::UX
#define VGREEN                    vec3::UY
#define VBLUE                     vec3::UZ
#define VYELLOW                   vec3(1, 1, 0)
#define VGRAY                     vec3(0.5,0.5,0.5)
#define VBLACK                    vec3(0, 0, 0)
#define VWHITE                    vec3(1, 1, 1)

#define RED                       0xFF0000FF
#define DARKRED                   0xFF1818FF
#define ORANGE                    0xFF0080FF
#define GREEN                     0xFF00FF00
#define DARKGREEN                 0xFF008000
#define YELLOW                    0xFF00FFFF
#define GRAYYELLOW                0xFF008080
#define DARKYELLOW                0xFF001818
#define BLUE                      0xFFFF0000
#define GRAYBLUE                  0xFFFF8080
#define LIGHTBLUE                 0xFFFF8080
#define PURPLE                    0xFFFF00FF
#define CYAN                      0xFFFFFF80
#define VIOLET                    0xFF8000FF
#define WHITE                     0xFFFFFFFF
#define BLACK                     0xFF000000
#define GRAY                      0xFF888888
#define DARKGRAY                  0xFF585858

// **********************************************************************
inline real radians(real deg) { return deg / 180.0 * PI; }
inline bool iszero(real x, real tor = 1e-5) { return fabs(x) < tor; }
DEVICE_CALLABLE inline bool is_equal(real a, real b, const real c_EPSILON = EPSILON)
{
    return (::fabs(a - b) <= c_EPSILON);
}
inline real fract(real x) { return x - floor(x); }

// HASH
template <typename T>
inline void hash_combine(std::size_t& seed, const T& v)
{
    seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
inline void hash_combine(std::size_t& seed, const real& v)
{
    seed ^= std::hash<real>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
inline void hash_combine(std::size_t& seed, const char* v)
{
    std::size_t hash = 0;
    while (*v)
    {
        hash = (hash << 5) + hash + *v;
        v++;
    }
    seed ^= hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// ----------------------------------------------------------------------
#include "vector.hpp"
#include "quaternion.hpp"

#include "coord2.hpp"
#include "coord.hpp"
