/**
 * coordinate_system_binding.cpp
 *
 * Python bindings for PMSYS Coordinate System Library
 * Global version - English interface
 *
 * **Authors:** Pan Guojun
 * **DOI:** https://doi.org/10.5281/zenodo.14435613
 * Version: 3.0.0 (Right-Handed Coordinate System)
 *
 * Build:
 *   python setup.py build
 *   python setup.py install
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

// Include minimal PMSYS header (vec3, vec2, quat, coord3 only)
// Note: real is defined as float in pmsys_minimal.hpp
#include "pmsys_minimal.hpp"

namespace py = pybind11;

PYBIND11_MODULE(coordinate_system, m) {
    m.doc() = "PMSYS Coordinate System Library - High-performance 3D math for Python";

    // ===========================================================================
    // vec3 class - 3D Vector
    // ===========================================================================
    py::class_<vec3>(m, "vec3", "3D vector class with comprehensive operations")
        .def(py::init<>(), "Default constructor (0, 0, 0)")
        .def(py::init<real, real, real>(), "Construct from x, y, z components",
            py::arg("x"), py::arg("y"), py::arg("z"))

        // Properties
        .def_readwrite("x", &vec3::x, "X component")
        .def_readwrite("y", &vec3::y, "Y component")
        .def_readwrite("z", &vec3::z, "Z component")

        // String representation
        .def("__repr__", [](const vec3& v) {
            return "<vec3(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ")>";
        })
        .def("to_string", [](const vec3& v) {
            return std::to_string(v.x) + "," + std::to_string(v.y) + "," + std::to_string(v.z);
        }, "Convert to comma-separated string")

        // Arithmetic operators
        .def("__add__", &vec3::operator+, "Vector addition")
        .def("__iadd__", [](vec3& self, const vec3& other) { return self += other; }, "In-place addition")
        .def("__sub__", [](const vec3& self, const vec3& other) { return self - other; }, "Vector subtraction")
        .def("__isub__", [](vec3& self, const vec3& other) { return self -= other; }, "In-place subtraction")
        .def("__neg__", [](const vec3& self) { return -self; }, "Negation")
        .def("__mul__", [](const vec3& self, real scalar) { return self * scalar; }, "Scalar multiplication")
        .def("__mul__", [](const vec3& self, const vec3& other) { return self * other; }, "Component-wise multiplication")
        .def("__mul__", [](const vec3& v, const coord3& c) { return v * c; }, "Transform by coordinate system")
        .def("__mul__", [](const vec3& v, const quat& q) { return v * q; }, "Rotate by quaternion")
        .def("__rmul__", [](const vec3& self, real scalar) { return self * scalar; }, "Right scalar multiplication")
        .def("__truediv__", [](const vec3& self, const vec3& other) { return self / other; }, "Component-wise division")
        .def("__truediv__", [](const vec3& self, real scalar) { return self / scalar; }, "Scalar division")
        .def("__truediv__", [](const vec3& v, const coord3& c) { return v / c; }, "Inverse transform by coordinate system")
        .def("__itruediv__", [](vec3& self, real scalar) { return self /= scalar; }, "In-place scalar division")
        .def("__eq__", &vec3::operator==, "Equality comparison")
        .def("__ne__", &vec3::operator!=, "Inequality comparison")

        // Vector operations
        .def("dot", &vec3::dot, "Dot product", py::arg("other"))
        .def("cross", &vec3::cross, "Cross product (right-handed)", py::arg("other"))
        .def("cross_right", &vec3::cross_right, "Cross product (right-handed, alias)", py::arg("other"))
        .def("len", &vec3::len, "Length of vector")
        .def("length", &vec3::length, "Length of vector (alias)")
        .def("lenxy", &vec3::lenxy, "Length in XY plane")
        .def("lenxz", &vec3::lenxz, "Length in XZ plane")
        .def("sqrlen", &vec3::sqrlen, "Squared length")
        .def("sqrlenxy", &vec3::sqrlenxy, "Squared length in XY plane")
        .def("sqrlenxz", &vec3::sqrlenxz, "Squared length in XZ plane")
        .def("len_squared", &vec3::len_squared, "Squared length")
        .def("abslen", &vec3::abslen, "Sum of absolute component values")
        .def("normalize", [](vec3& v) -> vec3& { return vec3_norm(v); }, "Normalize in-place (returns self)")
        .def("normalized", [](const vec3& v) { vec3 copy = v; return vec3_norm(copy); }, "Return normalized copy")
        .def("normcopy", [](const vec3& v) { vec3 copy = v; return vec3_norm(copy); }, "Return normalized copy")
        .def("project", &vec3::project, "Project onto another vector", py::arg("onto"))
        .def("reflect", [](const vec3& v, const vec3& normal) { return vec3_reflect(v, normal); }, "Reflect across normal", py::arg("normal"))
        .def("distance", [](const vec3& v, const vec3& other) { return vec3_distance(v, other); }, "Distance to another point", py::arg("other"))
        // lerp is defined as static method below
        .def("scale", &vec3::scale, "Scale along direction", py::arg("v"), py::arg("scale"))
        .def("mean", &vec3::mean, "Mean of components")
        .def("volum", &vec3::volum, "Volume (product of absolute values)")
        .def("hash", &vec3::hash, "Hash value", py::arg("precision_level") = 0)

        // Boolean checks
        .def("isINF", &vec3::isINF, "Check if any component is infinite")

        // Axis flipping
        .def("flipX", &vec3::flipX, "Flip X component")
        .def("flipY", &vec3::flipY, "Flip Y component")
        .def("flipZ", &vec3::flipZ, "Flip Z component")

        // Swizzling operations
        .def("xxx", &vec3::xxx) .def("xxy", &vec3::xxy) .def("xxz", &vec3::xxz)
        .def("xyx", &vec3::xyx) .def("xyy", &vec3::xyy) .def("xyz", &vec3::xyz)
        .def("xzx", &vec3::xzx) .def("xzy", &vec3::xzy) .def("xzz", &vec3::xzz)
        .def("yxx", &vec3::yxx) .def("yxy", &vec3::yxy) .def("yxz", &vec3::yxz)
        .def("yyx", &vec3::yyx) .def("yyy", &vec3::yyy) .def("yyz", &vec3::yyz)
        .def("yzx", &vec3::yzx) .def("yzy", &vec3::yzy) .def("yzz", &vec3::yzz)
        .def("zxx", &vec3::zxx) .def("zxy", &vec3::zxy) .def("zxz", &vec3::zxz)
        .def("zyx", &vec3::zyx) .def("zyy", &vec3::zyy) .def("zyz", &vec3::zyz)
        .def("zzx", &vec3::zzx) .def("zzy", &vec3::zzy) .def("zzz", &vec3::zzz)
        .def("xyo", &vec3::xyo) .def("xoz", &vec3::xoz) .def("oyz", &vec3::oyz)

        // Static methods
        .def_static("min3", &vec3::min3, "Component-wise minimum", py::arg("a"), py::arg("b"))
        .def_static("max3", &vec3::max3, "Component-wise maximum", py::arg("a"), py::arg("b"))
        .def_static("rnd", &vec3::rnd, "Random vector")
        .def_static("lerp", [](const vec3& a, const vec3& b, real t) { return vec3::lerp(a, b, t); },
            "Linear interpolation", py::arg("a"), py::arg("b"), py::arg("t"))
        .def_static("angle", [](const vec3& a, const vec3& b) { return vec3::angle(a, b); },
            "Angle between two vectors", py::arg("a"), py::arg("b"))
        .def_static("angle", [](const vec3& a, const vec3& b, const vec3& axis) { return vec3::angle(a, b, axis); },
            "Angle between two vectors with axis", py::arg("a"), py::arg("b"), py::arg("axis"))

        // Static constants
        .def_readonly_static("ZERO", &vec3::ZERO, "Zero vector (0, 0, 0)")
        .def_readonly_static("UX", &vec3::UX, "Unit X axis (1, 0, 0)")
        .def_readonly_static("UY", &vec3::UY, "Unit Y axis (0, 1, 0)")
        .def_readonly_static("UZ", &vec3::UZ, "Unit Z axis (0, 0, 1)");

    // ===========================================================================
    // vec2 class - 2D Vector
    // ===========================================================================
    py::class_<vec2>(m, "vec2", "2D vector class")
        .def(py::init<>())
        .def(py::init<real, real>(), py::arg("x"), py::arg("y"))
        .def(py::init<real>(), py::arg("value"))

        .def_readwrite("x", &vec2::x)
        .def_readwrite("y", &vec2::y)

        .def("__repr__", [](const vec2& v) {
            return "<vec2(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ")>";
        })
        .def("to_string", [](const vec2& v) {
            return std::to_string(v.x) + "," + std::to_string(v.y);
        })

        .def("__add__", &vec2::operator+)
        .def("__iadd__", [](vec2& self, const vec2& other) { return self += other; })
        .def("__sub__", [](const vec2& self, const vec2& other) { return self - other; })
        .def("__neg__", [](const vec2& self) { return -self; })
        .def("__mul__", [](const vec2& self, real s) { return self * s; })
        .def("__mul__", [](const vec2& self, const vec2& other) { return self * other; })
        .def("__rmul__", [](const vec2& self, real s) { return self * s; })
        .def("__truediv__", [](const vec2& self, real s) { return self / s; })
        .def("__truediv__", [](const vec2& self, const vec2& other) { return self / other; })
        .def("__eq__", &vec2::operator==)
        .def("__ne__", &vec2::operator!=)
        .def("__lt__", &vec2::operator<)
        .def("__le__", &vec2::operator<=)
        .def("__gt__", &vec2::operator>)
        .def("__ge__", &vec2::operator>=)

        .def("dot", &vec2::dot)
        .def("cross", &vec2::cross, "2D cross product (returns scalar)")
        .def("len", &vec2::len)
        .def("length", &vec2::length)
        .def("sqrlen", &vec2::sqrlen)
        .def("normalize", [](vec2& v) -> vec2& { return vec2_norm(v); })
        .def("normalized", [](const vec2& v) { vec2 copy = v; return vec2_norm(copy); })
        .def("normcopy", [](const vec2& v) { vec2 copy = v; return vec2_norm(copy); })
        .def("angle", [](const vec2& self) { return self.angle(); }, "Get angle")
        .def("rot", [](vec2& self, real angle) { self.rot(angle); }, "Rotate by angle")
        .def("rotcopy", [](const vec2& v, real angle) { return v.rotcopy(angle); }, "Return rotated copy", py::arg("angle"))
        .def("rotcopy", [](const vec2& v, real angle, const vec2& o) { return v.rotcopy(angle, o); }, "Return rotated copy around origin", py::arg("angle"), py::arg("origin"))
        .def("roted", [](const vec2& v, real angle) { return v.roted(angle); }, "Return rotated copy", py::arg("angle"))
        .def("roted", [](const vec2& v, real angle, const vec2& o) { return v.roted(angle, o); }, "Return rotated copy around origin", py::arg("angle"), py::arg("origin"))
        .def("isINF", &vec2::isINF)
        .def("distance", [](const vec2& v, const vec2& other) { return vec2_distance(v, other); }, py::arg("other"))
        .def("lerp", [](const vec2& v, const vec2& other, real t) { return vec2_lerp(v, other, t); }, py::arg("other"), py::arg("t"))

        .def("xx", &vec2::xx) .def("xy", &vec2::xy)
        .def("yx", &vec2::yx) .def("yy", &vec2::yy)

        .def_static("ang_len", &vec2::ang_len, "Create from angle and length", py::arg("angle"), py::arg("length"));

    // ===========================================================================
    // quaternion class
    // ===========================================================================
    py::class_<quat>(m, "quat", "Quaternion for 3D rotations")
        .def(py::init<>(), "Identity quaternion")
        .def(py::init<real, real, real, real>(), "Construct from w, x, y, z",
            py::arg("w"), py::arg("x"), py::arg("y"), py::arg("z"))
        .def(py::init<real, real, real>(), "Construct from Euler angles (pitch, yaw, roll)",
            py::arg("pitch"), py::arg("yaw"), py::arg("roll"))
        .def(py::init<real, const vec3&>(), "Construct from angle and axis",
            py::arg("angle"), py::arg("axis"))
        .def(py::init<const vec3&, const vec3&>(), "Construct from two vectors",
            py::arg("from"), py::arg("to"))

        .def_readwrite("w", &quat::w)
        .def_readwrite("x", &quat::x)
        .def_readwrite("y", &quat::y)
        .def_readwrite("z", &quat::z)

        .def("__repr__", [](const quat& q) {
            return "<quat(w=" + std::to_string(q.w) + ", x=" + std::to_string(q.x) +
                   ", y=" + std::to_string(q.y) + ", z=" + std::to_string(q.z) + ")>";
        })
        .def("to_string", [](const quat& q) {
            return std::to_string(q.w) + ", " + std::to_string(q.x) + ", " +
                   std::to_string(q.y) + ", " + std::to_string(q.z);
        })

        .def("__add__", &quat::operator+)
        .def("__sub__", [](const quat& self, const quat& other) { return self - other; })
        .def("__mul__", [](const quat& self, const quat& other) { return self * other; }, "Quaternion multiplication")
        .def("__mul__", [](const quat& self, const vec3& v) { return self * v; }, "Rotate vector")
        .def("__rmul__", [](const quat& self, const vec3& v) { return v * self; })
        .def("__truediv__", [](const quat& self, const quat& other) { return self / other; })
        .def("__eq__", &quat::operator==)
        .def("__ne__", &quat::operator!=)
        .def("__pow__", [](const quat& self, int n) { return self ^ n; }, "Integer power")
        .def("__pow__", [](quat& self, real n) { return self ^ n; }, "Real power")

        .def("normalize", [](quat& q) -> quat& { return quat_norm(q); }, "Normalize in-place")
        .def("normalized", [](const quat& q) { return quat_normcopy(q); }, "Return normalized copy")
        .def("angle", &quat::angle, "Get rotation angle")
        .def("axis", &quat::axis, "Get rotation axis")
        .def("conj", &quat::conj, "Conjugate in-place")
        .def("conjcopy", &quat::conjcopy, "Return conjugate")
        .def("inverse", &quat::inverse, "Return inverse")
        .def("length", &quat::length)
        .def("dot", &quat::dot, py::arg("other"))
        .def("angle_to", &quat::angle_to, "Angle to another quaternion", py::arg("other"))
        .def("rotate", [](const quat& q, const vec3& v) { return quat_rotate(q, v); }, "Rotate a vector", py::arg("v"))

        .def("to_eulers", [](const quat& self) {
            vec3 eulers = self.toeulers();
            return py::make_tuple(eulers.x, eulers.y, eulers.z);
        }, "Convert to Euler angles (pitch, yaw, roll)")
        .def("to_angle_axis", [](const quat& self) {
            real angle;
            vec3 axis;
            self.to_angle_axis(angle, axis);
            return py::make_tuple(angle, axis);
        }, "Get angle and axis")

        .def("xyz", &quat::xyz, "Get XYZ components as vec3")
        .def("set_angle", &quat::set_angle, "Set rotation angle", py::arg("angle"))
        .def("rotate_x", &quat::rotate_x, "Rotate around X axis", py::arg("angle"))
        .def("rotate_y", &quat::rotate_y, "Rotate around Y axis", py::arg("angle"))
        .def("rotate_z", &quat::rotate_z, "Rotate around Z axis", py::arg("angle"))
        .def("is_finite", &quat::is_finite)
        .def("from_vectors", &quat::from_vectors, "Create from two vectors", py::arg("from"), py::arg("to"))
        .def("from_eulers", &quat::from_eulers, "Set from Euler angles",
            py::arg("pitch"), py::arg("yaw"), py::arg("roll"))
        .def("log", &quat::log, "Logarithm")

        .def_static("slerp", [](const quat& a, const quat& b, real t) { return quat::slerp(a, b, t); },
            "Spherical linear interpolation", py::arg("a"), py::arg("b"), py::arg("t"))
        .def_static("nlerp", [](real t, const quat& a, const quat& b, bool shortest = true) {
            return quat::nlerp(t, a, b, shortest);
        }, "Normalized linear interpolation", py::arg("t"), py::arg("a"), py::arg("b"), py::arg("shortest") = true)
        .def_static("from_euler", [](real pitch, real yaw, real roll) {
            quat q;
            q.from_eulers(pitch, yaw, roll);
            return q;
        }, "Create from Euler angles", py::arg("pitch"), py::arg("yaw"), py::arg("roll"))
        .def_static("from_axis_angle", [](const vec3& axis, real angle) {
            return quat(angle, axis);
        }, "Create from axis and angle", py::arg("axis"), py::arg("angle"));

    // ===========================================================================
    // coord3 class - 3D Coordinate System
    // ===========================================================================
    py::class_<coord3>(m, "coord3", "3D coordinate system (frame)")
        .def(py::init<>(), "Identity coordinate system")
        .def(py::init<real, real, real>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def(py::init<real, real, real, real, real, real>(), "From position and rotation angles")
        .def(py::init<real, real, real, real, real, real, real>())
        .def(py::init<const vec3&>(), "From position", py::arg("position"))
        .def(py::init<const vec3&, const vec3&, const vec3&, const vec3&>(),
            "From origin and axes", py::arg("o"), py::arg("ux"), py::arg("uy"), py::arg("uz"))
        .def(py::init<const vec3&, const vec3&, const vec3&, const vec3&, const vec3&>())
        .def(py::init<const vec3&, const vec3&, const vec3&>(), "From three axes")
        .def(py::init<real, const vec3&>(), "From angle and axis", py::arg("angle"), py::arg("axis"))
        .def(py::init<const quat&>(), "From quaternion", py::arg("q"))
        .def(py::init<const vec3&, const quat&>(), "From position and quaternion",
            py::arg("position"), py::arg("rotation"))
        .def(py::init<const vec3&, const quat&, const vec3&>(), "From position, rotation, and scale")

        // Properties
        .def_readwrite("o", &coord3::o, "Origin/position")
        .def_readwrite("p", &coord3::o, "Origin/position (alias)")
        .def_readwrite("x", &coord3::x, "Origin X component")
        .def_readwrite("y", &coord3::y, "Origin Y component")
        .def_readwrite("z", &coord3::z, "Origin Z component")
        .def_readwrite("ux", &coord3::ux, "X axis direction")
        .def_readwrite("uy", &coord3::uy, "Y axis direction")
        .def_readwrite("uz", &coord3::uz, "Z axis direction")
        .def_readwrite("s", &coord3::s, "Scale")

        .def("__repr__", [](const coord3& c) {
            return "<coord3(origin=" + std::to_string(c.o.x) + "," + std::to_string(c.o.y) + "," + std::to_string(c.o.z) + ")>";
        })
        .def("to_string", [](const coord3& c) {
            quat q = c.Q();
            return std::to_string(c.o.x) + ", " + std::to_string(c.o.y) + ", " + std::to_string(c.o.z) + ", " +
                   std::to_string(q.w) + ", " + std::to_string(q.x) + ", " + std::to_string(q.y) + ", " + std::to_string(q.z) + ", " +
                   std::to_string(c.s.x) + ", " + std::to_string(c.s.y) + ", " + std::to_string(c.s.z);
        })

        // Operators
        .def("__add__", [](const coord3& self, const coord3& other) { return self + other; })
        .def("__add__", [](const coord3& self, const vec3& v) { return self + v; })
        .def("__iadd__", [](coord3& self, const coord3& other) { return self += other; })
        .def("__iadd__", [](coord3& self, const vec3& v) { return self += v; })
        .def("__sub__", [](const coord3& self, const coord3& other)  { return self - other; })
        .def("__sub__", [](const coord3& self, const vec3& v) { return self - v; })
        .def("__isub__", [](coord3& self, const coord3& other) { return self -= other; })
        .def("__isub__", [](coord3& self, const vec3& v) { return self -= v; })
        .def("__neg__", [](const coord3& self) { return -self; })
        .def("__mul__", [](const coord3& self, const coord3& other) { return self * other; })
		.def("__mul__", [](const coord3& self, const real& scaler) { return self * scaler; })
        .def("__mul__", [](const coord3& self, const vec3& v) { return self * v; })
        .def("__mul__", [](const vec3& v, const coord3& c) { return v * c; })
        .def("__mul__", [](const coord3& self, const quat& q) { return self * q; })
        .def("__rmul__", [](const coord3& self, const vec3& v) { return self * v; })
        .def("__imul__", [](coord3& self, const coord3& other) { return self *= other; })
        .def("__imul__", [](coord3& self, const quat& q) { return self *= q; })
        .def("__truediv__", [](const coord3& self, const coord3& other) { return self / other; })
        .def("__truediv__", [](const coord3& self, const vec3& v) { return self / v; })
		.def("__truediv__", [](const coord3& self, const real& scaler) { return self / scaler; })
        .def("__itruediv__", [](coord3& self, const coord3& other) { return self /= other; })
        .def("__eq__", &coord3::operator==)
        .def("__ne__", &coord3::operator!=)

        // Methods
        .def("Q", [](const coord3& self) { return self.Q(); }, "Get rotation as quaternion")
		.def("R", [](const coord3& self) { return coord3(self.R()); }, "Get rotation as coord3")
        .def("P", &coord3::pos, "Get position")
        .def("V", &coord3::tovec, "Convert to vector")
        .def("normalize", &coord3::normalize, "Normalize axes")
        .def("to_eulers", &coord3::coord2eulers, "Convert to Euler angles")
        .def("rot", [](coord3& self, real angle, const vec3& axis) { self.rot(angle, axis); },
            "Rotate by angle and axis", py::arg("angle"), py::arg("axis"))
        .def("rot", [](coord3& self, const quat& q) { self.rot(q); }, "Rotate by quaternion", py::arg("q"))
        .def("equal_dirs", &coord3::equal_dirs, "Check if directions are equal", py::arg("other"))
        // .def("hash", &coord3::hash, "Get hash")  // hash() method not available in coord3
        .def("dump", &coord3::dump, "Dump information")
        .def("lie_cross", &coord3::lie_cross, py::arg("other"))
        .def("grad", &coord3::grad, py::arg("other"))		
		.def("inverse", &coord3::reverse, "Inverse")
		.def("inversed", &coord3::reversed, "Get Inversed copy")
        .def("reverse", &coord3::reverse, "Reverse")
        .def("reversed", &coord3::reversed, "Get reversed copy")
        .def("distance_to", &coord3::distance_to, "Distance to another coordinate system", py::arg("other"))
        .def("rotation_distance_to", &coord3::rotation_distance_to, "Rotation angle difference", py::arg("other"))
        .def("pose", &coord3::pose, "Get pose (rotation + position)")
        .def("pos", &coord3::pos, "Get position")
        .def("to_world", [](const coord3& c, const vec3& local) { return coord3_to_world(c, local); }, "Transform to world coordinates", py::arg("local"))
        .def("to_local", [](const coord3& c, const vec3& world) { return coord3_to_local(c, world); }, "Transform to local coordinates", py::arg("world"))

        // Advanced accessors
        .def("VX", static_cast<vec3(coord3::*)() const>(&coord3::VX), "Get scaled X axis")
		.def("VY", static_cast<vec3(coord3::*)() const>(&coord3::VY), "Get scaled Y axis")
		.def("VZ", static_cast<vec3(coord3::*)() const>(&coord3::VZ), "Get scaled Z axis")
        .def("X", &coord3::X, "Get X axis with position offset")
        .def("Y", &coord3::Y, "Get Y axis with position offset")
        .def("Z", &coord3::Z, "Get Z axis with position offset")
        .def("ucoord", [](const coord3& self) { return self.ucoord(); }, "Get rotation part")
        .def("UC", [](const coord3& self) { return self.UC(); }, "Get rotation part (alias)")
        .def("R", [](const coord3& self) { return self.R(); }, "Get rotation part (alias)")
        .def("VC", &coord3::VC, "Get rotation+scale part")

        // Correct metric determinant calculation for intrinsic gradient operator
        .def("compute_metric_det", [](const coord3& self) {
            // Correct formula: det(G) = det(F^T * Λ^2 * F)
            // where F = [ux, uy, uz] and Λ = diag(s)

            // For 2D surface embedded in 3D, we only use first two components
            real E = self.s.x * self.s.x;  // |∂r/∂u|^2
            real G = self.s.y * self.s.y;  // |∂r/∂v|^2
            real F = self.s.x * self.s.y * self.ux.dot(self.uy);  // ∂r/∂u · ∂r/∂v

            return E * G - F * F;
        }, "Compute metric determinant for surface (2D manifold)")

        // Static methods
        .def_static("from_axes", &coord3::from_axes, "Create from three axes")
        .def_static("from_angle", &coord3::from_angle, "Create from angle and axis")
        .def_static("look_at", &coord3::look_at, "Create look-at coordinate system",
            py::arg("eye"), py::arg("target"), py::arg("up") = vec3::UY)
        .def_static("from_forward", &coord3::from_forward, "Create from position and forward direction",
            py::arg("pos"), py::arg("forward"), py::arg("up") = vec3::UY)
        .def_static("from_eulers", &coord3::from_eulers, "Create from Euler angles")
        .def_static("lerp", &coord3::lerp, "Linear interpolation", py::arg("a"), py::arg("b"), py::arg("t"))
        .def_static("slerp", &coord3::slerp, "Spherical linear interpolation", py::arg("a"), py::arg("b"), py::arg("t"))

        // Common constructors
        .def_static("identity", []() {
            return coord3(vec3::ZERO, vec3::UX, vec3::UY, vec3::UZ);
        }, "Create identity coordinate system")
        .def_static("zero", []() {
            return coord3(vec3::ZERO, vec3::UX, vec3::UY, vec3::UZ);
        }, "Create zero coordinate system (same as identity)")
        .def_static("from_position", [](const vec3& pos) {
            return coord3(pos, vec3::UX, vec3::UY, vec3::UZ);
        }, "Create coordinate system at position with identity rotation", py::arg("position"))
        .def_static("from_rotation", [](const quat& q) {
            return coord3(vec3::ZERO, q);
        }, "Create coordinate system at origin with rotation", py::arg("rotation"))

        // Differential Geometry
		.def("metric", &coord3::metric,
            "Compute metric tensor determinant: det(g) = E*G - F²\n\n"
            "For a parametric surface, the first fundamental form is g = [[E, F], [F, G]].\n"
            "Returns det(g) = E*G - F² which appears in Gaussian curvature formula.")  
			
        .def("metric_det", &coord3::metric_det,
            "Compute metric tensor determinant: det(g) = E*G - F²\n\n"
            "For a parametric surface, the first fundamental form is g = [[E, F], [F, G]].\n"
            "Returns det(g) = E*G - F² which appears in Gaussian curvature formula.")        
        
        .def_static("lie_bracket", &coord3::lie_bracket,
            "Compute Lie bracket [A, B] = A * B - B * A\n\n"
            "The Lie bracket measures the non-commutativity of the connection operators.\n"
            "It appears in the curvature tensor formula as the Lie derivative term.\n\n"
            "Args:\n"
            "    A: connection operator in u-direction\n"
            "    B: connection operator in v-direction\n\n"
            "Returns:\n"
            "    [A, B] (Lie bracket as coord3)",
            py::arg("A"), py::arg("B"));

    // ===========================================================================
    // Utility Functions
    // ===========================================================================

    // Blend functions
    m.def("blend", [](const vec3& v1, const vec3& v2, real alpha, real power) {
        return blend(v1, v2, alpha, power);
    }, "Blend two vectors", py::arg("v1"), py::arg("v2"), py::arg("alpha"), py::arg("power") = 1.0f);

    // Slerp functions
    m.def("slerp", [](const vec3& v1, const vec3& v2, real t) {
        return blender::slerp(v1, v2, t);
    }, "Spherical linear interpolation (vec3)", py::arg("v1"), py::arg("v2"), py::arg("t"));

    m.def("slerp", [](const quat& q1, const quat& q2, real t) {
        return blender::slerp(q1, q2, t);
    }, "Spherical linear interpolation (quat)", py::arg("q1"), py::arg("q2"), py::arg("t"));

    m.def("slerp", [](const coord3& a, const coord3& b, real t) {
        return blender::slerp(a, b, t);
    }, "Spherical linear interpolation (coord3)", py::arg("a"), py::arg("b"), py::arg("t"));

    // Lerp functions
    m.def("lerp", [](const coord3& a, const coord3& b, real t) {
        // Use linear interpolation (position + rotation + scale all linear).
        vec3 pos = a.o * (1.0 - t) + b.o * t;
        quat q = quat::slerp(a.Q(), b.Q(), t);  // Rotation uses slerp.
        vec3 scale = a.s * (1.0 - t) + b.s * t;
        return coord3(pos, q, scale);
    }, "Linear interpolation (coord3)", py::arg("a"), py::arg("b"), py::arg("t"));

    m.def("lerpPQ", [](const coord3& a, const coord3& b, real t) {
        return blender::lerpPQ(a, b, t);
    }, "Lerp with position and quaternion", py::arg("a"), py::arg("b"), py::arg("t"));

    // Vector operations
    m.def("dot", [](const vec3& a, const vec3& b) { return a.dot(b); },
        "Dot product", py::arg("a"), py::arg("b"));
    m.def("cross", [](const vec3& a, const vec3& b) { return a.cross(b); },
        "Cross product (right-handed)", py::arg("a"), py::arg("b"));
    m.def("cross_right", [](const vec3& a, const vec3& b) { return a.cross_right(b); },
        "Cross product (right-handed, alias)", py::arg("a"), py::arg("b"));
    m.def("distance", [](const vec3& a, const vec3& b) { return (b - a).len(); },
        "Distance between two vectors", py::arg("a"), py::arg("b"));
    m.def("angle_between", [](const vec3& a, const vec3& b) { return vec3::angle(a, b); },
        "Angle between two vectors", py::arg("a"), py::arg("b"));
    m.def("angle_between", [](const vec3& a, const vec3& b, const vec3& axis) {
        return vec3::angle(a, b, axis);
    }, "Angle between two vectors with axis", py::arg("a"), py::arg("b"), py::arg("axis"));

    // Spherical coordinates
    m.def("vec3_from_spherical", [](real theta, real phi, real r) {
        return vec3(r * sin(phi) * cos(theta),
                    r * sin(phi) * sin(theta),
                    r * cos(phi));
    }, "Create vec3 from spherical coordinates", py::arg("theta"), py::arg("phi"), py::arg("r") = 1.0f);

    m.def("vec3_to_spherical", [](const vec3& v) {
        real r = v.len();
        real theta = atan2(v.y, v.x);
        real phi = acos(v.z / r);
        return py::make_tuple(theta, phi, r);
    }, "Convert vec3 to spherical coordinates", py::arg("v"));

    // Coordinate system transforms
    m.def("transform_point", [](const vec3& point, const coord3& coord) {
        return point * coord;
    }, "Transform point by coordinate system", py::arg("point"), py::arg("coord"));

    m.def("inverse_transform_point", [](const vec3& point, const coord3& coord) {
        return point / coord;
    }, "Inverse transform point", py::arg("point"), py::arg("coord"));

    m.def("transform_vector", [](const vec3& vec, const coord3& coord) {
        vec3 result = vec;
        result = result.x * coord.ux + result.y * coord.uy + result.z * coord.uz;
        return result;
    }, "Transform vector (no translation)", py::arg("vector"), py::arg("coord"));

    // Quaternion utilities
    m.def("quat_from_two_vectors", [](const vec3& from, const vec3& to) {
        quat q;
        q.from_vectors(from, to);
        return q;
    }, "Create quaternion from two vectors", py::arg("from"), py::arg("to"));

    m.def("quat_look_rotation", [](const vec3& forward, const vec3& up) {
        vec3 zaxis = forward.normcopy();
        vec3 xaxis = up.cross(zaxis).normcopy();
        vec3 yaxis = zaxis.cross(xaxis);
        ucoord3 uc(xaxis, yaxis, zaxis);
        return uc.Q();
    }, "Create quaternion for look rotation", py::arg("forward"), py::arg("up") = vec3::UY);

    // Constants
    m.attr("PI") = PI;
    m.attr("EPSILON") = (real)EPSILON;
    m.attr("VERSION") = GCU_VERSION;

    m.attr("__version__") = "4.0.1";
    m.attr("__author__") = "PanGuoJun";
    m.attr("__coordinate_system__") = "right-handed";

    // ===========================================================================
    // Coordinate System Handedness Control
    // ===========================================================================
    m.def("set_handedness", [](const std::string& handedness) {
        if (handedness == "left" || handedness == "LEFT") {
            set_coordinate_handedness(true);
        } else if (handedness == "right" || handedness == "RIGHT") {
            set_coordinate_handedness(false);
        } else {
            throw std::invalid_argument("handedness must be 'left' or 'right'");
        }
    }, "Set coordinate system handedness ('left' or 'right')", py::arg("handedness"));

    m.def("get_handedness", []() -> std::string {
        return get_coordinate_handedness() ? "left" : "right";
    }, "Get current coordinate system handedness");

    m.def("is_left_handed", []() -> bool {
        return get_coordinate_handedness();
    }, "Check if using left-handed coordinate system");

    m.def("is_right_handed", []() -> bool {
        return !get_coordinate_handedness();
    }, "Check if using right-handed coordinate system");

}
