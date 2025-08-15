/****************************************************
*           Coordinate-System Python API
* 
****************************************************/
namespace python_interface
{
    static PyObject* pModule = nullptr;

    bool _init_python()
    {
        Py_Initialize();
        if (!Py_IsInitialized()) {
            std::cerr << "Python initialization failed!" << std::endl;
            return false;
        }

        // Import necessary modules
        PyRun_SimpleString("import sys");
        PyRun_SimpleString("import numpy as np");

        // Create a module object for the main module
        pModule = PyImport_AddModule("__main__");
        if (!pModule) {
            std::cerr << "Failed to create main module!" << std::endl;
            return false;
        }

        return true;
    }

    void _destroy_python()
    {
        if (Py_IsInitialized()) {
            Py_Finalize();
        }
    }

    std::string dostring(const std::string& pystr)
    {
        if (!pModule)
        {
            if (!_init_python()) {
                return "";
            }
        }

        // Execute the Python string
        PyRun_SimpleString(pystr.c_str());

        // Get the attribute 'cc_out' from the main module
        PyObject* pRet = PyObject_GetAttrString(pModule, "cc_out");
        if (pRet == nullptr) {
            std::cerr << "Attribute 'cc_out' not found!" << std::endl;
            return "";
        }

        // Ensure it's a string
        if (!PyUnicode_Check(pRet)) {
            std::cerr << "'cc_out' is not a string!" << std::endl;
            Py_DECREF(pRet);
            return "";
        }

        // Convert Python string to C++ string
        PyObject* pStr = PyUnicode_AsEncodedString(pRet, "utf-8", "Error");
        if (pStr == nullptr) {
            std::cerr << "Failed to convert Python string to C++ string!" << std::endl;
            Py_DECREF(pRet);
            return "";
        }

        std::string ret = PyBytes_AS_STRING(pStr); // Get the C string
        Py_DECREF(pStr);
        Py_DECREF(pRet);
        return ret;
    }
}

// Coordnate System Export
PYBIND11_MODULE(coordinate_system, m) {
    // Binding for vec3 class
    py::class_<vector3>(m, "vec3")
        .def(py::init<>()) // Default constructor
        .def(py::init<real, real, real>(), "Create a vec3 from x, y, z") // Constructor with x, y, z

        .def_readwrite("x", &vector3::x, "X component") // X component
        .def_readwrite("y", &vector3::y, "Y component") // Y component
        .def_readwrite("z", &vector3::z, "Z component") // Z component

        .def("__repr__", [](const vector3& v) {
                return "<vec3 (" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ")>";
            })

        .def("to_string", [](const vector3& v) {
                return std::to_string(v.x) + "," + std::to_string(v.y) + "," + std::to_string(v.z);
            })
        // Operator Overloads
        .def("__add__", &vector3::operator+, "Addition") // Addition
        .def("__iadd__", (vector3(vector3::*)(const vector3&)) & vector3::operator+=, "In-place addition") // In-place addition
        .def("__sub__", (vector3(vector3::*)(const vector3&) const) & vector3::operator-, "Subtraction") // Subtraction
        .def("__isub__", (vector3(vector3::*)(const vector3&)) & vector3::operator-=, "In-place subtraction") // In-place addition
        .def("__neg__", (vector3(vector3::*)() const) & vector3::operator-, "Negation") // Negation operator
        .def("__mul__", (vector3(vector3::*)(real) const) & vector3::operator*, "Scalar multiplication") // Scalar multiplication
        .def("__mul__", (vector3(vector3::*)(const vector3&) const) & vector3::operator*, "Component-wise multiplication") // Component-wise multiplication
        .def("__mul__", [](const vec3& v, const coord3& coord) {return (v * coord); }, "Multiply vec3 with coord3") // Multiplication with coord3
		.def("__mul__", [](const vec3& v, const quat& q) {return (v * q); }, "Multiply vec3 with quaternion") // Multiplication with quaternion
        .def("__rmul__", (vector3(*)(real, const vector3&)) & operator*, "Right scalar multiplication") // Right scalar multiplication
        .def("__mul__", [](const vec3& v, real s) {return (v * s); }, "Scalar multiplication") // Scalar multiplication
        .def("__truediv__", (vector3(vector3::*)(const vector3&) const) & vector3::operator/, "Component-wise division") // Component-wise division
        .def("__truediv__", [](const vec3& v, const coord3& coord) {return (v / coord); }, "Division vec3 with crd3")
        .def("__truediv__", [](const vector3& v, real s) {return v / s; }, "Scalar division") // Scalar division
        .def("__itruediv__", (vector3(vector3::*)(real)) & vector3::operator/=, "In-place scalar division") // In-place scalar division
        .def("__itruediv__", (vector3(vector3::*)(const vector3&)) & vector3::operator/=, "In-place component-wise division") // In-place component-wise division
        .def("__eq__", (bool (vector3::*)(const vector3&) const) & vector3::operator==, "Equality") // Equality
        .def("__ne__", (bool (vector3::*)(const vector3&) const) & vector3::operator!=, "Inequality") // Inequality

        .def("dot", &vector3::dot, "Dot product with another vec3") // Dot product
        .def("cross", &vector3::cross, "Cross product with another vec3") // Cross product
        .def("cross_right", &vector3::cross_right, "Right Cross product with another vec3") // Cross product
        .def("length", &vector3::length, "Length of the vector") // Length of the vector
        .def("normalize", &vector3::normalize, "Normalize the vector") // Normalize
        .def("isINF", &vector3::isINF, "Check if the vector has any infinite component") // Check for infinite components
        .def("flipX", &vector3::flipX, "Flip the X component") // Flip X component
        .def("flipY", &vector3::flipY, "Flip the Y component") // Flip Y component
        .def("flipZ", &vector3::flipZ, "Flip the Z component") // Flip Z component
        .def("serialise", &vector3::serialise, "Serialize the vector to a string") // Serialize
        .def("len", &vector3::len, "Length of the vector") // Length of the vector
        .def("lenxy", &vector3::lenxy, "Length in the XY plane") // Length in XY plane
        .def("sqrlen", &vector3::sqrlen, "Squared length of the vector") // Squared length
        .def("abslen", &vector3::abslen, "Sum of absolute lengths") // Absolute length
        .def("normalize", &vector3::normalize, "Normalize the vector in place") // Normalize in place
        .def("normalized", &vector3::normalized, "Return a normalized version of the vector") // Normalized version

        // Static methods
        .def_static("min3", &vector3::min3, "Return the minimum of two vec3") // Return minimum of two vec3
        .def_static("max3", &vector3::max3, "Return the maximum of two vec3") // Return maximum of two vec3
        .def_static("rnd", &vector3::rnd, "Generate a random vec3") // Generate random vec3
        .def_static("lerp", (vector3(*)(const vector3&, const vector3&, real)) & vector3::lerp, "Linear interpolation between two vec3") // Linear interpolation
        .def_static("angle", (real(*)(const vector3&, const vector3&)) & vector3::angle, "Angle between two vec3") // Angle between two vec3
        .def_static("angle", (real(*)(const vector3&, const vector3&, const vector3&)) & vector3::angle, "Angle between two vec3 considering an axis") // Angle with axis

        // Accessor methods for vector components
        .def("xxx", &vector3::xxx, "Return vec3 with x, x, x") // xxx
        .def("xxy", &vector3::xxy, "Return vec3 with x, x, y") // xxy
        .def("xxz", &vector3::xxz, "Return vec3 with x, x, z") // xxz
        .def("xyx", &vector3::xyx, "Return vec3 with x, y, x") // xyx
        .def("xyy", &vector3::xyy, "Return vec3 with x, y, y") // xyy
        .def("xyz", &vector3::xyz, "Return vec3 with x, y, z") // xyz
        .def("xzx", &vector3::xzx, "Return vec3 with x, z, x") // xzx
        .def("xzy", &vector3::xzy, "Return vec3 with x, z, y") // xzy
        .def("xzz", &vector3::xzz, "Return vec3 with x, z, z") // xzz
        .def("yxx", &vector3::yxx, "Return vec3 with y, x, x") // yxx
        .def("yxy", &vector3::yxy, "Return vec3 with y, x, y") // yxy
        .def("yxz", &vector3::yxz, "Return vec3 with y, x, z") // yxz
        .def("yyx", &vector3::yyx, "Return vec3 with y, y, x") // yyx
        .def("yyy", &vector3::yyy, "Return vec3 with y, y, y") // yyy
        .def("yyz", &vector3::yyz, "Return vec3 with y, y, z") // yyz
        .def("yzx", &vector3::yzx, "Return vec3 with y, z, x") // yzx
        .def("yzy", &vector3::yzy, "Return vec3 with y, z, y") // yzy
        .def("yzz", &vector3::yzz, "Return vec3 with y, z, z") // yzz
        .def("zxx", &vector3::zxx, "Return vec3 with z, x, x") // zxx
        .def("zxy", &vector3::zxy, "Return vec3 with z, x, y") // zxy
        .def("zxz", &vector3::zxz, "Return vec3 with z, x, z") // zxz
        .def("zyx", &vector3::zyx, "Return vec3 with z, y, x") // zyx
        .def("zyy", &vector3::zyy, "Return vec3 with z, y, y") // zyy
        .def("zyz", &vector3::zyz, "Return vec3 with z, y, z") // zyz
        .def("zzx", &vector3::zzx, "Return vec3 with z, z, x") // zzx
        .def("zzy", &vector3::zzy, "Return vec3 with z, z, y") // zzy
        .def("zzz", &vector3::zzz, "Return vec3 with z, z, z") // zzz
        .def("xyo", &vector3::xyo, "Return vec3 with x, y, 0") // xyo
        .def("xoz", &vector3::xoz, "Return vec3 with x, 0, z") // xoz
        .def("oyz", &vector3::oyz, "Return vec3 with 0, y, z"); // oyz

    // Binding for quaternion class
    py::class_<quaternion>(m, "quat")
        .def(py::init<>()) // Default constructor
        .def(py::init<float, float, float, float>()) // Constructor with w, x, y, z
		.def(py::init<float, float, float>()) // Constructor with pitch,yaw,roll
        .def(py::init<float, const vec3&>()) // Constructor with angle and axis
        .def(py::init<const vec3&, const vec3&>()) // Constructor from two vectors

        .def_readwrite("w", &quaternion::w, "W component") // W component
        .def_readwrite("x", &quaternion::x, "X component") // X component
        .def_readwrite("y", &quaternion::y, "Y component") // Y component
        .def_readwrite("z", &quaternion::z, "Z component") // Z component

        .def("__repr__", [](const quaternion& q) {
            return "<quat (" + std::to_string(q.w) + ", " + std::to_string(q.x) + ", " + std::to_string(q.y) + ", " + std::to_string(q.z) + ")>";
            })
        .def("to_string", [](const quaternion& q) {
                return std::to_string(q.w) + ", " + std::to_string(q.x) + ", " + std::to_string(q.y) + ", " + std::to_string(q.z);
            })

        // Operator Overloads
        .def("__add__", &quaternion::operator+) // Addition
        .def("__sub__", (quaternion(quaternion::*)(const quaternion&) const) & quaternion::operator-, "Subtraction") // Subtraction
        .def("__mul__", (quaternion(quaternion::*)(const quaternion&) const) & quaternion::operator*, "Quaternion multiplication") // Quaternion multiplication
        .def("__mul__", (vector3(quaternion::*)(const vector3&) const) & quaternion::operator*, "Quaternion and vector multiplication") // Quaternion and vector multiplication
        .def("__truediv__", (quaternion(quaternion::*)(const quaternion&) const) & quaternion::operator/, "Division") // Division
		.def("__eq__", (bool (quaternion::*)(const quaternion&) const) & quaternion::operator==, "Equality") // Equality
        .def("__ne__", (bool (quaternion::*)(const quaternion&) const) & quaternion::operator!=, "Inequality") // Inequality
		
        // Member Functions
        .def("normalize", &quaternion::normalize, "Normalize the quaternion") // Normalize
        .def("normalized", &quaternion::normalized, "Return normalized quaternion") // Return normalized quaternion
        .def("angle", &quaternion::angle, "Angle of the quaternion") // Angle
        .def("axis", &quaternion::axis, "Axis of the quaternion") // Axis
        .def("conj", &quaternion::conj, "Conjugate the quaternion") // Conjugate
        .def("conjcopy", &quaternion::conjcopy, "Return conjugate quaternion") // Return conjugate quaternion
        .def("length", &quaternion::length, "Length of quaternion") // Length of quaternion
        .def("dot", &quaternion::dot, "Dot product") // Dot product
        .def("angle_to", &quaternion::angle_to, "Angle to another quaternion") // Angle to another quaternion
        .def("to_eulers", [](const quaternion& q) {
                vec3 eulers = q.toeulers();
                return std::make_tuple(eulers.x, eulers.y, eulers.z);
            }, "Convert to Euler angles") // Convert to Euler angles
        .def_static("slerp", (quaternion(*)(const quaternion&, const quaternion&, float)) & quaternion::slerp, "Spherical linear interpolation") // Spherical linear interpolation
        .def_static("nlerp", (quaternion(*)(real, const quaternion&, const quaternion&, bool)) & quaternion::nlerp, "Normalized linear interpolation") // Normalized linear interpolation
        .def("spherical_cubic_interpolate", &quaternion::spherical_cubic_interpolate, "Spherical cubic interpolation") // Spherical cubic interpolation
        
        .def("is_finite", &quaternion::is_finite, "Check if quaternion components are finite") // Check if quaternion components are finite
        .def("from_vectors", &quaternion::from_vectors, "Create quaternion from two vectors") // Create quaternion from two vectors
        .def("from_eulers", &quaternion::from_eulers, "Create quaternion from Euler angles") // Create quaternion from Euler angles
        .def("from_ang_axis", &quaternion::ang_axis, "Create quaternion from angle and axis") // Create quaternion from angle and axis
        .def("exp", &quaternion::exp, "Exponential of quaternion") // Exponential of quaternion
        .def("log", &quaternion::log, "Logarithm of quaternion"); // Logarithm of quaternion

    // coord3 class binding
    py::class_<coord3>(m, "coord3")
        .def(py::init<>()) // Default constructor
        .def(py::init<real, real, real>()) // Constructor with x, y, z
        .def(py::init<real, real, real, real, real, real>()) // Constructor with six real parameters
        .def(py::init<real, real, real, real, real, real, real>()) // Constructor with six real parameters
        .def(py::init<const vec3&>()) // Constructor from vec3
        .def(py::init<const vec3&, const vec3&, const vec3&, const vec3&>()) // Constructor with multiple vec3 parameters
        .def(py::init<const vec3&, const vec3&, const vec3&, const vec3&, const vec3&>())
        .def(py::init<const vec3&, const vec3&, const vec3&>()) // Constructor with three axis vec3
        .def(py::init<real, const vec3&>()) // Constructor with angle and axis
        .def(py::init<const quaternion&>()) // Constructor from quaternion
        .def(py::init<const vec3&, const quaternion&>()) // Constructor with position, quaternion
        .def(py::init<const vec3&, const quaternion&, const vec3&>()) // Constructor with position, quaternion, and scale

        // Static methods
        .def_static("from_axes",  &coord3::from_axes, "Create a coordinate system from three axes.")
        .def_static("from_angle", &coord3::from_angle, "Create a coordinate system from an angle and axis.")

        // Properties
        .def_readwrite("o",  &coord3::o, "Origin") // Origin
        .def_readwrite("p",  &coord3::o, "Origin") // Origin
        .def_readwrite("x",  &coord3::x, "OriginX") // OriginX
        .def_readwrite("y",  &coord3::y, "OriginY") // OriginY
        .def_readwrite("z",  &coord3::z, "OriginZ") // OriginZ
        .def_readwrite("ux", &coord3::ux, "X-axis") // X-axis
        .def_readwrite("uy", &coord3::uy, "Y-axis") // Y-axis
        .def_readwrite("uz", &coord3::uz, "Z-axis") // Z-axis
        .def_readwrite("s",  &coord3::s, "Scale") // Scale

        .def("__repr__", [](const coord3& c) {
            return "<coord3 (" + std::to_string(c.o.x) + ", " + std::to_string(c.o.y) + ", " + std::to_string(c.o.z) +
                ", " + std::to_string(c.ux.x) + ", " + std::to_string(c.ux.y) + ", " + std::to_string(c.ux.z) +
                ", " + std::to_string(c.uy.x) + ", " + std::to_string(c.uy.y) + ", " + std::to_string(c.uy.z) +
                ", " + std::to_string(c.uz.x) + ", " + std::to_string(c.uz.y) + ", " + std::to_string(c.uz.z) +
                ", " + std::to_string(c.s.x) + ", " + std::to_string(c.s.y) + ", " + std::to_string(c.s.z) + ")>";
                })

        .def("to_string", [](const coord3& c) {
                    quat q = c.Q();
                    return  std::to_string(c.o.x) + ", " + std::to_string(c.o.y) + ", " + std::to_string(c.o.z) + ", " +
                        std::to_string(q.w) + ", " + std::to_string(q.x) + ", " + std::to_string(q.y) + ", " + std::to_string(q.z) + ", " +
                        std::to_string(c.s.x) + ", " + std::to_string(c.s.y) + ", " + std::to_string(c.s.z);
                })

        // Operator Overloads
        .def("__add__", (coord3(coord3::*)(const coord3&) const)& coord3::operator+, "Addition") // Addition
        .def("__add__", (coord3(coord3::*)(const vec3&) const)& coord3::operator+, "Addition with vec3") // Multiplication with vec3
        .def("__iadd__", (coord3(coord3::*)(const coord3&))& coord3::operator+=, "In-place addition") // In-place addition
        .def("__iadd__", (coord3(coord3::*)(const vec3&))& coord3::operator+=, "In-place addition") // In-place addition
        .def("__sub__", (coord3(coord3::*)(const coord3&) const)& coord3::operator-, "Subtraction") // Subtraction
        .def("__sub__", (coord3(coord3::*)(const vec3&) const)& coord3::operator-, "Subtraction with vec3") // Multiplication with vec3
        .def("__isub__", (coord3(coord3::*)(const coord3&) const)& coord3::operator-=, "In-place subtraction") // In-place subtraction
        .def("__isub__", (coord3(coord3::*)(const vec3&))& coord3::operator-=, "In-place subtraction") // In-place addition
        .def("__neg__", (coord3(coord3::*)() const)& coord3::operator-, "Negation") // Negation operator
        .def("__mul__", (coord3(coord3::*)(const coord3&) const)& coord3::operator*, "Multiplication with another coord3") // Multiplication with another coord3
        .def("__mul__", (coord3(coord3::*)(const vec3&) const)& coord3::operator*, "Multiplication with vec3") // Multiplication with vec3
        .def("__mul__", [](const vec3& vec, const coord3& coord) { return operator*(vec, coord); }, "Left multiplication with vec3") // Multiplication with vec3
        .def("__mul__", (coord3(coord3::*)(const quaternion&) const)& coord3::operator*, "Multiply coord3 by quaternion") // Multiplication with quaternion
        .def("__rmul__", [](const vec3& v, const coord3& c) { return c * v; }, "Right multiplication with vec3") // Right multiplication with vec3
        .def("__imul__", (coord3(coord3::*)(const coord3&))& coord3::operator*=, "In-place multiplication") // In-place multiplication
        .def("__imul__", (coord3(coord3::*)(const quaternion&))& coord3::operator*=, "In-place multiplication of coord3 by quaternion") // In-place multiplication with quaternion
        .def("__truediv__", (coord3(coord3::*)(const coord3&) const)& coord3::operator/, "Division by another coord3") // Division by another coord3
        .def("__truediv__", (coord3(coord3::*)(const vec3&) const)& coord3::operator/, "Division by vec3") // Division by vec3
        .def("__itruediv__", (coord3(coord3::*)(const coord3&))& coord3::operator/=, "In-place division") // In-place division
        .def("__eq__", (bool (coord3::*)(const coord3&) const)& coord3::operator==, "Equality") // Equality
        .def("__ne__", (bool (coord3::*)(const coord3&) const)& coord3::operator!=, "Inequality") // Inequality

        // Additional member functions
        .def("Q", [](const coord3& self) { return self.Q(); }, "Convert to quaternion") // Convert to quaternion
        .def("P", &coord3::pos, "Get position") // Get position
        .def("V", &coord3::tovec, "Convert to vec3") // Convert to vec3
        .def("normalize", &coord3::normalize, "Normalize the coordinate system") // Normalize
        .def("to_eulers", &coord3::coord2eulers, "Convert coordinate system to Euler angles") // Convert to Euler angles
        .def("rot", (void(coord3::*)(real, const vec3&))& coord3::rot, "Rotate by angle and axis") // Rotate by angle and axis
        .def("rot", (void(coord3::*)(const quaternion&))& coord3::rot, "Rotate by quaternion") // Rotate by quaternion
        .def("equal_dirs", &coord3::equal_dirs, "Check if directions are equal") // Check if directions are equal
        .def("hash", &coord3::hash, "Get hash") // Get hash
        .def("serialise", &coord3::serialise, "Serialize") // Serialize
        .def("dump", &coord3::dump, "Dump information") // Dump information
        .def("lie_cross", &coord3::lie_cross, "Lie cross product") // Lie cross product
        .def("grad", &coord3::grad, "Gradient") // Gradient
        .def("reverse", &coord3::reverse, "Reverse") // Reverse
        .def("reversed", &coord3::reversed, "Get reversed coord3"); // Get reversed coord3

    // Add bindings for blend function
    m.def("blend", [](const vector3& v1, const vector3& v2, real alpha, real power = 1.0) {
        return blend(v1, v2, alpha, power);
        }, "Blend two vec3 with given alpha and power");

    // Add bindings for slerp functions for vec3
    m.def("slerp", [](const vector3& v1, const vector3& v2, real t, real power = 1) {
        return blender::slerp(v1, v2, t, power);
        }, "Spherical linear interpolation between two vec3");

    // Add bindings for slerp functions for quaternion
    m.def("slerp", [](const quaternion& q1, const quaternion& q2, real t) {
        return blender::slerp(q1, q2, t);
        }, "Spherical linear interpolation between two quaternions");

    m.def("slerp", [](const quaternion& q1, const quaternion& q2, const vec3& ax, real t) {
        return blender::slerp(q1, q2, ax, t);
        }, "Spherical linear interpolation between two quaternions considering an axis");

    // Add bindings for lerp function
    m.def("lerp", [](const coord3& a, const coord3& b, real t) {
        return blender::lerp(a, b, t);
        }, "Linear interpolation between two coord3");

    // Add bindings for lerpPQ function
    m.def("lerpPQ", [](const coord3& a, const coord3& b, real t) {
        return blender::lerpPQ(a, b, t);
        }, "Lerp with position and quaternion interpolation between two coord3");

    // Add bindings for slerp function
    m.def("slerp", [](const coord3& a, const coord3& b, real t) {
        return blender::slerp(a, b, t);
        }, "Spherical linear interpolation between two coord3");

}


