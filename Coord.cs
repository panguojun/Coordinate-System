/// <summary>
/// ***************************************************************************************
/// Represents a three-dimensional coordinate system in Unity.
/// ***************************************************************************************
/// </summary>
[System.Serializable]
public struct UCoord3
{
    public static readonly UCoord3 zero = new UCoord3(Vector3.zero, Vector3.zero, Vector3.zero);
    public static readonly UCoord3 one = new UCoord3(Vector3.right, Vector3.up, Vector3.forward);

    public Vector3 right;
    public Vector3 up;
    public Vector3 forward;

    /// <summary>
    /// Creates a new UCoord3 instance with the specified vectors as the axes.
    /// </summary>
    /// <param name="right">The right-axis vector.</param>
    /// <param name="up">The up-axis vector.</param>
    /// <param name="forward">The forward-axis vector.</param>
    public UCoord3(Vector3 right, Vector3 up, Vector3 forward)
    {
        this.right = right;
        this.up = up;
        this.forward = forward;
    }

    /// <summary>
    /// Creates a new UCoord3 instance with the specified vectors as the axes and calculates the forward-axis vector.
    /// </summary>
    /// <param name="right">The right-axis vector.</param>
    /// <param name="up">The up-axis vector.</param>
    public UCoord3(Vector3 right, Vector3 up)
    {
        this.right = right;
        this.up = up;
        this.forward = Vector3.Cross(right, up);
    }

    /// <summary>
    /// Creates a new UCoord3 instance with the specified rotation angle around the specified axis.
    /// </summary>
    /// <param name="angle">The rotation angle in degrees.</param>
    /// <param name="axis">The rotation axis vector.</param>
    public UCoord3(float angle, Vector3 axis)
    {
        Quaternion q = Quaternion.AngleAxis(angle, axis);

        this.right = q * Vector3.right;
        this.up = q * Vector3.up;
        this.forward = q * Vector3.forward;
    }
    public UCoord3(Quaternion q)
    {
        this.right = q * Vector3.right;
        this.up = q * Vector3.up;
        this.forward = q * Vector3.forward;
    }
    /// <summary>
    /// Sets the UCoord3 instance based on the rotation between two vectors.
    /// </summary>
    /// <param name="v1">The first vector.</param>
    /// <param name="v2">The second vector.</param>
    public void FromVecsR(Vector3 v1, Vector3 v2)
    {
        Quaternion q = Quaternion.FromToRotation(v1, v2);
        this.right = q * Vector3.right;
        this.up = q * Vector3.up;
        this.forward = q * Vector3.forward;
    }

    /// <summary>
    /// Converts the UCoord3 instance to a Quaternion.
    /// </summary>
    /// <returns>The Quaternion representation of the UCoord3 instance.</returns>
    public Quaternion ToQuaternion()
    {
        Vector3 eulers = ToEulerAngles();
        Quaternion q = Quaternion.Euler(eulers);
        return q;
    }

    /// <summary>
    /// Checks if the current UCoord3 instance has the same direction as the specified UCoord3.
    /// </summary>
    /// <param name="c">The UCoord3 to compare with.</param>
    /// <returns>True if the directions are the same, false otherwise.</returns>
    public bool SameDirections(UCoord3 c)
    {
        return this.right == c.right && this.up == c.up && this.forward == c.forward;
    }

    /// <summary>
    /// Transforms the specified vector using the UCoord3 instance.
    /// </summary>
    /// <param name="p">The vector to transform.</param>
    /// <param name="c">The UCoord3 instance.</param>
    /// <returns>The transformed vector.</returns>
    public static Vector3 operator *(Vector3 p, UCoord3 c)
    {
        return c.right * p.x + c.up * p.y + c.forward * p.z;
    }

    /// <summary>
    /// Multiplies two UCoord3 instances together.
    /// </summary>
    /// <param name="c1">The first UCoord3.</param>
    /// <param name="c2">The second UCoord3.</param>
    /// <returns>The result of the multiplication.</returns>
    public static UCoord3 operator *(UCoord3 c1, UCoord3 c2)
    {
        UCoord3 rc = new UCoord3();
        rc.right = new Vector3(c1.right.x * c2.right.x + c1.right.y * c2.up.x + c1.right.z * c2.forward.x,
                               c1.right.x * c2.right.y + c1.right.y * c2.up.y + c1.right.z * c2.forward.y,
                               c1.right.x * c2.right.z + c1.right.y * c2.up.z + c1.right.z * c2.forward.z);
        rc.up = new Vector3(c1.up.x * c2.right.x + c1.up.y * c2.up.x + c1.up.z * c2.forward.x,
                            c1.up.x * c2.right.y + c1.up.y * c2.up.y + c1.up.z * c2.forward.y,
                            c1.up.x * c2.right.z + c1.up.y * c2.up.z + c1.up.z * c2.forward.z);
        rc.forward = new Vector3(c1.forward.x * c2.right.x + c1.forward.y * c2.up.x + c1.forward.z * c2.forward.x,
                                 c1.forward.x * c2.right.y + c1.forward.y * c2.up.y + c1.forward.z * c2.forward.y,
                                 c1.forward.x * c2.right.z + c1.forward.y * c2.up.z + c1.forward.z * c2.forward.z);
        return rc;
    }

    /// <summary>
    /// Multiplies a Quaternion and a UCoord3 together.
    /// </summary>
    /// <param name="q">The Quaternion.</param>
    /// <param name="c">The UCoord3.</param>
    /// <returns>The result of the multiplication.</returns>
    public static Quaternion operator *(Quaternion q, UCoord3 c)
    {
        Quaternion q1 = c.ToQuaternion();
        return q * q1;
    }

    /// <summary>
    /// Multiplies a UCoord3 and a Quaternion together.
    /// </summary>
    /// <param name="c">The UCoord3.</param>
    /// <param name="q">The Quaternion.</param>
    /// <returns>The result of the multiplication.</returns>
    public static UCoord3 operator *(UCoord3 c, Quaternion q)
    {
        UCoord3 rc = new UCoord3();
        rc.right = q * c.right;
        rc.up = q * c.up;
        rc.forward = q * c.forward;
        return rc;
    }

    /// <summary>
    /// Divides two UCoord3 instances.
    /// </summary>
    /// <param name="c1">The numerator UCoord3.</param>
    /// <param name="c2">The denominator UCoord3.</param>
    /// <returns>The result of the division.</returns>
    public static UCoord3 operator /(UCoord3 c1, UCoord3 c2)
    {
        UCoord3 rc = new UCoord3();
        rc.right = new Vector3(c1.right.x / c2.right.x + c1.right.y / c2.up.x + c1.right.z / c2.forward.x,
                               c1.right.x / c2.right.y + c1.right.y / c2.up.y + c1.right.z / c2.forward.y,
                               c1.right.x / c2.right.z + c1.right.y / c2.up.z + c1.right.z / c2.forward.z);
        rc.up = new Vector3(c1.up.x / c2.right.x + c1.up.y / c2.up.x + c1.up.z / c2.forward.x,
                            c1.up.x / c2.right.y + c1.up.y / c2.up.y + c1.up.z / c2.forward.y,
                            c1.up.x / c2.right.z + c1.up.y / c2.up.z + c1.up.z / c2.forward.z);
        rc.forward = new Vector3(c1.forward.x / c2.right.x + c1.forward.y / c2.up.x + c1.forward.z / c2.forward.x,
                                 c1.forward.x / c2.right.y + c1.forward.y / c2.up.y + c1.forward.z / c2.forward.y,
                                 c1.forward.x / c2.right.z + c1.forward.y / c2.up.z + c1.forward.z / c2.forward.z);
        return rc;
    }

    /// <summary>
    /// Divides a vector by a UCoord3.
    /// </summary>
    /// <param name="v">The vector to divide.</param>
    /// <param name="c">The UCoord3.</param>
    /// <returns>The result of the division.</returns>
    public static Vector3 operator /(Vector3 v, UCoord3 c)
    {
        return new Vector3(Vector3.Dot(v, c.right), Vector3.Dot(v, c.up), Vector3.Dot(v, c.forward));
    }

    /// <summary>
    /// Calculates the reverse of the UCoord3 instance.
    /// </summary>
    /// <returns>The reverse of the UCoord3 instance.</returns>
    public UCoord3 Reversed()
    {
        return one / this;
    }

    /// <summary>
    /// Converts the UCoord3 instance to Euler angles.
    /// </summary>
    /// <returns>The Euler angles representation of the UCoord3 instance.</returns>
    public Vector3 ToEulerAngles()
    {
        Quaternion q = Quaternion.LookRotation(forward, up);
        Vector3 eulerAngles = q.eulerAngles;
        float pitch = eulerAngles.x;
        float yaw = eulerAngles.y;
        float roll = eulerAngles.z;
        return new Vector3(pitch, yaw, roll);
    }
}

/// <summary>
/// ***************************************************************************************
/// Represents a three-dimensional coordinate system with scaling and translation in Unity.
/// ***************************************************************************************
/// </summary>
[System.Serializable]
public struct Coord3
{
    public static readonly Coord3 zero = new Coord3(UCoord3.one, Vector3.one, Vector3.zero);
    public static readonly Coord3 one = new Coord3(UCoord3.one, Vector3.one, Vector3.zero);

    public Vector3 origin;
    public Vector3 scale;
    public UCoord3 uCoord;

    /// <summary>
    /// Creates a new Coord3 instance with the specified UCoord3, scale, and origin.
    /// </summary>
    /// <param name="uCoord">The UCoord3 instance.</param>
    /// <param name="scale">The scale vector.</param>
    /// <param name="origin">The origin vector.</param>
    public Coord3(UCoord3 uCoord, Vector3 scale, Vector3 origin)
    {
        this.uCoord = uCoord;
        this.scale = scale;
        this.origin = origin;
    }
    /// <summary>
    /// Creates a new Coord3 instance with the specified UCoord3 and origin.
    /// </summary>
    /// <param name="uCoord">The UCoord3 instance.</param>
    /// <param name="origin">The origin vector.</param>
    public Coord3(UCoord3 uCoord, Vector3 origin)
    {
        this.uCoord = uCoord;
        this.scale = Vector3.one;
        this.origin = origin;
    }

    /// <summary>
    /// Creates a new Coord3 instance with the specified UCoord3.
    /// </summary>
    /// <param name="uCoord">The UCoord3 instance.</param>
    public Coord3(UCoord3 uCoord)
    {
        this.uCoord = uCoord;
        this.scale = Vector3.one;
        this.origin = Vector3.zero;
    }

    /// <summary>
    /// Creates a new Coord3 instance with the specified origin.
    /// </summary>
    /// <param name="origin">The origin vector.</param>
    public Coord3(Vector3 origin)
    {
        this.uCoord = UCoord3.one;
        this.scale = Vector3.one;
        this.origin = origin;
    }

    /// <summary>
    /// Creates a new Coord3 instance with the specified x-axis, y-axis, and z-axis vectors.
    /// </summary>
    /// <param name="ux">The x-axis vector.</param>
    /// <param name="uy">The y-axis vector.</param>
    /// <param name="uz">The z-axis vector.</param>
    public Coord3(Vector3 ux, Vector3 uy, Vector3 uz)
    {
        this.uCoord = new UCoord3(ux, uy, uz);
        this.scale = Vector3.one;
        this.origin = Vector3.zero;
    }

    /// <summary>
    /// Creates a new Coord3 instance with the specified x-axis and y-axis vectors.
    /// </summary>
    /// <param name="ux">The x-axis vector.</param>
    /// <param name="uy">The y-axis vector.</param>
    public Coord3(Vector3 ux, Vector3 uy)
    {
        this.uCoord = new UCoord3(ux, uy);
        this.scale = Vector3.one;
        this.origin = Vector3.zero;
    }

    /// <summary>
    /// Creates a new Coord3 instance with the specified x-axis, y-axis, z-axis, and origin vectors.
    /// </summary>
    /// <param name="ux">The x-axis vector.</param>
    /// <param name="uy">The y-axis vector.</param>
    /// <param name="uz">The z-axis vector.</param>
    /// <param name="origin">The origin vector.</param>
    public Coord3(Vector3 ux, Vector3 uy, Vector3 uz, Vector3 origin)
    {
        this.uCoord = new UCoord3(ux, uy, uz);
        this.scale = Vector3.one;
        this.origin = origin;
    }

    /// <summary>
    /// Creates a new Coord3 instance with the specified x-axis, y-axis, z-axis, scale, and origin vectors.
    /// </summary>
    /// <param name="ux">The x-axis vector.</param>
    /// <param name="uy">The y-axis vector.</param>
    /// <param name="uz">The z-axis vector.</param>
    /// <param name="scale">The scale vector.</param>
    /// <param name="origin">The origin vector.</param>
    public Coord3(Vector3 ux, Vector3 uy, Vector3 uz, Vector3 scale, Vector3 origin)
    {
        this.uCoord = new UCoord3(ux, uy, uz);
        this.scale = scale;
        this.origin = origin;
    }

    /// <summary>
    /// Transforms the specified vector using the Coord3 instance.
    /// </summary>
    /// <param name="p">The vector to transform.</param>
    /// <param name="c">The Coord3 instance.</param>
    /// <returns>The transformed vector.</returns>
    public static Vector3 operator *(Vector3 p, Coord3 c)
    {
        return p.x * c.uCoord.right + p.y * c.uCoord.up + p.z * c.uCoord.forward + c.origin;
    }

    /// <summary>
    /// Multiplies two Coord3 instances together.
    /// </summary>
    /// <param name="c1">The first Coord3.</param>
    /// <param name="c2">The second Coord3.</param>
    /// <returns>The result of the multiplication.</returns>
    public static Coord3 operator *(Coord3 c1, Coord3 c2)
    {
        Coord3 rc = new Coord3();
        rc.uCoord = c1.uCoord * c2.uCoord;
        rc.scale = new Vector3(c1.scale.x * c2.scale.x, c1.scale.y * c2.scale.y, c1.scale.z * c2.scale.z);
        rc.origin = c2.origin + c1.origin.x * c2.scale.x * c2.uCoord.right + c1.origin.y * c2.scale.y * c2.uCoord.up + c1.origin.z * c2.scale.z * c2.uCoord.forward;
        return rc;
    }

    /// <summary>
    /// Multiplies a Quaternion and a Coord3 together.
    /// </summary>
    /// <param name="q">The Quaternion.</param>
    /// <param name="c">The Coord3.</param>
    /// <returns>The result of the multiplication.</returns>
    public static Coord3 operator *(Quaternion q, Coord3 c)
    {
        Coord3 rc = c;
        rc.uCoord = rc.uCoord * q;
        rc.origin = q * rc.origin;
        return rc;
    }

    /// <summary>
    /// Multiplies a Coord3 and a Quaternion together.
    /// </summary>
    /// <param name="c">The Coord3.</param>
    /// <param name="q">The Quaternion.</param>
    /// <returns>The result of the multiplication.</returns>
    public static Coord3 operator *(Coord3 c, Quaternion q)
    {
        return new Coord3(c.uCoord * q);
    }

    /// <summary>
    /// Divides two Coord3 instances.
    /// </summary>
    /// <param name="c1">The numerator Coord3.</param>
    /// <param name="c2">The denominator Coord3.</param>
    /// <returns>The result of the division.</returns>
    public static Coord3 operator /(Coord3 c1, Coord3 c2)
    {
        Coord3 rc = new Coord3();
        rc.uCoord = c1.uCoord / c2.uCoord;
        rc.scale = new Vector3(c1.scale.x / c2.scale.x, c1.scale.y / c2.scale.y, c1.scale.z / c2.scale.z);
        rc.origin = c1.origin / c2;
        return rc;
    }

    /// <summary>
    /// Divides a Coord3 by a Quaternion.
    /// </summary>
    /// <param name="c">The Coord3.</param>
    /// <param name="q">The Quaternion.</param>
    /// <returns>The result of the division.</returns>
    public static Coord3 operator /(Coord3 c, Quaternion q)
    {
        return c * Quaternion.Inverse(q);
    }

    /// <summary>
    /// Divides a vector by a Coord3.
    /// </summary>
    /// <param name="p">The vector to divide.</param>
    /// <param name="c">The Coord3.</param>
    /// <returns>The result of the division.</returns>
    public static Vector3 operator /(Vector3 p, Coord3 c)
    {
        Vector3 v = p - c.origin;
        return new Vector3(Vector3.Dot(v, c.uCoord.right) / c.scale.x, Vector3.Dot(v, c.uCoord.up) / c.scale.y, Vector3.Dot(v, c.uCoord.forward) / c.scale.z);
    }

    /// <summary>
    /// Calculates the reverse of the Coord3 instance.
    /// </summary>
    /// <returns>The reverse of the Coord3 instance.</returns>
    public Coord3 Reversed()
    {
        return one / this;
    }
}
