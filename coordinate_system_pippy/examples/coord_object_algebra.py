"""Basic coord3 object algebra: transform, inverse transform, composition."""

import math

from coordinate_system import coord3, quat, vec3


C = coord3(
    vec3(10.0, 0.0, 0.0),
    quat.from_axis_angle(vec3.UZ, math.pi / 2.0),
    vec3(2.0, 2.0, 2.0),
)

local = vec3(1.0, 0.0, 0.0)
world = C.to_world(local)
round_trip = C.to_local(world)

print("C:", C)
print("local:", local)
print("world:", world)
print("round trip:", round_trip)
