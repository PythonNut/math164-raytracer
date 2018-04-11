#pragma once
#include <optional>
#include "../util/linear_algebra.h"
#include "../raytrace/ray.h"

class Geometry {
public:
     virtual Vector3d get_normal(const Vector3d& point) const = 0;
     virtual optional<double> intersect(const Ray& ray) const = 0;
};
