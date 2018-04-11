#pragma once
#include "geometry.h"

class Sphere: public Geometry {
private:
     double radius;
     Vector3d origin;

     double get_radius() const;

     Vector3d get_origin() const;

public:
     Sphere(double r, Vector3d orig);

     virtual optional<double> intersect(const Ray& ray) const override;

     virtual Vector3d get_normal(const Vector3d& point) const override;
};
