#include "sphere.h"

Sphere::Sphere(double r, Vector3d orig)
     : radius(r)
     , origin(orig) { }


double Sphere::get_radius() const {
     return this->radius;
}


Vector3d Sphere::get_origin() const {
     return this->origin;
}


optional<double> Sphere::intersect(const Ray& ray) const {
     Vector3d op = this->origin - ray.origin;
     double epsilon = 1e-4;
     double b = op.dot(ray.direction);
     double r2 = this->radius * this->radius;
     double descriminant = b*b - op.dot(op) + r2;

     if (descriminant < 0) {
          return {};
     }

     descriminant = sqrt(descriminant);

     if (b - descriminant > epsilon) {
          return b - descriminant;
     }
     if (b + descriminant > epsilon) {
          return b + descriminant;
     }

     return {};
}


Vector3d Sphere::get_normal(const Vector3d& point) const {
     return (point - this->get_origin()).normalized();
}
