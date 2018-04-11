#include "radiance.h"

const Color radiance(const vector<Object>& scene,
                     const Ray& ray,
                     int depth,
                     Random& rand) {
     auto hit_check = intersect(ray, scene);
     if (!hit_check.has_value() || depth > 5) {
          return Color(0, 0, 0);
     }

     double distance;
     vector<Object>::const_iterator obj;
     tie(distance, obj) = hit_check.value();

     Vector3d intersect_point = ray.origin + ray.direction * distance;

     Vector3d wo = -ray.direction.normalized();
     Vector3d normal = obj->shape->get_normal(intersect_point);
     Vector3d oriented_normal = normal;
     if (normal.dot(wo) < 0) {
          oriented_normal *= -1;
     }

     Vector3d wi = obj->mat->sample(wo, oriented_normal, rand);
     Color r = radiance(scene, Ray(intersect_point, wi), depth+1, rand);
     double pdf = obj->mat->pdf(wi, wo, oriented_normal);

     Color emit = obj->mat->emission(wi, wo, oriented_normal);

     if (emit.maxCoeff()) {
          return emit;
     }

     return obj->mat->eval(wi, wo, oriented_normal) * r/pdf;
}
