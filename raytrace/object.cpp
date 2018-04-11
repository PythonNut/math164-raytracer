#include "object.h"

const optional<pair<double, vector<Object>::const_iterator>> intersect(const Ray& ray,
                                                                       const vector<Object>& scene) {
     double nearest_distance = std::numeric_limits<double>::infinity();
     auto nearest_object = scene.end();

     for (auto it=scene.begin(); it != scene.end(); ++it) {
          auto hit_check = it->shape->intersect(ray);
          if (hit_check.has_value()) {
               double distance = hit_check.value();
               if (distance < nearest_distance) {
                    nearest_object = it;
                    nearest_distance = distance;
               }
          }
     }

     if (nearest_object == scene.end()) {
          return {};
     }

     return make_pair(nearest_distance, nearest_object);
}
