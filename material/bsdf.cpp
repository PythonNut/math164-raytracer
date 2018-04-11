#include "bsdf.h"

Vector3d BSDF::to_world(double x,
                        double y,
                        double z,
                        const Vector3d& normal) const {
     Vector3d majorAxis;

     const double inv_sqrt3 = 1/sqrt(3);
     if (abs(normal.x()) < inv_sqrt3) {
          majorAxis = Vector3d(1, 0, 0);
     } else if (abs(normal.y()) < inv_sqrt3) {
          majorAxis = Vector3d(0, 1, 0);
     } else {
          majorAxis = Vector3d(0, 0, 1);
     }

     const Vector3d u = normal.cross(majorAxis).normalized();
          const Vector3d v = normal.cross(u);
          const Vector3d w = normal;

          return u * x + v * y + w * z;
}
