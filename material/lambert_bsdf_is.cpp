#include "lambert_bsdf_is.h"

LambertBSDF_IS::LambertBSDF_IS(Color c) : LambertBSDF(c) {}

double LambertBSDF_IS::pdf(const Vector3d& wi,
                           const Vector3d& wo,
                           const Vector3d& normal) const {
     return normal.dot(wi) * M_1_PI;
};

Vector3d LambertBSDF_IS::sample(const Vector3d& wo,
                                const Vector3d& normal,
                                Random& rand) const  {
     double r = sqrt(rand.unit_rand());
     double theta = 2 * M_PI * rand.unit_rand();
     double x = r * cos(theta);
     double y = r * sin(theta);
     double z = sqrtf(1.0 - r * r);

     return this->to_world(x, y, z, normal).normalized();
};
