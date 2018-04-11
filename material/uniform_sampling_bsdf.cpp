#include "uniform_sampling_bsdf.h"

double UniformSamplingBSDF::pdf(const Vector3d& wi,
                                const Vector3d& wo,
                                const Vector3d& normal) const {
     return M_1_PI/2;
};

Vector3d UniformSamplingBSDF::sample(const Vector3d& wo,
                                     const Vector3d& normal,
                                     Random& rand) const {
     double theta = 2 * M_PI * rand.unit_rand();
     double phi = acos(1 - 2 * rand.unit_rand());
     double x = sin(phi) * cos(theta);
     double y = sin(phi) * sin(theta);
     double z = abs(cos(phi));

     return this->to_world(x, y, z, normal).normalized();
};
