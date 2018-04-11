#include "ggx_conductor_bsdf_is.h"

GGXConductorBSDF_IS::GGXConductorBSDF_IS(Color c, Color n, Color k, double a) :
     GGXConductorBSDF(c, n, k, a) {  }


double GGXConductorBSDF_IS::pdf(const Vector3d& wi,
                                const Vector3d& wo,
                                const Vector3d& normal) const {
     Vector3d m = (wi + wo).normalized();

     double p_m = microfacet_distribution(m, normal) * abs(normal.dot(m));
     return p_m/abs(4 * wi.dot(m) * wi.dot(normal));
};


Vector3d GGXConductorBSDF_IS::sample(const Vector3d& wo,
                                     const Vector3d& normal,
                                     Random& rand) const {
     double xi1 = rand.unit_rand();
     double xi2 = rand.unit_rand();

     double alpha2 = this->alpha * this->alpha;
     double phi = 2 * M_PI * xi1;
     double cos_theta = sqrt((1 - xi2)/(xi2 * (alpha2 - 1) + 1));
     double sin_theta = sqrt(1 - cos_theta * cos_theta);

     double x = sin_theta * cos(phi);
     double y = sin_theta * sin(phi);
     double z = cos_theta;

     Vector3d m = this->to_world(x, y, z, normal).normalized();
     Vector3d wi = 2 * abs(wo.dot(m)) * m - wo;

     return wi;
};
