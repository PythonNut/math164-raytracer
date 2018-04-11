#include "ggx_conductor_bsdf.h"

GGXConductorBSDF::GGXConductorBSDF(Color c, Color n, Color k, double a) :
     color(c), eta(n), kappa(k), alpha(a) { }


double GGXConductorBSDF::microfacet_distribution(const Vector3d& m,
                                                 const Vector3d& n) const {
     double n_dot_m = n.dot(m);
     double sin_theta_m = sqrt(1 - n_dot_m * n_dot_m);
     double tan_theta_m = sin_theta_m / n_dot_m;
     if (n_dot_m <= 0) {
          return 0;
     }

     double alpha2 = this->alpha * this->alpha;

     double denom = n_dot_m * n_dot_m * (alpha2 + tan_theta_m * tan_theta_m);

     return alpha2/(M_PI * denom * denom);
}


double GGXConductorBSDF::shadow_masking(const Vector3d& v,
                                        const Vector3d& m,
                                        const Vector3d& n) const {
     double alpha2 = this->alpha * this->alpha;
     double n_dot_v=n.dot(v);
     if (v.dot(m)/n_dot_v <= 0) {
          return 0;
     }

     double tan_theta_v = tan(acos(n_dot_v));

     return 2/(1 + sqrt(1 + alpha2 * tan_theta_v * tan_theta_v));
}


Color GGXConductorBSDF::eval(const Vector3d& wi,
                             const Vector3d& wo,
                             const Vector3d& normal) const {
     Vector3d h = (wi + wo).normalized();
     // Paper implies this condition, but the results were bad.
     // if (wi.dot(wo) < 0) {
     //     h *= -1;
     // }

     Color F = this->fresnel(h.dot(wi), this->eta, this->kappa);
     double D = this->microfacet_distribution(h, normal);
     double G = this->shadow_masking(wi, h, normal) * this->shadow_masking(wo, h, normal);

     return this->color * F * D * G/ (4 * abs(normal.dot(wi)) * abs(normal.dot(wo)));
};


Color GGXConductorBSDF::emission(const Vector3d& wi,
                                 const Vector3d& wo,
                                 const Vector3d& normal) const {
     return Color(0,0,0);
};
