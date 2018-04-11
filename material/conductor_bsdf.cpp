#include "conductor_bsdf.h"

Color ConductorBSDF::fresnel(double cos_theta_i,
                             Color eta,
                             Color kappa) const {
     Color eta2kappa2 = eta * eta + kappa * kappa;
     Color pm = 2 * eta * cos_theta_i;
     double cos_theta_i2 = cos_theta_i * cos_theta_i;

     Color r_para2_base = eta2kappa2 * cos_theta_i2 + 1;
     Color r_para2 = (r_para2_base - pm)/(r_para2_base + pm);

     Color r_perp2_base = eta2kappa2 + cos_theta_i2;
     Color r_perp2 = (r_perp2_base - pm)/(r_perp2_base + pm);

     return clamp_intensity((r_para2 + r_perp2)/2);
};
