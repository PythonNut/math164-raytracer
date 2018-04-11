#pragma once
#include "uniform_sampling_bsdf.h"

class ConductorBSDF : public UniformSamplingBSDF {
protected:
     Color fresnel(double cos_theta_i,
                   Color eta,
                   Color kappa) const;
};
