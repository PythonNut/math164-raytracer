#pragma once
#include "ggx_conductor_bsdf.h"

class GGXConductorBSDF_IS : public GGXConductorBSDF {
public:
     GGXConductorBSDF_IS(Color c, Color n, Color k, double a);

     virtual double pdf(const Vector3d& wi,
                        const Vector3d& wo,
                        const Vector3d& normal) const override;

     virtual Vector3d sample(const Vector3d& wo,
                             const Vector3d& normal,
                             Random& rand) const override;
};
