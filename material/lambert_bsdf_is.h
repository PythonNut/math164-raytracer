#pragma once
#include "lambert_bsdf.h"

class LambertBSDF_IS : public LambertBSDF {
public:
     LambertBSDF_IS(Color c);

     virtual double pdf(const Vector3d& wi,
                        const Vector3d& wo,
                        const Vector3d& normal) const override;

     virtual Vector3d sample(const Vector3d& wo,
                             const Vector3d& normal,
                             Random& rand) const override;
};
