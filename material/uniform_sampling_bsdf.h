#pragma once
#include "bsdf.h"

class UniformSamplingBSDF : public BSDF {
     virtual double pdf(const Vector3d& wi,
                        const Vector3d& wo,
                        const Vector3d& normal) const override;

     virtual Vector3d sample(const Vector3d& wo,
                             const Vector3d& normal,
                             Random& rand) const override;
};
