#pragma once
#include "uniform_sampling_bsdf.h"

class LambertBSDF : public UniformSamplingBSDF {
protected:
     Color color;
public:
     LambertBSDF(Color c);

     virtual Color eval(const Vector3d& wi,
                        const Vector3d& wo,
                        const Vector3d& normal) const override;

     virtual Color emission(const Vector3d& wi,
                            const Vector3d& wo,
                            const Vector3d& normal) const override;
};
