#pragma once
#include "lambert_bsdf_is.h"

class EmissionBSDF : public LambertBSDF_IS {
public:
     EmissionBSDF(Color c);

     virtual Color emission(const Vector3d& wi,
                            const Vector3d& wo,
                            const Vector3d& normal) const override;
};
