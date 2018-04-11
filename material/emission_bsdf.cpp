#include "emission_bsdf.h"

EmissionBSDF::EmissionBSDF(Color c) : LambertBSDF_IS(c) {}


Color EmissionBSDF::emission(const Vector3d& wi,
                             const Vector3d& wo,
                             const Vector3d& normal) const {
     return this->color;
};
