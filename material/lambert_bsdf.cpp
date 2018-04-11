#include "lambert_bsdf.h"

LambertBSDF::LambertBSDF(Color c) : color(c) {}

Color LambertBSDF::eval(const Vector3d& wi,
                        const Vector3d& wo,
                        const Vector3d& normal) const {
     return this->color * M_1_PI * normal.dot(wi);
};

Color LambertBSDF::emission(const Vector3d& wi,
                            const Vector3d& wo,
                            const Vector3d& normal) const {
     return Color(0,0,0);
};
