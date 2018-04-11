#pragma once
#include "conductor_bsdf.h"

class GGXConductorBSDF : public ConductorBSDF {
protected:
     Color color, eta, kappa;
     double alpha;

     double microfacet_distribution(const Vector3d& m,
                                    const Vector3d& n) const;

     double shadow_masking(const Vector3d& v,
                           const Vector3d& m,
                           const Vector3d& n) const;

public:
     GGXConductorBSDF(Color c, Color n, Color k, double a);

     virtual Color eval(const Vector3d& wi,
                        const Vector3d& wo,
                        const Vector3d& normal) const override;

     virtual Color emission(const Vector3d& wi,
                            const Vector3d& wo,
                            const Vector3d& normal) const override;
};
