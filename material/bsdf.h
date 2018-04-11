#pragma once
#include "../util/linear_algebra.h"
#include "../util/util.h"

class BSDF {
protected:
     /* to_world, convert a vector from the normal coordinate system to
      * the global coordinate system. */
     Vector3d to_world(double x,
                       double y,
                       double z,
                       const Vector3d& normal) const;

public:
     /* Sampling probability density function. */
     virtual double pdf(const Vector3d& wi,
                        const Vector3d& wo,
                        const Vector3d& normal) const = 0;

     /* Generate a reflected ray with the distribution dictated by the
      * pdf() */
     virtual Vector3d sample(const Vector3d& wo,
                             const Vector3d& normal,
                             Random& rand) const = 0;

     /* Calculate the BSDF, as a function of the incidence angles. */
     virtual Color eval(const Vector3d& wi,
                        const Vector3d& wo,
                        const Vector3d& normal) const = 0;

     /* Calculate the emission as a function of the incidence
      * angles. */
     virtual Color emission(const Vector3d& wi,
                            const Vector3d& wo,
                            const Vector3d& normal) const = 0;
};
