#pragma once
#include "../material/lambert_bsdf.h"
#include "../material/lambert_bsdf_is.h"
#include "../material/ggx_conductor_bsdf.h"
#include "../material/ggx_conductor_bsdf_is.h"
#include "../material/emission_bsdf.h"
#include "../geometry/sphere.h"
#include "../raytrace/object.h"

typedef LambertBSDF_IS Diffuse;
typedef GGXConductorBSDF_IS Glossy;

extern const vector<Object> scene;

extern const Ray camera;
