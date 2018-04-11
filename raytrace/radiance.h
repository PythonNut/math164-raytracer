#pragma once
#include "ray.h"
#include "object.h"

const Color radiance(const vector<Object>& scene,
                     const Ray& ray,
                     int depth,
                     Random& rand);
