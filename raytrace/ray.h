#pragma once
#include "../util/linear_algebra.h"

struct Ray {
     Vector3d origin, direction;
     Ray(Vector3d orig, Vector3d dir)
          : origin(orig)
          , direction(dir) { }
};
