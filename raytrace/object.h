#pragma once
#include <optional>
#include <memory>
#include "../geometry/geometry.h"
#include "../material/bsdf.h"

struct Object {
     unique_ptr<Geometry> shape;
     unique_ptr<BSDF> mat;

     Object(unique_ptr<Geometry>(s), unique_ptr<BSDF>(m))
          : shape(move(s))
          , mat(move(m)) { }
};


const optional<pair<double, vector<Object>::const_iterator>> intersect(
     const Ray& ray,
     const vector<Object>& scene);
