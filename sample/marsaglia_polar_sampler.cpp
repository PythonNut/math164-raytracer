#include "marsaglia_polar_sampler.h"

MarsagliaPolarSampler::MarsagliaPolarSampler(Random& r) : PixelSampler(r) {  };

pair<double, double> MarsagliaPolarSampler::sample() {
     double v1, v2, s;
     do {
          v1 = 2.0 * this->rand.unit_rand() - 1.0;
          v2 = 2.0 * this->rand.unit_rand() - 1.0;
          s = v1 * v1 + v2 * v2;
     } while (s >= 1.0 || s == 0);

     s = sqrt((-2.0 * log(s)) / s);

     return make_pair(v1 * s, v2 * s);
};
