#include "box_muller_sampler.h"

BoxMullerSampler::BoxMullerSampler(Random& r) : PixelSampler(r) {  };

pair<double, double> BoxMullerSampler::sample() {
     double u1 = this->rand.unit_rand();
     double u2 = this->rand.unit_rand();
     double r = sqrt(-2*log(u1));
     double theta = 2 * M_PI * u2;
     return make_pair(r*cos(theta), r*sin(theta));
};
