#pragma once
#include "pixel_sampler.h"

class BoxMullerSampler : public PixelSampler {
public:
     BoxMullerSampler(Random& r);

     virtual pair<double, double> sample() override;
};
