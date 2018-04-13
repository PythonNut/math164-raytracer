#pragma once
#include "pixel_sampler.h"

class MarsagliaPolarSampler : public PixelSampler {
public:
     MarsagliaPolarSampler(Random& r);

     virtual pair<double, double> sample() override;
};
