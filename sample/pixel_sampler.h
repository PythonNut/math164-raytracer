#pragma once
#include "../util/util.h"

class PixelSampler {
protected:
     Random& rand;

public:
     PixelSampler(Random& rand);

     virtual pair<double, double> sample() = 0;
};
