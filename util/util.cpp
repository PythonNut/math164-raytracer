#include "util.h"

double clamp_intensity(const double x) {
     return min(max(0.0, x), 1.0);
}

Color clamp_intensity(Color c) {
     return Color(clamp_intensity(c.x()),
                  clamp_intensity(c.y()),
                  clamp_intensity(c.z()));
}

template <typename T> int sgn(T val) {
     return (T(0) < val) - (val < T(0));
}


int toPPM(const double x) {
     return pow(clamp_intensity(x), 1/2.2) * 255 + 0.5;
}

Random::Random() : uniform_rand(0, 1) {
     random_device rand_dev;
     this->rand_engine = default_random_engine(rand_dev());
}

double Random::unit_rand() {
     return this->uniform_rand(this->rand_engine);
}
