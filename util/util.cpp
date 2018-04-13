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

Xoroshiro128::Xoroshiro128() {
     // http://www.iro.umontreal.ca/~lecuyer/myftp/papers/lfsr04.pdf
     mt19937_64 rng(random_device{}());
     rng.discard(700000);
     this->state[0] = rng();
     this->state[1] = rng();
}

Xoroshiro128::result_type Xoroshiro128::max() {
     return numeric_limits<Xoroshiro128::result_type>::max();
}

Xoroshiro128::result_type Xoroshiro128::min() {
     return numeric_limits<Xoroshiro128::result_type>::min();
}

Xoroshiro128::result_type Xoroshiro128::operator()() {
     uint64_t s0 = this->state[0];
     uint64_t s1 = this->state[1];
     uint64_t result = s0 + s1;
     s1 ^= s0;
     this->state[0] = ((s0 << 55) | (s0 >> 9)) ^ s1 ^ (s1 << 14);
     this->state[1] = (s1 << 36) | (s1 >> 28);
     return result;
}

Random::Random() : rand_engine(), uniform_rand(0, 1) { }

double Random::unit_rand() {
     return this->uniform_rand(this->rand_engine);
}
