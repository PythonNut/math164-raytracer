#pragma once
#include <random>
#include "linear_algebra.h"

using namespace std;

double clamp_intensity(const double x);

Color clamp_intensity(Color c);

template <typename T> int sgn(T val);

int toPPM(const double x);

class Random {
     default_random_engine rand_engine;
     uniform_real_distribution<> uniform_rand;
public:
     Random();

     double unit_rand();
};
