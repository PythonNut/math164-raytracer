#pragma once
#include <random>
#include <array>
#include "linear_algebra.h"

using namespace std;

double clamp_intensity(const double x);

Color clamp_intensity(Color c);

template <typename T> int sgn(T val);

int toPPM(const double x);

class Xoroshiro128 {
private:
     array<uint64_t, 2> state;

public:
     using result_type = uint64_t;

     Xoroshiro128();

     result_type max();
     result_type min();
     result_type operator()();
};

class Random {
     Xoroshiro128 rand_engine;
     uniform_real_distribution<> uniform_rand;
public:
     Random();

     double unit_rand();
};
