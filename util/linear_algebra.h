#pragma once
#include <eigen3/Eigen/Dense>

/* Comment this line out to get a small speed boost, at the cost of
 * some features */
#define USE_EIGEN
using namespace std;

class Vec {
private:
     double a, b, c;

public:
     Vec(double x=0, double y=0, double z=0);

     Vec operator+(const Vec &o) const;
     Vec operator-(const Vec &o) const;
     Vec operator-() const;
     Vec operator*(double o) const;
     Vec operator*(const Vec &o) const;
     Vec operator/(double o) const;
     Vec operator/(const Vec &o) const;
     Vec& operator*=(double o);
     Vec& operator+=(const Vec& o );

     Vec normalized () const;
     double dot(const Vec &o) const;
     Vec cross(const Vec& o) const;

     double x() const;
     double y() const;
     double z() const;

     double maxCoeff() const;
};

Vec operator* (double k, const Vec& v);

#ifdef USE_EIGEN
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Array3d Color;
#else
typedef Vec Vector3d;
typedef Vec Color;
#endif
