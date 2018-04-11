#include "linear_algebra.h"

Vec::Vec(double x, double y, double z) : a(x), b(y), c(z) { }

Vec Vec::operator+(const Vec &o) const {
     return Vec(a + o.a, b + o.b, c + o.c);
}

Vec Vec::operator-(const Vec &o) const {
     return Vec(a - o.a, b - o.b, c - o.c);
}

Vec Vec::operator-() const {
     return Vec(-a, -b, -c);
}

Vec Vec::operator*(double o) const {
     return Vec(a*o, b*o, c*o);
}

Vec Vec::operator*(const Vec &o) const {
     return Vec(a*o.a, b*o.b, c*o.c);
}

Vec Vec::operator/(double o) const {
     return Vec(a/o, b/o, c/o);
}

Vec Vec::operator/(const Vec &o) const {
     return Vec(a/o.a, b/o.b, c/o.c);
}

Vec& Vec::operator*=(double o) {
     this->a *= o;
     this->b *= o;
     this->c *= o;
     return *this;
}

Vec& Vec::operator+=(const Vec& o ) {
     this->a += o.x();
     this->b += o.y();
     this->c += o.z();
     return *this;
}

Vec Vec::normalized () const {
     return this->operator/(sqrt(a*a + b*b + c*c));
}

double Vec::dot(const Vec &o) const {
     return a*o.a + b*o.b + c*o.c;
}

Vec Vec::cross(const Vec& o) const{
     return Vec(b*o.c - c*o.b, c*o.a - a*o.c, a*o.b - b*o.a);
}

double Vec::x() const {return this->a;}
double Vec::y() const {return this->b;}
double Vec::z() const {return this->c;}

double Vec::maxCoeff() const{
     return max(this->a, max(this->b, this->c));
}

Vec operator* (double k, const Vec& v) {
     return v * k;
}
