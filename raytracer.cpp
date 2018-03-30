#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>
#include <random>
#include <functional>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

/* Comment this line out to get a small speed boost, at the cost of
 * some features */
#define USE_EIGEN
using namespace std;

class Vec {
private:
    double a, b, c;

public:
    Vec(double x=0, double y=0, double z=0) : a(x), b(y), c(z) { }

    Vec operator+(const Vec &o) const {
        return Vec(a + o.a, b + o.b, c + o.c);
    }

    Vec operator-(const Vec &o) const {
        return Vec(a - o.a, b - o.b, c - o.c);
    }

    Vec operator*(double o) const {
        return Vec(a*o, b*o, c*o);
    }

    Vec operator*(const Vec &o) const {
        return Vec(a*o.a, b*o.b, c*o.c);
    }

    Vec operator/(double o) const {
        return Vec(a/o, b/o, c/o);
    }

    Vec operator/(const Vec &o) const {
        return Vec(a/o.a, b/o.b, c/o.c);
    }

    Vec& operator*=(double o) {
        this->a *= o;
        this->b *= o;
        this->c *= o;
        return *this;
    }

    Vec& operator+=(const Vec& o ) {
        this->a += o.x();
        this->b += o.y();
        this->c += o.z();
        return *this;
    }

    Vec normalized () const {
        return this->operator/(sqrt(a*a + b*b + c*c));
    }

    double dot(const Vec &o) const {
        return a*o.a + b*o.b + c*o.c;
    }

    Vec cross(const Vec& o) const{
        return Vec(b*o.c - c*o.b, c*o.a - a*o.c, a*o.b - b*o.a);
    }

    double x() const {return this->a;}
    double y() const {return this->b;}
    double z() const {return this->c;}

    double maxCoeff() const{
        return max(this->a, max(this->b, this->c));
    }
};

#ifdef USE_EIGEN
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Array3d Color;
#else
typedef Vec Vector3d;
typedef Vec Color;
#endif

class Random {
    default_random_engine rand_engine;
    uniform_real_distribution<> uniform_rand;
public:
    Random() : uniform_rand(0, 1) {
        random_device rand_dev;
        this->rand_engine = default_random_engine(rand_dev());
    }

    double unit_rand() {
        return this->uniform_rand(this->rand_engine);
    }
};

struct Ray {
    Vector3d origin, direction;
    Ray(Vector3d orig, Vector3d dir)
        : origin(orig)
        , direction(dir) { }
};

class BSDF {
protected:
    /* to_world, convert a vector from the normal coordinate system to
     * the global coordinate system. */
    Vector3d to_world(double x,
                      double y,
                      double z,
                      const Vector3d& normal) const {
        Vector3d majorAxis;

        const double inv_sqrt3 = 1/sqrt(3);
        if (abs(normal.x()) < inv_sqrt3) {
            majorAxis = Vector3d(1, 0, 0);
        } else if (abs(normal.y()) < inv_sqrt3) {
            majorAxis = Vector3d(0, 1, 0);
        } else {
            majorAxis = Vector3d(0, 0, 1);
        }

        const Vector3d u = normal.cross(majorAxis).normalized();
        const Vector3d v = normal.cross(u);
        const Vector3d w = normal;

        return u * x + v * y + w * z;
    }

public:
    /* Sampling probability density function. */
    virtual double pdf(const Vector3d& wi,
                       const Vector3d& wo,
                       const Vector3d& normal) const = 0;

    /* Generate a reflected ray with the distribution dictated by the
     * pdf() */
    virtual Vector3d sample(const Vector3d& wo,
                            const Vector3d& normal,
                            Random& rand) const = 0;

    /* Calculate the BSDF, as a function of the incidence angles. */
    virtual Color eval(const Vector3d& wi,
                       const Vector3d& wo,
                       const Vector3d& normal) const = 0;

    /* Calculate the emission as a function of the incidence
     * angles. */
    virtual Color emission(const Vector3d& wi,
                           const Vector3d& wo,
                           const Vector3d& normal) const = 0;
};

class LambertBSDF : public BSDF {
protected:
    Color color;
public:
    LambertBSDF(Color c) : color(c) {}

    virtual double pdf(const Vector3d& wi,
                       const Vector3d& wo,
                       const Vector3d& normal) const override {
        return M_1_PI/2;
    };

    virtual Vector3d sample(const Vector3d& wo,
                            const Vector3d& normal,
                            Random& rand) const override {
        double theta = 2 * M_PI * rand.unit_rand();
        double phi = acos(1 - 2 * rand.unit_rand());
        double x = sin(phi) * cos(theta);
        double y = sin(phi) * sin(theta);
        double z = abs(cos(phi));

        return this->to_world(x, y, z, normal).normalized();
    };

    virtual Color eval(const Vector3d& wi,
                       const Vector3d& wo,
                       const Vector3d& normal) const override {
        return this->color * M_1_PI * normal.dot(wi);
    };

    virtual Color emission(const Vector3d& wi,
                           const Vector3d& wo,
                           const Vector3d& normal) const override {
        return Color(0,0,0);
    };
};

class LambertBSDF_IS : public LambertBSDF {
public:
    LambertBSDF_IS(Color c) : LambertBSDF(c) {}

    virtual double pdf(const Vector3d& in,
                       const Vector3d& out,
                       const Vector3d& normal) const override {
        return normal.dot(in) * M_1_PI;
    };

    virtual Vector3d sample(const Vector3d& wo,
                            const Vector3d& normal,
                            Random& rand) const override {
        double r = sqrt(rand.unit_rand());
        double theta = 2 * M_PI * rand.unit_rand();
        double x = r * cos(theta);
        double y = r * sin(theta);
        double z = sqrtf(1.0 - r * r);

        return this->to_world(x, y, z, normal).normalized();
    };
};

class EmissionBSDF : public LambertBSDF_IS {
public:
    EmissionBSDF(Color c) : LambertBSDF_IS(c) {}

    virtual Color emission(const Vector3d& wi,
                           const Vector3d& wo,
                           const Vector3d& normal) const override {
        return this->color;
    };
};

class Geometry {
public:
    virtual Vector3d get_normal(const Vector3d& point) const = 0;
    virtual optional<double> intersect(const Ray& ray) const = 0;
};

class Sphere: public Geometry {
private:
    double radius;
    Vector3d origin;

    double get_radius() const {
        return this->radius;
    }

    Vector3d get_origin() const {
        return this->origin;
    }

public:
    Sphere(double r, Vector3d orig)
        : radius(r)
        , origin(orig)
        { }

    virtual optional<double> intersect(const Ray& ray) const override {
        Vector3d op = this->origin - ray.origin;
        double epsilon = 1e-4;
        double b = op.dot(ray.direction);
        double r2 = this->radius * this->radius;
        double descriminant = b*b - op.dot(op) + r2;

        if (descriminant < 0) {
            return {};
        }

        descriminant = sqrt(descriminant);

        if (b - descriminant > epsilon) {
            return b - descriminant;
        }
        if (b + descriminant > epsilon) {
            return b + descriminant;
        }

        return {};
    }

    virtual Vector3d get_normal(const Vector3d& point) const override {
        return (point - this->get_origin()).normalized();
    }
};

struct Object {
    unique_ptr<Geometry> shape;
    unique_ptr<BSDF> mat;

    Object(unique_ptr<Geometry>(s), unique_ptr<BSDF>(m))
        : shape(move(s))
        , mat(move(m)) { }
};

typedef LambertBSDF_IS Diffuse;

Object scene_temp[] = {
    // left
    Object(move(make_unique<Sphere>(Sphere(1e5, Vector3d(1e5+1, 40.8, 81.6)))),
           move(make_unique<Diffuse>(Diffuse(Color(0.75, 0.25, 0.25))))),

    // right
    Object(move(make_unique<Sphere>(Sphere(1e5, Vector3d(-1e5+99, 40.8, 81.6)))),
           move(make_unique<Diffuse>(Diffuse(Color(0.25, 0.25, 0.75))))),

    // back
    Object(move(make_unique<Sphere>(Sphere(1e5, Vector3d(50,40.8, 1e5)))),
           move(make_unique<Diffuse>(Diffuse(Color(0.75, 0.75, 0.75))))),

    // front
    Object(move(make_unique<Sphere>(Sphere(1e5, Vector3d(50, 40.8, -1e5+170)))),
           move(make_unique<Diffuse>(Diffuse(Color(0, 0, 0))))),

    // bottom
    Object(move(make_unique<Sphere>(Sphere(1e5, Vector3d(50, 1e5, 81.6)))),
           move(make_unique<Diffuse>(Diffuse(Color(0.75, 0.75, 0.75))))),

    // top
    Object(move(make_unique<Sphere>(Sphere(1e5, Vector3d(50, -1e5+81.6, 81.6)))),
           move(make_unique<Diffuse>(Diffuse(Color(0.75, 0.75, 0.75))))),

    // mirror
    Object(move(make_unique<Sphere>(Sphere(16.5, Vector3d(27, 16.5, 47)))),
           move(make_unique<Diffuse>(Diffuse(Color(1, 1, 1)*.999)))),

    // glass
    Object(move(make_unique<Sphere>(Sphere(16.5, Vector3d(73, 16.5, 78)))),
           move(make_unique<Diffuse>(Diffuse(Color(1, 1, 1)*.999)))),

    // light
    Object(move(make_unique<Sphere>(Sphere(600, Vector3d(50, 681.6-.27, 81.6)))),
           move(make_unique<EmissionBSDF>(EmissionBSDF(Color(12, 12, 12)))))
};

const vector<Object> scene{make_move_iterator(begin(scene_temp)), make_move_iterator(end(scene_temp))};


const optional<pair<double, vector<Object>::const_iterator>> intersect(const Ray& ray,
                                                                       const vector<Object>& scene) {
    double nearest_distance = std::numeric_limits<double>::infinity();
    auto nearest_object = scene.end();

    for (auto it=scene.begin(); it != scene.end(); ++it) {
        auto hit_check = it->shape->intersect(ray);
        if (hit_check.has_value()) {
            double distance = hit_check.value();
            if (distance < nearest_distance) {
                nearest_object = it;
                nearest_distance = distance;
            }
        }
    }

    if (nearest_object == scene.end()) {
        return {};
    }

    return make_pair(nearest_distance, nearest_object);
}

const Color radiance(const vector<Object>& scene,
                     const Ray& ray,
                     int depth,
                     Random& rand) {
    auto hit_check = intersect(ray, scene);
    if (!hit_check.has_value() || depth > 5) {
        return Color(0, 0, 0);
    }

    double distance;
    vector<Object>::const_iterator obj;
    tie(distance, obj) = hit_check.value();

    Vector3d intersect_point = ray.origin + ray.direction * distance;

    Vector3d wo = -ray.direction.normalized();
    Vector3d normal = obj->shape->get_normal(intersect_point);
    Vector3d oriented_normal = normal;
    if (normal.dot(wo) < 0) {
        oriented_normal *= -1;
    }

    Vector3d wi = obj->mat->sample(wo, oriented_normal, rand);
    Color r = radiance(scene, Ray(intersect_point, wi), depth+1, rand);
    double pdf = obj->mat->pdf(wi, wo, oriented_normal);

    Color emit = obj->mat->emission(wi, wo, oriented_normal);

    if (emit.maxCoeff()) {
        return emit;
    }

    return obj->mat->eval(wi, wo, oriented_normal) * r/pdf;
}

double clamp_intensity(const double x) {
    return min(max(0.0, x), 1.0);
}

Color clamp_intensity(Color c) {
    return Color(clamp_intensity(c.x()),
                 clamp_intensity(c.y()),
                 clamp_intensity(c.z()));
}

int toPPM(const double x) {
    return pow(clamp_intensity(x), 1/2.2) * 255 + 0.5;
}

sf::Color toPPM(Color c) {
    return sf::Color(toPPM(c.x()), toPPM(c.y()), toPPM(c.z()));
}

int main (int argc, char *argv[])
{
    cout << fixed << setw(2) << setprecision(2) ;

    int width = 1024, height = 768;
    int samples = argc==2 ? atoi(argv[1])/4 : 1;

    Random rand;

    Ray cam(Vector3d(50,52,295.6), Vector3d(0,-0.042612,-1).normalized());
    Vector3d cx(width*.5135/height,0,0);
    Vector3d cy=cx.cross(cam.direction).normalized()*.5135;
    vector<Color> c(width * height);

    sf::Image image;
    image.create(width, height, sf::Color::Black);

    // This is a pretty straightforward port from smallpt
    // TODO: Make this less horrible.
#pragma omp parallel for schedule(dynamic) private(rand)
    for (int y=0; y<height; y++) {
#pragma omp critical
        cout << "\rRendering (" << samples*4 << " spp) " << 100.*y/(height-1) << "%" << flush;
        for (int x=0; x<width; x++) {
            int i=(height-y-1)*width+x;
            for (int sy=0; sy<2; sy++) {
                for (int sx=0; sx<2; sx++){
                    Color r(0, 0, 0);
                    for (int s=0; s<samples; s++) {
                        double r1=2*rand.unit_rand();
                        double dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
                        double r2=2*rand.unit_rand();
                        double dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);

                        Vector3d d = cx*(((sx+.5 + dx)/2 + x)/width - .5) + cy*(((sy+.5 + dy)/2 + y)/height - .5) + cam.direction;
                        d = d.normalized();
                        Ray ray(cam.origin + d*140, d);
                        r += radiance(scene, ray, 0, rand)/samples;
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] += clamp_intensity(r)*.25;
                }
            }
            image.setPixel(x, (height-1)-y, toPPM(c[i]));
        }
    }

    if (!image.saveToFile("image.png")) {
        return -1;
    }
}
// Local Variables:
// irony-additional-clang-options: ("-std=c++17")
// End:
