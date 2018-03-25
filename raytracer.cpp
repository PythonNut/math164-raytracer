#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <tuple>
#include <random>
#include <functional>
#include <eigen3/Eigen/Dense>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

using namespace Eigen;
using namespace std;

typedef Array3d Color;

struct Ray {
    Vector3d origin, direction;
    Ray(Vector3d orig, Vector3d dir)
        : origin(orig)
        , direction(dir) { }
};

class Material {
public:
    enum class Reflection {
        DIFFUSE,
        SPECULAR,
        REFRACTIVE
    };

private:
    Color color, emission;
    Reflection reflection;

public:
    Material(Color emit, Color col, Reflection ref)
        : color(col)
        , emission(emit)
        , reflection(ref) {}

    Color get_color() const {
        return this->color;
    }

    Color get_emission() const {
        return this->emission;
    }


    Reflection get_reflection() const {
        return this->reflection;
    }
};

class Sphere {
private:
    double radius;
    Vector3d origin;
    Material material;
public:
    Sphere(double r, Vector3d orig, Material mat)
        : radius(r)
        , origin(orig)
        , material(mat) { }

    optional<double> intersect(const Ray &ray) const {
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

    Vector3d get_normal(Vector3d point) const {
        return (point - this->get_origin()).normalized();
    }

    double get_radius() const {
        return this->radius;
    }

    Vector3d get_origin() const {
        return this->origin;
    }

    Material get_material() const {
        return this->material;
    }
};

const vector<Sphere> spheres = {
    // left
    Sphere(1e5,
           Vector3d(1e5+1, 40.8, 81.6),
           Material(Color(0, 0, 0),
                    Color(0.75, 0.25, 0.25),
                    Material::Reflection::DIFFUSE)),

    // right
    Sphere(1e5,
           Vector3d(-1e5+99, 40.8, 81.6),
           Material(Color(0, 0, 0),
                    Color(0.25, 0.25, 0.75),
                    Material::Reflection::DIFFUSE)),

    // back
    Sphere(1e5,
           Vector3d(50,40.8, 1e5),
           Material(Color(0, 0, 0),
                    Color(0.75, 0.75, 0.75),
                    Material::Reflection::DIFFUSE)),

    // front
    Sphere(1e5,
           Vector3d(50, 40.8, -1e5+170),
           Material(Color(0, 0, 0),
                    Color(0, 0, 0),
                    Material::Reflection::DIFFUSE)),

    // bottom
    Sphere(1e5,
           Vector3d(50, 1e5, 81.6),
           Material(Color(0, 0, 0),
                    Color(0.75, 0.75, 0.75),
                    Material::Reflection::DIFFUSE)),

    // top
    Sphere(1e5,
           Vector3d(50, -1e5+81.6, 81.6),
           Material(Color(0, 0, 0),
                    Color(0.75, 0.75, 0.75),
                    Material::Reflection::DIFFUSE)),

    // mirror
    Sphere(16.5,
           Vector3d(27, 16.5, 47),
           Material(Color(0, 0, 0),
                    Color(1,1,1)*.999,
                    Material::Reflection::SPECULAR)),

    // glass
    Sphere(16.5,
           Vector3d(73, 16.5, 78),
           Material(Color(0, 0, 0),
                    Color(1,1,1)*.999,
                    Material::Reflection::REFRACTIVE)),

    // light
    Sphere(600,
           Vector3d(50, 681.6-.27, 81.6),
           Material(Color(12, 12, 12),
                    Color(0, 0, 0),
                    Material::Reflection::DIFFUSE))
};

double clamp_intensity(const double x) {
    return fmin(fmax(0, x), 1);
}

int toPPM(const double x) {
    return pow(clamp_intensity(x), 1/2.2) * 255 + 0.5;
}

const optional<pair<double, vector<Sphere>::const_iterator>> intersect(const Ray& ray,
                                                                       const vector<Sphere>& scene) {
    double nearest_distance = std::numeric_limits<double>::infinity();
    auto nearest_object = scene.end();

    for (auto it=scene.begin(); it != scene.end(); ++it) {
        auto hit_check = it->intersect(ray);
        if (hit_check.has_value()) {
            double distance = hit_check.value();
            if (distance < nearest_distance) {
                nearest_object = it;
                nearest_distance = distance;
            }
        }
    }

    if (nearest_object == spheres.end()) {
        return {};
    }

    return make_pair(nearest_distance, nearest_object);
}

const Color radiance(const Ray& ray,
                     int depth,
                     default_random_engine rand,
                     uniform_real_distribution<> uniform_rand) {
    auto hit_check = intersect(ray, spheres);
    if (!hit_check.has_value() || depth > 10) {
        return Color(0, 0, 0);
    }

    double distance;
    vector<Sphere>::const_iterator sphere;
    tie(distance, sphere) = hit_check.value();

    Vector3d intersect_point = ray.origin + ray.direction * distance;
    Vector3d normal = sphere->get_normal(intersect_point);
    Vector3d oriented_normal = normal;
    if (normal.dot(ray.direction) > 0) {
        oriented_normal *= -1;
    }

    auto mat = sphere->get_material();

    Color f = mat.get_color();

    double p = f.maxCoeff();
    if (++depth>5 || !p) {
        if (uniform_rand(rand) < p) {
            f /= p;
        } else {
            return mat.get_emission();
        }
    }

    switch (mat.get_reflection()) {
    case Material::Reflection::DIFFUSE: {
        double r1 = uniform_rand(rand) * 2 * M_PI;
        double r2 = uniform_rand(rand), r2s = sqrt(r2);
        Vector3d w = oriented_normal;
        Vector3d u = fabs(w.x()) > 0.1 ? Vector3d(0, 1, 0) : Vector3d(1, 0, 0);
        u = u.cross(w).normalized();
        Vector3d v = w.cross(u);
        Vector3d d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).normalized();

        Ray r(intersect_point, d);
        return mat.get_emission() + f * radiance(r, depth, rand, uniform_rand);
    }

    case Material::Reflection::SPECULAR: {
        Ray r(intersect_point, ray.direction - normal*2*normal.dot(ray.direction));
        return mat.get_emission() + f * radiance(r, depth, rand, uniform_rand);
    }

    case Material::Reflection::REFRACTIVE: {
        Ray refl_ray(intersect_point, ray.direction - normal*2*normal.dot(ray.direction));
        bool into = normal.dot(oriented_normal) > 0;

        // compute the IOR ratio
        double nc=1, nt=1.5, nnt=into ? nc/nt : nt/nc;
        double ddn = ray.direction.dot(oriented_normal);
        double cos2t = 1- nnt * nnt * (1 - ddn * ddn);
        if (cos2t < 0) {
            return mat.get_emission() + f * radiance(refl_ray, depth, rand, uniform_rand);
        }

        // fresnel math makes me sad
        Vector3d tdir = (ray.direction*nnt - normal*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).normalized();
        Ray trans_ray(intersect_point, tdir);
        double a=nt-nc, b=nt+nc, R0=(a*a)/(b*b), c=1-(into?-ddn:tdir.dot(normal));
        // TODO: Make this less awful

        double Re=R0+(1-R0)*pow(c, 5);
        Color result;
        if (depth > 2) {
            double Tr=1-Re, P=.25+.5*Re, RP=Re/P, TP=Tr/(1-P);
            if (uniform_rand(rand) < P) {
                result = radiance(refl_ray, depth, rand, uniform_rand) * RP;
            } else {
                result = radiance(trans_ray, depth, rand, uniform_rand) * TP;
            }
        } else {
            if (uniform_rand(rand) < Re) {
                result = radiance(refl_ray, depth, rand, uniform_rand);
            } else {
                result = radiance(trans_ray, depth, rand, uniform_rand);
            }
        }
        return mat.get_emission() + f * result;
    }
    default: return mat.get_emission();
    }
}

int main (int argc, char *argv[])
{
    cout << fixed << setw(2) << setprecision(2) ;

    int width = 1024, height = 768;
    int samples = argc==2 ? atoi(argv[1])/4 : 1;

    random_device rand_dev;
    default_random_engine rand_engine(rand_dev());
    uniform_real_distribution<> uniform_rand(0, 1);

    Ray cam(Vector3d(50,52,295.6), Vector3d(0,-0.042612,-1).normalized());
    Vector3d cx(width*.5135/height,0,0);
    Vector3d cy=cx.cross(cam.direction).normalized()*.5135;
    vector<Color> c(width * height);

    sf::Image image;
    image.create(width, height, sf::Color::Black);

    // This is a pretty straightforward port from smallpt
    // TODO: Make this less horrible.
#pragma omp parallel for schedule(dynamic)
    for (int y=0; y<height; y++) {
#pragma omp critical
        cout << "\rRendering (" << samples*4 << " spp) " << 100.*y/(height-1) << "%" << flush;
        for (int x=0; x<width; x++) {
            int i=(height-y-1)*width+x;
            for (int sy=0; sy<2; sy++) {
                for (int sx=0; sx<2; sx++){
                    Color r(0, 0, 0);
                    for (int s=0; s<samples; s++) {
                        double r1=2*uniform_rand(rand_engine);
                        double dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
                        double r2=2*uniform_rand(rand_engine);
                        double dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);

                        Vector3d d = cx*(((sx+.5 + dx)/2 + x)/width - .5) + cy*(((sy+.5 + dy)/2 + y)/height - .5) + cam.direction;
                        Ray ray(cam.origin + d*140, d.normalized());
                        r += radiance(ray, 0, rand_engine, uniform_rand)/samples;
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] += r.unaryExpr(ptr_fun(clamp_intensity))*.25;
                }
            }
            sf::Color color(toPPM(c[i].x()), toPPM(c[i].y()), toPPM(c[i].z()));
            image.setPixel(x, (height-1)-y, color);
        }
    }

    if (!image.saveToFile("image.png")) {
        return -1;
    }
}
// Local Variables:
// irony-additional-clang-options: ("-std=c++17")
// End:
