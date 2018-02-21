#include <iostream>
#include <math.h>
#include <vector>
#include <tuple>
#include <random>
#include <eigen3/Eigen/Dense>

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

     double intersect(const Ray &ray) const {
          auto op = this->origin - ray.origin;
          double epsilon = 1e-4;
          double b = op.dot(ray.direction);
          double r2 = this->radius * this->radius;
          double descriminant = b*b - op.dot(op) + r2;

          if (descriminant < 0) {
               return 0;
          }

          descriminant = sqrt(descriminant);

          if (b - descriminant > epsilon) {
               return b - descriminant;
          }
          if (b + descriminant > epsilon) {
               return b + descriminant;
          }

          return 0;
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

double clamp(double x) {
     return fmin(fmax(0, x), 1);
}

int toPPM(double x) {
     return pow(clamp(x), 1/2.2) * 255 + 0.5;
}

const optional<pair<double, vector<Sphere>::const_iterator>> intersect(const Ray& ray,
                                                                       const vector<Sphere>& spheres) {
     double inf = std::numeric_limits<double>::infinity();
     double nearest_distance = inf;
     auto nearest_sphere = spheres.cend();

     for (auto it=spheres.begin(); it != spheres.end(); it++) {
          double distance = it->intersect(ray);
          if (distance < nearest_distance) {
               nearest_sphere = it;
               nearest_distance = distance;
          }
     }

     if (nearest_sphere == spheres.end()) {
          return {};
     }

     return make_pair(nearest_distance, nearest_sphere);
}

Color radiance(const Ray& ray, int depth, default_random_engine rand) {
     auto hit_check = intersect(ray, spheres);
     if (!hit_check.has_value() || depth > 10) {
          return Color(0, 0, 0);
     }

     double distance;
     vector<Sphere>::const_iterator sphere;
     tie(distance, sphere) = hit_check.value();

     Vector3d intersect_point = ray.origin + ray.direction * distance;
     Vector3d normal = (intersect_point - sphere->get_origin()).normalized();
     Vector3d oriented_normal = normal.dot(ray.direction) < 0 ? normal : -normal;

     Material mat = sphere->get_material();

     Color f = mat.get_color();

     // this might be slow...
     uniform_real_distribution<> uniform_rand(0, 1);

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
          return mat.get_emission() + f * radiance(r, depth, rand);
     }

     case Material::Reflection::SPECULAR: {
          Ray r(intersect_point, ray.direction - normal*2*normal.dot(ray.direction));
          return mat.get_emission() + f * radiance(r, depth, rand);
     }
     case Material::Reflection::REFRACTIVE: {
          Ray refl_ray(intersect_point, ray.direction - normal*2*normal.dot(ray.direction));
          bool into = normal.dot(oriented_normal) > 0;

          // compute the IOR ratio
          double nc=1, nt=1.5, nnt=into ? nc/nt : nt/nc;
          double ddn = ray.direction.dot(oriented_normal);
          double cos2t = 1- nnt * nnt * (1 - ddn * ddn);
          if (cos2t < 0) {
               return mat.get_emission() + f * radiance(refl_ray, depth, rand);
          }

          // fresnel math makes me sad
          Vector3d tdir = (ray.direction * nnt - normal*((into?1:-1)*(ddn * nnt * sqrt(cos2t)))).normalized();
          Ray trans_ray(intersect_point, tdir);
          double a=nt - nc, b=nt + nc, R0 = (a*a)/(b*b), c = 1-(into?-ddn:tdir.dot(normal));
          // TODO: Make this less awful
          double Re=R0+(1-R0)*pow(c, 5), Tr=1-Re, P=.25+.5*Re, RP=Re/P, TP=Tr/(1-P);
          if (depth > 2) {
               if (uniform_rand(rand) < P) {
                    return mat.get_emission() + radiance(refl_ray, depth, rand) * RP;
               } else {
                    return mat.get_emission() + radiance(trans_ray, depth, rand) * TP;
               }
          } else {
               return mat.get_emission() + radiance(trans_ray, depth, rand) * Tr;
          }
     }
     }
}

int main ()
{
     cout << "Hello World!" << endl;
}
// Local Variables:
// irony-additional-clang-options: ("-std=c++17")
// End:
