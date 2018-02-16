#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

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
     Vector3d color, emission;
     Reflection reflection;

public:
     Material(Vector3d col, Vector3d emit, Reflection ref)
          : color(col)
          , emission(emit)
          , reflection(ref) {}
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
};

const vector<Sphere> spheres = {
     // left
     Sphere(1e5,
            Vector3d(1e5+1, 40.8, 81.6),
            Material(Vector3d(0, 0, 0),
                     Vector3d(0.75, 0.25, 0.25),
                     Material::Reflection::DIFFUSE)),

     // right
     Sphere(1e5,
            Vector3d(-1e5+99, 40.8, 81.6),
            Material(Vector3d(0, 0, 0),
                     Vector3d(0.25, 0.25, 0.75),
                     Material::Reflection::DIFFUSE)),

     // back
     Sphere(1e5,
            Vector3d(50,40.8, 1e5),
            Material(Vector3d(0, 0, 0),
                     Vector3d(0.75, 0.75, 0.75),
                     Material::Reflection::DIFFUSE)),

     // front
     Sphere(1e5,
            Vector3d(50, 40.8, -1e5+170),
            Material(Vector3d(0, 0, 0),
                     Vector3d(0, 0, 0),
                     Material::Reflection::DIFFUSE)),

     // bottom
     Sphere(1e5,
            Vector3d(50, 1e5, 81.6),
            Material(Vector3d(0, 0, 0),
                     Vector3d(0.75, 0.75, 0.75),
                     Material::Reflection::DIFFUSE)),

     // top
     Sphere(1e5,
            Vector3d(50, -1e5+81.6, 81.6),
            Material(Vector3d(0, 0, 0),
                     Vector3d(0.75, 0.75, 0.75),
                     Material::Reflection::DIFFUSE)),

     // mirror
     Sphere(16.5,
            Vector3d(27, 16.5, 47),
            Material(Vector3d(0, 0, 0),
                     Vector3d(1,1,1)*.999,
                     Material::Reflection::SPECULAR)),

     // glass
     Sphere(16.5,
            Vector3d(73, 16.5, 78),
            Material(Vector3d(0, 0, 0),
                     Vector3d(1,1,1)*.999,
                     Material::Reflection::REFRACTIVE)),

     // light
     Sphere(600,
            Vector3d(50, 681.6-.27, 81.6),
            Material(Vector3d(12, 12, 12),
                     Vector3d(0, 0, 0),
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



int main ()
{
     cout << "Hello World!" << endl;
}
