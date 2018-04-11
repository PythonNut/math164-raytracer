#include "scene.h"

typedef LambertBSDF_IS Diffuse;
typedef GGXConductorBSDF_IS Glossy;

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
            move(make_unique<Glossy>(Glossy(Color(1, 1, 1)*.999,
                                            Color(0.155, 0.424, 1.383),
                                            Color(3.602, 2.472, 1.916),
                                            0.001)))),

     // glass
     Object(move(make_unique<Sphere>(Sphere(16.5, Vector3d(73, 16.5, 78)))),
            move(make_unique<Glossy>(Glossy(Color(1, 1, 1)*.999,
                                            Color(0.155, 0.424, 1.383),
                                            Color(3.602, 2.472, 1.916),
                                            0.03)))),

     // light
     Object(move(make_unique<Sphere>(Sphere(600, Vector3d(50, 681.6-.27, 81.6)))),
            move(make_unique<EmissionBSDF>(EmissionBSDF(Color(12, 12, 12)))))
};

const vector<Object> scene{make_move_iterator(begin(scene_temp)), make_move_iterator(end(scene_temp))};

const Ray camera(Vector3d(50,52,295.6), Vector3d(0,-0.042612,-1).normalized());
