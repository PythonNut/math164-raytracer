#include <iostream>
#include <iomanip>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include "util/linear_algebra.h"
#include "util/util.h"
#include "raytrace/scene.h"
#include "raytrace/radiance.h"
#include "sample/marsaglia_polar_sampler.h"

sf::Color toPPM(Color c) {
     return sf::Color(toPPM(c.x()), toPPM(c.y()), toPPM(c.z()));
}

int main (int argc, char *argv[])
{
     cout << fixed << setw(2) << setprecision(2) ;

     int width = 1024, height = 768;
     int samples = argc==2 ? atoi(argv[1]) : 1;

     Random rand;
     MarsagliaPolarSampler sampler(rand);

     Vector3d cx(width*.5135/height,0,0);
     Vector3d cy=cx.cross(camera.direction).normalized()*.5135;
     vector<Color> c(width * height);

     sf::Image image;
     image.create(width, height, sf::Color::Black);

     // This is a pretty straightforward port from smallpt
     // TODO: Make this less horrible.
#pragma omp parallel for schedule(dynamic) private(rand)
     for (int y=0; y<height; y++) {
#pragma omp critical
          cout << "\rRendering (" << samples << " spp) " << 100.*y/(height-1) << "%" << flush;
          for (int x=0; x<width; x++) {
               Color r(0, 0, 0);
               for (int s=0; s<samples; s++) {
                    double dx, dy;
                    tie(dx, dy) = sampler.sample();
                    Vector3d d = cx*((dx/2 + x)/width - .5) + cy*((dy/2 + y)/height - .5) + camera.direction;
                    d = d.normalized();
                    Ray ray(camera.origin + d*140, d);
                    r += radiance(scene, ray, 0, rand)/samples;
               } // Camera rays are pushed ^^^^^ forward to start in interior
               image.setPixel(x, (height-1)-y, toPPM(r ));
          }
     }

     if (!image.saveToFile("image_progressive.png")) {
          return -1;
     }
}
