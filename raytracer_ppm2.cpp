#include <iostream>
#include <iomanip>
#include <fstream>
#include "util/linear_algebra.h"
#include "util/util.h"
#include "raytrace/scene.h"
#include "raytrace/radiance.h"
#include "sample/marsaglia_polar_sampler.h"

#define USE_EIGEN
using namespace std;
using namespace Eigen;

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
     typedef Matrix<Color, Dynamic, Dynamic> Matrix;
     Matrix m(height, width);


     // This is a pretty straightforward port from smallpt
     // TODO: Make this less horrible.
for (int s=0; s<samples; s++) {
#pragma omp parallel for schedule(dynamic) private(rand)
     for (int y=0; y<height; y++) {
#pragma omp critical
          cout << "\rRendering (" << samples << " spp) " << 100.*y/(height-1) << "%" << flush;
          for (int x=0; x<width; x++) {
            //int i=(height-y-1)*width+x;
               Color r(0, 0, 0);
               //for (int s=0; s<samples; s++) {
                    double dx, dy;
                    tie(dx, dy) = sampler.sample();
                    Vector3d d = cx*((dx/2 + x)/width - .5) + cy*((dy/2 + y)/height - .5) + camera.direction;
                    d = d.normalized();
                    Ray ray(camera.origin + d*140, d);
                    r += radiance(scene, ray, 0, rand)/s;
               //} // Camera rays are pushed ^^^^^ forward to start in interior
               //c[i] += clamp_intensity(r)*.25;
               m(height-y-1, x) += clamp_intensity(r);
          }
     }
     ofstream out_file;
     out_file.open("imagePPM1.ppm");
     out_file << "P3" << endl << width << " " << height << endl << 255 << endl;
     for (int j=0; j<height; j++){
        for (int i=0; i<width; i++) {
          out_file << toPPM((m(j,i)).x()) << " " << toPPM((m(j,i)).y()) << " " << toPPM((m(j,i)).z()) << " ";
        }
     }
     out_file.close();

}
/*
ofstream out_file;
out_file.open("imageDsP.ppm");
out_file << "P3" << endl << width << " " << height << endl << 255 << endl;
for (int j=0; j<height; j++){
   for (int i=0; i<width; i++) {
     out_file << toPPM((m(j,i)).x()) << " " << toPPM((m(j,i)).y()) << " " << toPPM((m(j,i)).z()) << " ";
   }
}
out_file.close();
*/


}
