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

     for (int s=0; s<samples; s++) {
          cout << "\rRendering (" << samples << "%) " << 100.0*(s+1)/samples << "%" << flush;
          #pragma omp parallel for schedule(dynamic) private(rand)
          for (int y=0; y<height; y++) {
               for (int x=0; x<width; x++) {
                    double dx, dy;
                    tie(dx, dy) = sampler.sample();
                    Vector3d d = cx*((dx/2 + x)/width - .5) + cy*((dy/2 + y)/height - .5) + camera.direction;
                    d.normalize();
                    Ray ray(camera.origin + d*140, d);
                    Color r = radiance(scene, ray, 0, rand);
                    m(height-y-1, x) += r;
               }
          }

          string math = "output";
          string var = to_string(s+1);
          ofstream myfile;
          myfile.open(math+var+".csv");
          for(int count = 0; count < width*height; count ++){
               myfile << m(count).x()/(s+1) << "\t"<< m(count).y()/(s+1) <<"\t" <<(m(count)).z()/(s+1)<<endl;
          }
          myfile.close();
     }

     ofstream out_file;
     out_file.open("imageDsP.ppm");
     out_file << "P3" << endl << width << " " << height << endl << 255 << endl;
     for (int j=0; j<height; j++){
          for (int i=0; i<width; i++) {
               out_file << toPPM((m(j,i)/samples).x()) << " " << toPPM((m(j,i)/samples).y()) << " " << toPPM((m(j,i)/samples).z()) << " ";
          }
     }
     out_file.close();
     return 0;
}
