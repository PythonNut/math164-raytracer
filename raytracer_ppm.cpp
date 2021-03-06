#include <iostream>
#include <iomanip>
#include <fstream>
#include "util/linear_algebra.h"
#include "util/util.h"
#include "raytrace/scene.h"
#include "raytrace/radiance.h"

#define USE_EIGEN
using namespace std;
using namespace Eigen;

int main (int argc, char *argv[])
{
     cout << fixed << setw(2) << setprecision(2) ;

     int width = 1024, height = 768;
     int samples = argc==2 ? atoi(argv[1])/4 : 1;

     Random rand;

     Vector3d cx(width*.5135/height,0,0);
     Vector3d cy=cx.cross(camera.direction).normalized()*.5135;
     vector<Color> c(width * height);
     typedef Matrix<Color, Dynamic, Dynamic> Matrix;
     Matrix m(height, width);

     // This is a pretty straightforward port from smallpt
     // TODO: Make this less horrible.
#pragma omp parallel for schedule(dynamic) private(rand)
     for (int y=0; y<height; y++) {
//#pragma omp critical
          cout << "\rRendering (" << samples*4 << " spp) " << 100.*y/(height-1) << "%" << flush;
          for (int x=0; x<width; x++) {
               //int i=(height-y-1)*width+x;
               for (int sy=0; sy<2; sy++) {
                    for (int sx=0; sx<2; sx++){
                         Color r(0, 0, 0);
                         for (int s=0; s<samples; s++) {
                              double r1=2*rand.unit_rand();
                              double dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
                              double r2=2*rand.unit_rand();
                              double dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);

                              Vector3d d = cx*(((sx+.5 + dx)/2 + x)/width - .5) + cy*(((sy+.5 + dy)/2 + y)/height - .5) + camera.direction;
                              d = d.normalized();
                              Ray ray(camera.origin + d*140, d);
                              r += radiance(scene, ray, 0, rand)/samples;
                         } // Camera rays are pushed ^^^^^ forward to start in interior
                         m(height-y-1, x) += clamp_intensity(r)*.25;

                    }
               }
          }
     }


     ofstream myfile;
     string str = "IS10K";
     myfile.open(str+".csv");
     for(int count = 0; count < width*height; count ++){
          myfile << m(count).x() << "\t"<< m(count).y() <<"\t" <<(m(count)).z()<<endl;
     }
     myfile.close();

/*
     ofstream out_file;
     out_file.open("image.ppm");
     out_file << "P3" << endl << width << " " << height << endl << 255 << endl;
     for (int j=0; j<height; j++){
          for (int i=0; i<width; i++) {
               out_file << toPPM((m(j,i)).x()) << " " << toPPM((m(j,i)).y()) << " " << toPPM((m(j,i)).z()) << " ";
          }
     }
     out_file.close();
     */
}
