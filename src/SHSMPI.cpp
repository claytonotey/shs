#include "SphereHalfspaceMPI.h"
#include "FileUtils.h"
#include <vector>
using namespace std;

#include "Dielectrics.h"

int main(int argc, char **argv) {
   vector<Real> as;
   vector<Real> gaps;
   vector<Real> freqs;
   vector<string> opts;
   
   pushFileLines("a.txt",as);
   pushFileLines("gaps.txt",gaps);
   pushFileLines("freqs.txt",freqs);
   pushFileLines("options.txt",opts);
   
   Real kappa0 = stof(opts[0]);
   Real kappa1 = stof(opts[1]);
   Real kappa2 = stof(opts[2]);
   Real eta = stof(opts[3]);
   Dielectric *diel1 = Dielectrics::parseDielectric(opts[4]);
   Dielectric *diel2 = Dielectrics::parseDielectric(opts[5]);
   //const char *shapeFile = "shape.x3d";
   
   SHSMPI(as,gaps,freqs,kappa0,kappa1,kappa2,eta/*,shapeFile*/,diel1,diel2);
}
