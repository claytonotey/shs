#include "SphereHalfspace.h"
#include "FileUtils.h"
#include <vector>
using namespace std;

#include "Silica.h"
#include "SilicaMIT.h"
#include "SiliconCarbide.h"
#include "SiliconDoped.h"

int main(int argc, char **argv) {
   vector<Real> as;
   vector<Real> gaps;

   Real kappa0 = 8;
   Real kappa1 = 2.5;
   Real kappa2 = 0.8;
   Real eta = 3e-6;
   Real freq0;
   Real freq1;

   if(argc < 5) {
     fprintf(stderr,"Usage: shs a fmin(eV) fmax(eV) [gap kappa0 kappa1 kappa2 eta]\n");
     exit(-1);
   }
   as.push_back(atof(argv[1]));
   gaps.push_back(atof(argv[2]));

   freq0 = atof(argv[3]);
   freq1 = atof(argv[4]);

   if(argc > 5) {
     if(argc < 9) {
       fprintf(stderr,"Usage: shs a fmin(eV) fmax(eV) [gap kappa0 kappa1 kappa2 eta]\n");
       exit(-1);
     } else {
       kappa0 = atof(argv[5]);
       kappa1 = atof(argv[6]);
       kappa2 = atof(argv[7]);
       eta = atof(argv[8]);
     }
   }

   Silica diel1;
   Silica diel2;
   //SiliconCarbide diel1;
   //SiliconCarbide diel2;
   //SiliconDoped diel2(SiliconDopedP,opts[4]);
   //const char *shapeFile = "shape.x3d";

   Real T1 = 300.;
   Real T2 = 290.;

   SHS(as,gaps,freq0,freq1,kappa0,kappa1,kappa2,eta,&diel1,&diel2,T1,T2);
}
