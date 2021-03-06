#include "SphereHalfspace.h"
#include "FileUtils.h"
#include <vector>
using namespace std;

#include "Dielectrics.h"

void usage() 
{
  fprintf(stderr,"Usage: shs radius(m) gap(m) fmin(eV) fmax(eV) T1(K) T2(K) dielectric1 dielectric2 [kappa0 kappa1 kappa2 eta]\n  Dielectrics are specified as one of {SiC, SiO2, SiN, Au, Al, SiliconDoped(N=?carrier/m^3)}.");
}


int main(int argc, char **argv) {
   Real a;
   Real gap;
   Real kappa0 = 8;
   Real kappa1 = 2.5;
   Real kappa2 = 0.8;
   Real eta = 3e-6;
   Real freq0;
   Real freq1;
   Real T1;
   Real T2;

   if(argc < 7) {
     usage();
     exit(-1);
   }
   a = atof(argv[1]);
   gap = atof(argv[2]);
   freq0 = atof(argv[3]);
   freq1 = atof(argv[4]);
   T1 = atof(argv[5]);
   T2 = atof(argv[6]);
   
   Dielectric *diel1 = Dielectrics::parseDielectric(string(argv[7]));
   Dielectric *diel2 = Dielectrics::parseDielectric(string(argv[8]));

   if(argc > 9) {
     if(argc < 13) {
       usage();
       exit(-1);
     } else {
       kappa0 = atof(argv[9]);
       kappa1 = atof(argv[10]);
       kappa2 = atof(argv[11]);
       eta = atof(argv[12]);
     }
   }

   SHS(a,gap,freq0,freq1,kappa0,kappa1,kappa2,eta,diel1,diel2,T1,T2);
}
