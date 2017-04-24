#include "SphereHalfspace.h"
#include "FileUtils.h"
#include <vector>
using namespace std;

#include "Dielectrics.h"

void usage() 
{
  fprintf(stderr,"Usage: shs radius(m) gap(m) fmin(eV) fmax(eV) T1(K) T2(K) dielectric1 dielectric2 [kappa0 kappa1 kappa2 eta]\n  Dielectrics are specified asone of {SiC, SiN, Au, Al, SiliconDoped(N=?carrier/m^3)}.");
}


Params parseParams(const string &s) 
{
  Params params;
  int pos0 = 0;
  int pos1;
  int pos2;
  do {
    pos1 = s.find_first_of("=", pos0);
    pos2 = s.find_first_of(",", pos1);
    if(pos1 == string::npos) {
      usage();
      exit(-1);
    } else {
      string key = s.substr(pos0,pos1-pos0);
      string val = s.substr(pos1+1,pos2-pos1-1);
      params[key] = Params::parseParam(val);
    }
    if(pos2 == string::npos) {
      break;
    }
    pos0 = pos2+1;
  } while(true);
  return params;
}

Dielectric *parseDielectric(const string &s)
{
  int pos = s.find_first_of("(");
   if(pos == string::npos) {
     return Dielectrics::getDielectric(s);
   } else {
     string name = s.substr(0,pos);
     int pos2 = s.find_first_of(")");
     if(pos2 == string::npos) {
       usage();
       exit(-1);
     } else {
       string param = s.substr(pos+1, pos2-pos-1);
       Params params = parseParams(param);
       return Dielectrics::getDielectric(name, params);
     }
   }
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
   
   Dielectric *diel1 = parseDielectric(string(argv[7]));
   Dielectric *diel2 = parseDielectric(string(argv[8]));

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
