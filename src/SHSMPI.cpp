#include "SphereHalfspaceMPI.h"
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
   vector<Real> freqs;
   vector<Real> opts;
   
   pushFileLines("a.txt",as);
   pushFileLines("gaps.txt",gaps);
   pushFileLines("freqs.txt",freqs);
   pushFileLines("options.txt",opts);
   
   Silica diel1;
   Silica diel2;
   //SiliconCarbide diel1;
   //SiliconCarbide diel2;
   //SiliconDoped diel2(SiliconDopedP,opts[4]);

   Real kappa0 = opts[0];
   Real kappa1 = opts[1];
   Real kappa2 = opts[2];
   Real eta = opts[3];
   //const char *shapeFile = "shape.x3d";
   
   SHSMPI(as,gaps,freqs,kappa0,kappa1,kappa2,eta/*,shapeFile*/,&diel1,&diel2);
}
