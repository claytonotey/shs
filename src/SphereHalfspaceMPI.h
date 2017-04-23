#ifndef SPHEREHALFSPACEMPI_H
#define SPHEREHALFSPACEMPI_H

#include "Util.h"
#include "Dielectric.h"
#include <vector>
using namespace std;

void SHSMPI(vector<Real> &as, vector<Real> &gaps, vector<Real> &freqs, Real kappa0, Real kappa1, Real kappa2, Real eta, /*const char *shapeFile,*/ Dielectric *diel1, Dielectric *diel2);

#endif
