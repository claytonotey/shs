#ifndef SPHEREHALFSPACE_H
#define SPHEREHALFSPACE_H

#include "MathUtils.h"
#include "Dielectric.h"
#include <vector>
using namespace std;

void SHS(Real a, Real gap, Real freq0, Real freq1, Real kappa0, Real kappa1, Real kappa2, Real eta, Dielectric *diel1, Dielectric *diel2, double T1, double T2);

#endif
