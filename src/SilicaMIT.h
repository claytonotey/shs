#ifndef SILICAMIT_H
#define SILICAMIT_H

#include "Dielectric.h"

#define SILICAMIT_TABLE_SIZE 4666

class SilicaMIT : public Dielectric {
 public:
  Complex eps(Real eV);
  Complex mu(Real eV);
  static Real lambdatab[SILICAMIT_TABLE_SIZE];
  static Real epsrtab[SILICAMIT_TABLE_SIZE];
  static Real epsitab[SILICAMIT_TABLE_SIZE];
};

#endif
