#ifndef SILICONDOPED_H
#define SILICONDOPED_H

#include "Dielectric.h"

enum SiliconDopedType {
  SiliconDopedN,
  SiliconDopedP,
};

class SiliconDoped : public Dielectric {
 public:
  SiliconDoped(SiliconDopedType type, Real N);
  Complex eps(Real eV);
  Complex mu(Real eV);
  SiliconDopedType type;
  Real N;
};

#endif
