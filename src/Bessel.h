#ifndef BESSEL_H
#define BESSEL_H

#include "MathUtils.h"
#include "Array.h"

template<class R, class T>
  T SphericalBesselSeriesJ(const T &z, int N) 
{
  T b = 1.;
  T bacc = 1.;
  long long int n = 4 * N + 6;
  T z2 = -z * z;
  long long int m = 1;
  while(norm(b/bacc) > norm(numeric_limits<R>::epsilon())) {
    b *= z2 / (R)(n * m);
    bacc += b;
    n += 4;
    m++;
  }
  return bacc;
}

template<class R, class T>
  T SphericalBesselSeriesH(const T &z, int N) 
{
  T b = 1.;
  T bacc = 1.;
  long long int n = 4 * N - 2;
  T z2 = z * z;
  long long int m = 1;
  while(norm(b/bacc) > norm(numeric_limits<R>::epsilon())) {
    b *= z2 / (R)(n * m);
    bacc += b;
    n -= 4;
    m++;
  }
  return bacc;
}

template<class R, class T>
T SphericalBesselDJoverJ(const T &z,
                           int n) 
{
  return ((R)(2 * n + 1) * SphericalBesselSeriesJ<R,T>(z,n-1)) / (z * SphericalBesselSeriesJ<R,T>(z,n)) - (R)n / z;
}

template<class R, class T>
T SphericalBesselDHoverH(const T &z,
                         int n)
{
  return (z * SphericalBesselSeriesH<R,T>(z,n-1)) / ((R)(2 * n - 1) * SphericalBesselSeriesH<R,T>(z,n)) - (R)n / z;
}

template<class R, class T>
void SphericalBesselSignJ(const T &z,
                          int N,
                          int N0,
                          const T &jN0,
                          T *sj) 
{
  T s = jN0 / sqrt(norm(jN0));
  T den = SphericalBesselSeriesJ<R,T>(z,N0);
  T num;

  for(int n=N0+1; n<=N; n++) {
    num = SphericalBesselSeriesJ<R,T>(z,n);
    s *= z * num/den;
    s /= sqrt(norm(s));
    sj[n] = s;
    den = num;
  }
}

template<class R, class T>
void SphericalBesselHAbsJ(const T &z,
                          int N,
                          int N0,
                          const T &jN0,
                          const T &hN0,
                          T *haj) 
{
  T s = sqrt(norm(jN0)) * hN0;
  T denj = SphericalBesselSeriesJ<R,T>(z,N0);
  T numj;
  T denh = SphericalBesselSeriesH<R,T>(z,N0);;
  T numh;
  int nj = (2 * N0 + 1);
  int nh = (2 * N0 - 1);
  for(int n=N0+1; n<=N; n++) {
    nj += 2;
    nh += 2;
    numj = SphericalBesselSeriesJ<R,T>(z,n);
    numh = SphericalBesselSeriesH<R,T>(z,n);
    s *= (R)nh * numh / ((R)nj * denh * z) * sqrt(norm(numj / denj * z));
    haj[n] = s;
    denj = numj;
    denh = numh;
  }
}

void BesselJ(const Complex &z,
             double nu,
             int n,
             Complex *j);

void BesselH(const Complex &z,
             double nu,
             int n,
             Complex *h);

void SphericalBesselJ(const Complex &z,
                      Real nu,
                      int n,
                      Complex *j,
                      Complex *dj);

void SphericalBesselH(const Complex &z,
                      Real nu,
                      int n,
                      Complex *h,
                      Complex *dh);


inline bool isBesselAsymptotic(const Complex &z, int l)
{
  if(l < 4)
    return false;
  Complex j = 1.;
  for(int n=0; n<l; n++) {
    j *= z / (Real)(2 * n + 3);
    if(norm(j) < 1e-256) return true;
  }
  return false;
}

template<class R, class T>
class SphericalBessel {
 public:

  SphericalBessel(const T &k_, int p_, int n_, int n0_, bool bInverse_ = false) :
    k(k_),
    p(p_),
    n(n_),
    n0(n0_),
    z(n0_+1),
    x(n0_+1),
    r(n_+1),
    d(n_+1),
    bInverse(bInverse_)  
    {
    if(p==1) {  
      SphericalBesselJ(k, (Real)1., n0, z+1, x+1);
      for(int l=1; l<=n; l++) {
        if(l > n0) {
          r[l] = k / (R)(2 * l + 1) * SphericalBesselSeriesJ<R,T>(k,l) / SphericalBesselSeriesJ<R,T>(k,l-1);
          d[l] = SphericalBesselDJoverJ<R,T>(k,l);
        } else {
          r[l] = z[l] / z[l-1];
          d[l] = x[l] / z[l];
        }
      }
    } else if(p==3) {
      SphericalBesselH(k, (Real)1., n0, z+1, x+1);
      for(int l=1; l<=n; l++) {
        if(l > n0) {
          r[l] = (R)(2 * l - 1) / k * SphericalBesselSeriesH<R,T>(k,l) / SphericalBesselSeriesH<R,T>(k,l-1);
          d[l] = SphericalBesselDHoverH<R,T>(k,l);
        } else {
          r[l] = z[l] / z[l-1];
          d[l] = x[l] / z[l];
        }
      }
    } else {
      abort();
    }
    
    if(bInverse) {
      for(int l=1; l<=n0; l++) {
        z[l] = (Real)1. / z[l];
        x[l] = (Real)1. / x[l];
      }
      for(int l=1; l<=n; l++) {
        r[l] = (Real)1. / r[l];
        d[l] = (Real)1. / d[l];
      }
    }
  }

  T k;
  int p;
  int n;
  int n0;
  Array<T> z;
  Array<T> x;
  Array<T> r;
  Array<T> d;
  bool bInverse;
};

#endif
