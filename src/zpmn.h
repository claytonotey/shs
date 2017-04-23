#ifndef ZPMN_H
#define ZPMN_H

#include "MathUtils.h"
#include <map>
#include <vector>
#include <utility>
#include <iostream>
using namespace std;

inline double real(double x) { return x; }
inline double imag(double x) { return 0.; }

class ZPMN {
 public:

  static Real sn[16384];
  static Real sd[16384];

  Real expArg;
  bool bReal;
  map<Real, vector<Complex> > memo;
  map<Real, int > nminMap;
  int rN0;
  Complex *c,*d,*rc,*rd;

  ZPMN(Real expArg, Complex *c, Complex *d, Complex *rc, Complex *rd, int rN0, bool bReal);

  inline void zpmn(int N0, int N, int M, 
                   Complex &zc, Real zs,
                   Complex &p1, Complex &p2,
                   Complex &p, Complex &dp);
  inline void eval(vector<Complex> &v, int N, int M, Real z, Real zs, Complex &p, Complex &dp);
  inline void eval(int N, int M, Real z, Real zs, bool bDBessel, Complex &p, Complex &dp);

};

#define SQRT2I 0.70710678118654752440084436210485

inline
void ZPMN :: zpmn(int N0, int N, int M, 
                  Complex &zc, Real zs,
                  Complex &p1, Complex &p2,
                  Complex &p, Complex &dp)
{
  Real Nr = (Real)N;
  Complex s = 1.;
  Complex e = 1.;

  Real zs2 = zs * zs;

  if(N > rN0 && N0 > rN0 && rc) {
    s = rc[N0];
  }
  if(N != N0) {
    if(N > 0 && N0 == 0) {
      e = exp(Complex(0.,1.)*expArg*zc/Nr);
    }
    if(N0 == 0) {
      p1 = 1.;
      p2 = 0.;
    }
    for(int n=N0; n<N; n++) {
      p2 *= s;
      if(N > rN0 && n >= rN0 && rc) {
        s = rc[n+1];
      }
      if(M > 0 && n < M) {
        p1 *= -sn[2*n+1]*sd[2*n+2]*zs*e*s;
        p2 = 0.;
      } else if(n == M) {
        p2 = p1 * e;
        p1 = sn[2*M+1]*zc*p2*s;
      } else {
        Complex tmp = p1;
        p1 = ((Real)(2*n+1)*zc*p1-sn[n+M]*sn[n-M]*p2) * sd[n+1-M] * sd[n+1+M] * e * s;
        p2 = tmp * e;
      }
    }
  }

  if(M == N) {
    p = p1;
    dp = (zs == 0.0?p1:-zc/zs2*Nr*p1);
  } else {
    p = p1;
    dp = (sn[N+M]*sn[N-M]*p2*s - zc*Nr*p1) / zs2;
  }    

  s = SQRT2I * ((M%2==0)?sn[2*N+1]:-sn[2*N+1]);
  p *= s;
  dp *= s;
}

#define MEMOEPS 1e-15

inline
void ZPMN :: eval(vector<Complex> &v, int N, int M, Real z, Real zs,
                  Complex &p, Complex &dp) {    

  Complex p1;
  Complex p2;
  Complex zc;
  if(bReal)
    zc = Complex(z,0.);
  else
    zc = Complex(0.,z);
  
  int nPrepend = 0;
  int nmin;
  int nminPrev;
  if(v.empty()) {
    nmin = N - 1;
    nminMap[z] = nmin;
    nminPrev = -1;
    nPrepend = 2;
  } else {
    nminPrev = nminMap[z];
    if(N <= nminPrev) {
      nmin = N - 1;
      nminMap[z] = nmin;
      nPrepend = nminPrev - nmin;
    } else {
      nmin = nminPrev;
    }
  }
  
  if(nPrepend) {
    v.insert(v.begin(),nPrepend,p1);
    zpmn(0,N,M,zc,zs,p1,p2,p,dp);
    v[0] = p2;
    v[1] = p1;
    for(int n=N; n<nminPrev-1; n++) {
      Complex p0, dp0;
      zpmn(n,n+1,M,zc,zs,p1,p2,p0,dp0);
      v[n+1-nmin] = p1;
    }
  } else {
    int size = v.size();
    v.resize(max(size,N-nmin+1));
    int nmax = size + nmin - 1;
    int n1 = min(nmax,N);
    int n2 = max(0,n1-1);
    int n0 = min(N,n1+1);
    
    p1 = v[n1-nmin];
    p2 = v[n2-nmin];
    
    do {
      zpmn(n1,n0,M,zc,zs,p1,p2,p,dp);
      v[n0-nmin] = p1;
      n0++;
      n1++;
    } while(n0 <= N);
  }
}

void ZPMN :: eval(int N, int M, Real z, Real zs, bool bDBessel, Complex &p, Complex &dp) {    

  map<Real, vector<Complex> >::iterator i = memo.lower_bound(z);

  bool bNegativeZ = (z < 0.);
  z = fabs(z);

  if(i != memo.end()) {
    if(fabs(i->first - z) < MEMOEPS) {        
      eval(i->second,N,M,i->first,zs,p,dp);
    } else if(i != memo.begin()) {
      --i;
      if(i != memo.end()) {
        if(fabs(i->first - z) < MEMOEPS) {
          eval(i->second,N,M,i->first,zs,p,dp);
        }
      }
    }
  }    
  
  i = memo.insert(i, pair<Real,vector<Complex> >(z,vector<Complex>()));
  eval(i->second,N,M,z,zs,p,dp);

  if(bNegativeZ) {
    if((M+N)%2==0) {
      dp = -dp;     
    } else {
      p = -p;
    }
  }

  if(c) {
    if(N > rN0) {
      p *= c[rN0];
      dp *= c[rN0];
      if(bDBessel) {
        p *= rd[N];
        dp *= rd[N];
      }
    } else {
      if(bDBessel) {
        p *= d[N];
        dp *= d[N];
      } else {
        p *= c[N];
        dp *= c[N];
      }
    }
  }
  if(isnan(p) || isnan(dp)) {
    fprintf(stderr,"%p %g\n",c,real(c[rN0]));
  }

}

#undef MEMOEPS

#endif
