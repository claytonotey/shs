#include "Util.h"
#include "SphereHalfspaceMPI.h"
#include "SI.h"
#include "quadgk715.h"
#include "zpmn.h"
#include "Bessel.h"
#include "MathUtils.h"
#include "Matrix.h"
#include <time.h>
#include "SparseMatrix.h"
#include "ConjugateGradient.h"
#include "mpi.h"
#include "pthread.h"

#include <queue>
using namespace std;

static bool logProfile = false;
#define MPI_REALREAL MPI_LONG_DOUBLE

/*
master - queues fjobs which send mjobs to servants
servants - process mjobs and send results to master
MPI heap = nproc fjobs (nproc includes the master proc, which has a master and a servant)
fjob = [[real eh] [real I[nz]]]
*/

class LinearSystem {
public:
  LinearSystem() {}
  virtual ~LinearSystem() {}
  virtual void insert(int i, int j, Complex &a)=0;
  virtual void end()=0;
  virtual Real solve(ComplexDoubleDenseMatrix &B)=0;
};

class DenseLinearSystem : public LinearSystem {
public:
  DenseLinearSystem(int n) : def(n,n,0) {}
  ~DenseLinearSystem() {}

  void insert(int i, int j, Complex &a) {
    def(i,j) = a;
  }

  void end() {
    def.LU();
  }  
  
  Real solve(ComplexDoubleDenseMatrix &B) {
    solveLU(def,B);
    return 0.;
  }
  
  ComplexDoubleDenseMatrix def;
};

class ComplexSparseMatrixVectorMultiply : public MatrixVectorMultiply<ComplexDouble> {
public:
  SparseMatrix<ComplexDouble> &A;
  ComplexSparseMatrixVectorMultiply(SparseMatrix<ComplexDouble> &A_) : A(A_) {}
  void multiply(Array<ComplexDouble> &x, Array<ComplexDouble> &y) {
    y.zero();
    A.MatrixVectorMultiply(x,1,y,1);
  }
};

class SparseLinearSystem : public LinearSystem {
public:
  SparseLinearSystem(int n, Real _eps, int _iters) : def(n,n), diag(n), eps(_eps), iters(_iters) {}
  ~SparseLinearSystem() {}

  void insert(int i, int j, Complex &a) {
    def.insert(i,j,(ComplexDouble)a);
    if(i == j) {
      diag[i] = (ComplexDouble)a;
    }
  }

  void end() {
    def.end();
  }  
  
  Real solve(ComplexDoubleDenseMatrix &B) {
    DiagonalPreconditioner<ComplexDouble,ComplexDouble> P(diag);
    ComplexSparseMatrixVectorMultiply mvm(def);
    Array<ComplexDouble> x(B.m);
    Array<ComplexDouble> b(B.m);
    
    double maxres = 0.;
    for(int j=0; j<B.n; j++) {
      B.writeColumn(x,j);
      B.writeColumn(b,j);
      double res = conjugateGradientSquared<double,ComplexDouble,ComplexDouble>(mvm,P,x,b,eps,iters);
      maxres = max(res,maxres);
      B.readColumn(x,j);
    }
    return maxres;
  }
  
  SparseMatrix<ComplexDouble> def;
  Array<ComplexDouble> diag;
  Real eps;
  int iters;
};

class tabubtba_integrand : public Function<Complex> {
public:
  tabubtba_integrand(bool bVar, int p, int q, int m, int l, int v, Real kf, Complex &kb,
                     ZPMN &zpmnl_, ZPMN &zpmnv_) : zpmnl(zpmnl_), zpmnv(zpmnv_) {
    this->bVar = bVar;
    this->p = p;
    this->q = q;
    this->m = m;
    this->l = l;
    this->v = v;
    this->kf = kf;
    this->kb = kb;
  }
  
  void f(int nx, Real *x, Complex *y) {

    for(int i=0; i<nx; i++) {
      Real z;
      Real beta, beta2;
      Complex cosalpha;
      Complex dx;

      if(bVar) {
        beta = sin(x[i]);
        beta2 = beta * beta;
        z = sqrt(max((Real)0.,1. - beta2));
        cosalpha = z;
        dx = beta;
      } else {
        beta = x[i];
        beta2 = beta * beta;
        z = sqrt(max(Real(0.),beta2 - 1.));
        cosalpha = Complex((Real)0.,z);
        dx = beta / cosalpha;
      }
      
      Real lambda = beta * kf;
      Complex kb2 = kb * kb;
      Real kf2 = kf * kf;
      Complex hf = cosalpha * kf;
      Complex hb = sqrt(kb2 - lambda * lambda);
      
      if(imag(hb) < 0.) {
        hb = -hb;
      }
      
      Complex Pl, dPl, Pv, dPv;
      
      zpmnl.eval(l,m,z,beta,p==1,Pl,dPl);
      zpmnv.eval(v,m,-z,beta,q==1,Pv,dPv);

      Complex c = Itothe(m-v+1) * dPv;
      Complex a = Itothe(m-v-1) * Pv * ((Real)m / beta2);
      
      Complex AB1 = Itothe(l-m-1) * dPl * beta2;
      Complex AB2 = Itothe(l-m+1) * Pl * (Real)m;

      Complex ABc;
      Complex ABa;

      if((p+q)%2==0) {
        ABc = AB1;
        ABa = AB2;
      } else {
        ABc = AB2;
        ABa = AB1;
      }
      
      Complex u1 = (hb - hf) / (hf + hb);
      Complex u2 = (hb * kf2 - hf * kb2) / (hb * kf2 + hf * kb2);
      
      Complex uvc;
      Complex uva;

      if(q==0) {
        uvc = u1;
        uva = u2;
      } else {
        uvc = u2;
        uva = u1;
      }

      y[i] = (Real)2. * dx * (c * uvc * ABc + a * uva * ABa);
    }
  }
  
  int bVar,p,q,m,l,v;
  Real kf;
  Complex kb;
  ZPMN &zpmnl;
  ZPMN &zpmnv;
};
  
Complex tabubtba(int p, int q, int m, int l, int v, Real kf, Complex &kb, Real d,
                 ZPMN &zpmnla, ZPMN &zpmnlb,
                 ZPMN &zpmnva, ZPMN &zpmnvb,
                 Real maxmaxbeta,
                 quadgk715<Complex> &quad, Real reltol, Real abstol) {
  
  tabubtba_integrand f1(true,p,q,m,l,v,kf,kb,zpmnla,zpmnva);
  tabubtba_integrand f2(false,p,q,m,l,v,kf,kb,zpmnlb,zpmnvb);

  Real minalpha = 1e-5;
  Real maxbeta = (Real)(8 + l + v)/(kf*d);
  Real scale = pow(4.,-max(0.,floor(log2(maxmaxbeta / maxbeta) / 2.)));
  int betaN = (scale * maxmaxbeta > 2. * maxbeta)?4:2;
  maxbeta = scale * maxmaxbeta + 1.;
  Real c = (Real)(l * (l+1)) * (Real)(v * (v+1));
  Real atol = abstol * c;
  Complex r1 = quad.integrate(&f1, minalpha, PIOVER2, reltol, atol, 2, false);
  Complex r2 = quad.integrate(&f2, 1., maxbeta, reltol, atol, betaN, true);
  
  Complex r = r1 + r2;

  r /= sqrt(c);

  return r;
}


class Stba_integrand : public Function<Complex> {
public:
  Stba_integrand(bool bVar, int p, int q, int m, int l, int v, Real kf, Complex &kb,
                 ZPMN &zpmnl_, ZPMN &zpmnv_) : zpmnl(zpmnl_), zpmnv(zpmnv_) {
    this->bVar = bVar;
    this->p = p;
    this->q = q;
    this->m = m;
    this->l = l;
    this->v = v;
    this->kf = kf;
    this->kb = kb;
  }
  
  void f(int nx, Real *x, Complex *y) {

    for(int i=0; i<nx; i++) {
      Real beta, beta2;
      Real z;
      Complex cosalpha;
      Complex dx;

      if(bVar) {
        beta = sin(x[i]);
        beta2 = beta * beta;
        z = sqrt(max((Real)0.,1. - beta2));
        cosalpha = z;
        dx = beta;
      } else {
        beta = x[i];
        beta2 = beta * beta;
        z = sqrt(max((Real)0.,beta2 - 1.));
        cosalpha = Complex(0.,z);
        dx = beta / cosalpha;
      }

      Real lambda = beta * kf;
      Complex kb2 = kb * kb;
      Real kf2 = kf * kf;
      Complex hf = cosalpha * kf;
      Complex hb = sqrt(kb2 - lambda * lambda);
      
      if(imag(hb) < 0.) {
        hb = -hb;
      }
      
      Complex Pl, dPl, Pv, dPv;

      zpmnl.eval(l,m,-z,beta,p==1,Pl,dPl);
      zpmnv.eval(v,m,-z,beta,q==1,Pv,dPv);
      
      Complex lca1 = Itothe(m-l+1) * dPl * beta;
      Complex lca2 = Itothe(m-l-1) * Pl * ((Real)m / beta);

      Complex vca1 = Itothe(m-v+1) * dPv * beta;
      Complex vca2 = Itothe(m-v-1) * Pv * ((Real)m / beta);

      Complex lc, la, vc, va;
      
      if(p==0) {
        lc = lca1;
        la = lca2;
      } else {
        lc = lca2;
        la = lca1;
      }

      if(q==0) {
        vc = vca1;
        va = vca2;
      } else {
        vc = vca2;
        va = vca1;
      }

      Complex uc = (hb - hf) / (hf + hb);
      Complex ua = (hb * kf2 - hf * kb2) / (hb * kf2 + hf * kb2);
      
      Complex SMN = (dx * lc * ((Real)1. - uc)) * conj(vc * ((Real)1. + uc));
      Complex SNM = (la * ((Real)1. + ua)) * conj(dx * va * ((Real)1. - ua));

      y[i] = SMN + SNM;
    }
  }
  
  int bVar,p,q,m,l,v;
  Real kf;
  Complex kb;
  ZPMN &zpmnl;
  ZPMN &zpmnv;
};
  
Complex Stba(int p, int q, int m, int l, int v, Real kf, Complex &kb, Real d,
             ZPMN &zpmnla, ZPMN &zpmnlb,
             ZPMN &zpmnva, ZPMN &zpmnvb,
             Real maxmaxbeta,
             quadgk715<Complex> &quad, Real reltol, Real abstol) {

  Stba_integrand f1(true,p,q,m,l,v,kf,kb,zpmnla,zpmnva);
  Stba_integrand f2(false,p,q,m,l,v,kf,kb,zpmnlb,zpmnvb);

  Real minalpha = 1e-5;
  Real maxbeta = (Real)(8 + l + v)/(kf*d);
  Real scale = pow(4.,-max(0.,floor(log2(maxmaxbeta / maxbeta) / 2.)));
  int betaN = (scale * maxmaxbeta > 2. * maxbeta)?4:2;
  maxbeta = scale * maxmaxbeta + (Real)1.;
  Real c = (Real)(l * (l+1)) * (Real)(v * (v+1));
  Real atol = abstol * c;
  Complex r1 = quad.integrate(&f1, minalpha, PIOVER2, reltol, atol, 2, false);
  Complex r2 = quad.integrate(&f2, (Real)1., maxbeta, reltol, atol, betaN, true);
  Complex r = r1 + r2;

  r /= sqrt(c);
   
  return r;
}

Real SphereHalfSpaceEH(Real a, Real gap, Real omega, const Complex &eps1, const Complex &eps2, int m, int nmax, Real eta) 
{  
  char filename[32];
  sprintf(filename,"eh_%.5f_%d",omega * SI::hbar / SI::eV,m);
  FILE *log = fopen(filename,"w");

  time_t start = clock();
  Real elapsed;

  Real d = a + gap;  
  Real kf = omega / SI::c;
  Complex ka = omega / SI::c * sqrt(eps1);
  Complex kb = omega / SI::c * sqrt(eps2);
  Complex kaa = ka * a;
  Real kfa = kf * a;
  Real kfd = kf * d;

  Real reltol = norm(eta);
  Real abstol1 = 1e-16;
  Real abstol2 = Square(1e-9/(Real)nmax);
  int maxTinyCount = 6;
  Real tiny1 = 1e-16;
  Real tiny2 = Square(1e-7/(Real)nmax);
  int sparseThresh = 256;
  Real sparseEps = Square(1e-4/(Real)nmax);
  Real sparseIters = 2*nmax;
  Real betamax = (Real)(8 + nmax + nmax) / kfd;

  quadgk715<Complex> quad;


  Array<Complex> sjaa(nmax+1);
  Array<Complex> hajaa(nmax+1);
  Array<Complex> uan(nmax+1);
  Array<Complex> uad(nmax+1);
  Array<Complex> ua(nmax+1);
  Array<Complex> van(nmax+1);
  Array<Complex> vad(nmax+1);
  Array<Complex> va(nmax+1);

  int nmaxfa;
  int nmaxaa;
  for(nmaxfa = 1; nmaxfa < nmax; nmaxfa++) {
    if(isBesselAsymptotic(kfa,nmaxfa)) break;
  }
  for(nmaxaa = 1; nmaxaa < nmax; nmaxaa++) {
   if(isBesselAsymptotic(kaa,nmaxaa)) break;
  }

  SphericalBessel<Real,Complex> jfa(kfa,1,nmax,nmaxfa);
  SphericalBessel<Real,Complex> jaa(kaa,1,nmax,nmaxaa);
  SphericalBessel<Real,Complex> ihfa(kfa,3,nmax,nmaxfa,true);
  SphericalBessel<Real,Complex> ihaa(kaa,3,nmax,nmaxaa,true);

  for(int l=1; l<=nmax; l++) {
    if(l <= nmaxaa) {
      sjaa[l] = jaa.z[l] / sqrt(norm(jaa.z[l]));
      hajaa[l] = sqrt(norm(jaa.z[l])) / ihaa.z[l];
    }
  }
  SphericalBesselSignJ<Real,Complex>(kaa,nmax,nmaxaa,jaa.z[nmaxaa],sjaa);
  SphericalBesselHAbsJ<Real,Complex>(kaa,nmax,nmaxaa,jaa.z[nmaxaa],(Real)1. / ihaa.z[nmaxaa],hajaa);

  for(int l=1; l<=nmax; l++) {
    uan[l] = kaa * jaa.d[l] - kfa * jfa.d[l];
    uad[l] = kaa * jaa.d[l] - kfa / ihfa.d[l];
    ua[l] = uan[l] / uad[l];

    van[l] = kaa - kfa * jaa.d[l] / jfa.d[l];
    vad[l] = kaa - kfa * jaa.d[l] * ihfa.d[l];
    va[l] = van[l] / vad[l];
  }

  Real EHtotal = 0.;

  ZPMN zpmnjr(kfd,jfa.z,jfa.x,jfa.r,jfa.d,nmaxfa,true);
  ZPMN zpmnjz(kfd,jfa.z,jfa.x,jfa.r,jfa.d,nmaxfa,false);
  ZPMN zpmnhr(kfd,ihfa.z,ihfa.x,ihfa.r,ihfa.d,nmaxfa,true);
  ZPMN zpmnhz(kfd,ihfa.z,ihfa.x,ihfa.r,ihfa.d,nmaxfa,false);
    
  ZPMN *zpr1[2]; zpr1[0] = &zpmnjr; zpr1[1] = &zpmnjr;
  ZPMN *zpz1[2]; zpz1[0] = &zpmnjz; zpz1[1] = &zpmnjz;
  ZPMN *zqr1[2]; zqr1[0] = &zpmnhr; zqr1[1] = &zpmnhr;
  ZPMN *zqz1[2]; zqz1[0] = &zpmnhz; zqz1[1] = &zpmnhz;
  
  ZPMN *zpr2[2]; zpr2[0] = &zpmnhr; zpr2[1] = &zpmnhr;
  ZPMN *zpz2[2]; zpz2[0] = &zpmnhz; zpz2[1] = &zpmnhz;
  ZPMN *zqr2[2]; zqr2[0] = &zpmnhr; zqr2[1] = &zpmnhr;
  ZPMN *zqz2[2]; zqz2[0] = &zpmnhz; zqz2[1] = &zpmnhz;
  
  int nmin = max(1,m);
  int nn = nmax - nmin + 1;

  LinearSystem *Ap;
  if(nn > sparseThresh) {
    Ap = new SparseLinearSystem(2*nn,sparseEps,sparseIters);
  } else {
    Ap = new DenseLinearSystem(2*nn);
  }
      
  LinearSystem &A = *Ap;
  ComplexDoubleDenseMatrix CM(2*nn,nn,0);
  ComplexDoubleDenseMatrix CN(2*nn,nn,0);
  DenseMatrix<Real,Real> normCM(2*nn,nn);
  DenseMatrix<Real,Real> normCN(2*nn,nn);
  Array<Real> EH(nn); EH.zero();

  elapsed = ((Real)(clock()-start))/(Real)CLOCKS_PER_SEC;
  if(logProfile) fprintf(log,"init: %g\n",elapsed);
  start = clock();

  for(int n1=0; n1<nn; n1++) {
    int l = n1 + nmin;
    int v;
    
    Real xiM = imag(conj(kaa/kfa) * conj(jaa.d[l]));
    Real xiN = -imag(conj(kaa/kfa) * jaa.d[l]);
    
    CM(n1,n1) = sqrt(xiM) * Complex(0.,-1.) / (kfa*uad[l]*sjaa[l]);
    CN(n1+nn,n1) = sqrt(xiN) * Complex(0.,1.) / (kfa*vad[l]*sjaa[l]);
    
    int tinyCount;
    
    for(int p=0; p<2; p++) {
      int nnp = (p==0?0:nn);
      for(int q=0; q<2; q++) {
        int nnq = (q==0?0:nn);
        
        tinyCount = 0;          
        for(int n2=0; n2<nn-n1; n2++) {
          l = nmin + n1;
          v = nmin + n1 + n2;   
          Complex u = (p==0?-ua[l]:-va[l]);     
          Real nu = norm(u);
          Real atol = abstol1/(nu==0.?1.:nu);
          
          Complex aij = u * tabubtba(p,q,m,l,v,kf,kb,d,*(zpr1[p]),*(zpz1[p]),*(zqr1[q]),*(zqz1[q]),betamax,quad,reltol,atol);
          if(l == v && p == q) {
            aij += 1.;
          }
          A.insert(n1+nnp,n1+n2+nnq,aij);
          if(norm(aij) < tiny1) {
            tinyCount++;
            if(tinyCount > maxTinyCount) break;
          } else {
            tinyCount = 0;
          }
        }
        
        tinyCount = 0;
        for(int n2=1; n2<nn-n1; n2++) {
          l = nmin + n1 + n2;
          v = nmin + n1;        
          
          Complex u = (p==0?-ua[l]:-va[l]);     
          Real nu = norm(u);
          Real atol = abstol1/(nu==0.?1.:nu);
          
          Complex aij = u * tabubtba(p,q,m,l,v,kf,kb,d,*(zpr1[p]),*(zpz1[p]),*(zqr1[q]),*(zqz1[q]),betamax,quad,reltol,atol);
          A.insert(n1+n2+nnp,n1+nnq,aij);
          if(norm(aij) < tiny1) {
            tinyCount++;
            if(tinyCount > maxTinyCount) break;
          } else {
            tinyCount = 0;
          }
        }
      }
    }
  }

  elapsed = ((Real)(clock()-start))/(Real)CLOCKS_PER_SEC;
  if(logProfile) fprintf(log,"tabubtba: %g\n",elapsed);
  start = clock();

  A.end();

  /*
  cout << "CM = " << CM << "\n";
  cout << "CN = " << CN << "\n";

  cout << "[";
  for(int i=0; i<2*nn; i++) {
    cout << "[";
    for(int j=0; j<2*nn; j++) {
      Complex z = ((DenseLinearSystem*)&A)->def(i,j);
      cout << real(z) << "+" << imag(z) << "i ";
    }
    cout << "]\n";
  }
  cout << "];\n";
  */

  Real resM = A.solve(CM);
  Real resN = A.solve(CN);
  delete Ap;
     
  elapsed = ((Real)(clock()-start))/(Real)CLOCKS_PER_SEC;
  if(logProfile) fprintf(log,"solve: %g %g %g\n",elapsed,resM,resN);
  start = clock();


  elapsed = ((Real)(clock()-start))/(Real)CLOCKS_PER_SEC;
  //if(logProfile) fprintf(log,"vis: %g\n",elapsed);
  start = clock();
  Real elapsed2 = 0.;

  for(int i=0; i<2*nn; i++) {
    for(int j=0; j<nn; j++) {
      normCM(i,j) = norm(CM(i,j));
      normCN(i,j) = norm(CN(i,j));
    }
  }

  for(int n1=0; n1<nn; n1++) {
    int tinyCount12[2][2];
    int tinyCount21[2][2];
    Real normS120[2][2];
    Real normS210[2][2];
    
    for(int p=0; p<2; p++) {
      for(int q=0; q<2; q++) {
        tinyCount12[p][q] = 0;
        tinyCount21[p][q] = 0;
        normS120[p][q] = Infinity;
        normS210[p][q] = Infinity;
      }
    }
    
    for(int n2=0; n2<nn-n1; n2++) {
      
      int v1 = nmin + n1;
      int v2 = nmin + n1 + n2;
      
      Real maxnormMM = 0.;
      Real maxnormMN = 0.;
      Real maxnormNM = 0.;
      Real maxnormNN = 0.;

      Real start2 = clock();

      int im1 = v1-nmin;
      int im2 = v2-nmin;
      int in1 = v1-nmin+nn;
      int in2 = v2-nmin+nn;

      for(int n3=0; n3<nn; n3++) {
        
        Real normMM1 = normCM(im1,n3);
        Real normMN1 = normCM(in1,n3);
        Real normMM2 = normCM(im2,n3);
        Real normMN2 = normCM(in2,n3);
        
        Real normNM1 = normCN(im1,n3);
        Real normNN1 = normCN(in1,n3);
        Real normNM2 = normCN(im2,n3);
        Real normNN2 = normCN(in2,n3);
        
        Real x;
        
        x = normMM1 * normMM2; if(x > maxnormMM) maxnormMM = x;
        x = normNM1 * normNM2; if(x > maxnormMM) maxnormMM = x;

        x = normMM1 * normMN2; if(x > maxnormMN) maxnormMN = x;
        x = normNM1 * normNN2; if(x > maxnormMN) maxnormMN = x;

        x = normMN1 * normMM2; if(x > maxnormNM) maxnormNM = x;
        x = normNN1 * normNM2; if(x > maxnormNM) maxnormNM = x;

        x = normMN1 * normMN2; if(x > maxnormNN) maxnormNN = x;
        x = normNN1 * normNN2; if(x > maxnormNN) maxnormNN = x;
      }
      
      Real maxnorm[2][2];
      maxnorm[0][0] = maxnormMM;
      maxnorm[0][1] = max(maxnormMN,maxnormNM);
      maxnorm[1][0] = max(maxnormMN,maxnormNM);
      maxnorm[1][1] = maxnormNN;

      elapsed2 += ((Real)(clock()-start2))/(Real)CLOCKS_PER_SEC;
      
      for(int p=0; p<2; p++) {
        for(int q=0; q<2; q++) {
          
          Complex S12;
          Complex S21;
          bool b12 = false;
          bool b21 = false;
          
          Real atol = abstol2 / (maxnorm[p][q] == 0.?1.:maxnorm[p][q]);
          
          if(tinyCount12[p][q] <= maxTinyCount || maxnorm[p][q] * normS120[p][q] > tiny2) {
            b12 = true;
            S12 = Stba(p,q,m,v1,v2,kf,kb,d,*(zpr2[p]),*(zpz2[p]),*(zqr2[q]),*(zqz2[q]),betamax,quad,reltol,atol);
            Real normS = norm(S12);
            if(normS <= normS120[p][q]) {
              normS120[p][q] = normS;
              tinyCount12[p][q]++;
            } else {
              normS120[p][q] = Infinity;
              tinyCount12[p][q] = 0;
            }
          }
          
          if(v1 != v2) {
            if(tinyCount21[p][q] <= maxTinyCount || maxnorm[p][q] * normS210[p][q] > tiny2) {
              b21 = true;
              S21 = Stba(p,q,m,v2,v1,kf,kb,d,*(zpr2[p]),*(zpz2[p]),*(zqr2[q]),*(zqz2[q]),betamax,quad,reltol,atol);
              Real normS = norm(S21);
              if(normS <= normS210[p][q]) {
                normS210[p][q] = normS;
                tinyCount21[p][q]++;
              } else {
                normS210[p][q] = Infinity;
                tinyCount21[p][q] = 0;
              }
            }
          }
          
          if(b12 || b21) {
            int nnp = (p==0?0:nn);
            int nnq = (q==0?0:nn);
            for(int l=nmin; l<=nmax; l++) {
              Real EHmn = 0.;
              if(b12) {
                EHmn += real((Complex)CM(v1-nmin+nnp,l-nmin) * S12 * conj((Complex)CM(v2-nmin+nnq,l-nmin)));
                EHmn += real((Complex)CN(v1-nmin+nnp,l-nmin) * S12 * conj((Complex)CN(v2-nmin+nnq,l-nmin)));
              }
              if(b21) {
                EHmn += real((Complex)CM(v2-nmin+nnp,l-nmin) * S21 * conj((Complex)CM(v1-nmin+nnq,l-nmin)));
                EHmn += real((Complex)CN(v2-nmin+nnp,l-nmin) * S21 * conj((Complex)CN(v1-nmin+nnq,l-nmin)));
              }
              if(m > 0) EHmn *= 2.;
              EH[l-nmin] += EHmn;
              EHtotal += EHmn;
            }
          }
        }
      }
    }
  }

  elapsed = ((Real)(clock()-start))/(Real)CLOCKS_PER_SEC;
  if(logProfile) fprintf(log,"stba: %g/%g\n",elapsed2,elapsed);

  for(int l=nmin; l<=nmax; l++) {
    fprintf(log,"%d %g\n",l,EH[l-nmin]);
  }
  fclose(log);

  return EHtotal;
}

struct MJob {
  int m;
  int nmax;
  Real a;
  Real gap;
  Real omega;
  Real eps1r;
  Real eps1i;
  Real eps2r;
  Real eps2i;
  Real eta;
};

class FJob {
public:
  Real a;
  Real gap;
  Real feV;
  Complex eps1;
  Complex eps2;
  int nmax;
  Real omega;
  vector<Real> EH;
  vector<bool> returned;
  Real kappa0;
  Real kappa1;
  Real kappa2;
  Real eta;
  int m;
  //  ShapeMesh &mesh;
  //map<Point*,int,ltpoint> &zMap;
  //Array<Real> I;
  Real EHtotal;
  
  FJob(Real a_, Real gap_, Real feV_, Complex &eps1_, Complex &eps2_, 
       Real kappa0_,
       Real kappa1_,
       Real kappa2_,
       Real eta_
       //       ShapeMesh &mesh_, 
       //map<Point*,int,ltpoint> zMap_
       ) :
    a(a_),
    gap(gap_),
    feV(feV_),
    eps1(eps1_),
    eps2(eps2_),
    kappa0(kappa0_),
    kappa1(kappa1_),
    kappa2(kappa2_),
    eta(eta_),
    m(0),
    //mesh(mesh_),
    //zMap(zMap_),
    //I(zMap_.size(),0),
    EHtotal(0.){
    omega = feV * SI::eV / SI::hbar;
    Real kf = omega / SI::c;
    nmax = lrintf(kappa0 + kappa1 * a / gap + kappa2 * kf * (a + gap));
  }
    
  void getMJob(MJob *mjob) {
    mjob->m = m;
    mjob->nmax = nmax;
    mjob->eta = eta;
    mjob->a = a;
    mjob->gap = gap;
    mjob->omega = omega;
    mjob->eps1r = real(eps1);
    mjob->eps1i = imag(eps1);
    mjob->eps2r = real(eps2);
    mjob->eps2i = imag(eps2);

    if(m >= 0) {
      EH.resize(max((int)EH.size(),m+1));
      returned.resize(max((int)returned.size(),m+1));
      EH[m] = 0.;
      returned[m] = false;
      m++;
    }
    if(m >= nmax) {
      m = -1;
    }
  }

  bool isDone() {
    return (m<0);
  }

  void applyResult(int mm, Real eh, Real *I) {
    /*
      for(unsigned int i=0; i<zMap.size(); i++) {
      this->I[i] += I[i];
    }
    */

    EH[mm] = eh;
    returned[mm] = true;
    EHtotal += eh;
    int tinyCount = 0;
    int mmax = EH.size() - 1;
    for(int i=0; i<=mmax; i++) {
      if(!returned[i]) return;
      if(EH[i] * (Real)(i+1) < 1e-6 * EHtotal) {
        if(++tinyCount > 4) {
          m = -1;
          return;
        } 
      } else {
        tinyCount = 0;
      }
    }
  }

  Real getS() {
    return 2. / PI * Square(omega * a / SI::c) * EHtotal;
  }

  /*
  Array<Real> &getI()  {
    SHSVisualize :: finalize(I,omega,a);
    return I;
  }
  */

};

#define TAG_MJOB 1
#define TAG_EHI 2

class Servant {
public:
  
  int masterProc;
  MPI_Datatype &jobType;
  //map<Point*,int,ltpoint> &zMap;

  Servant(int masterProc_, MPI_Datatype &jobType_/*,map<Point*,int,ltpoint> &zMap_*/) : 
    masterProc(masterProc_),
    jobType(jobType_) {}
  //zMap(zMap_) {}
  
  void start() {
    MJob mjob;
    MPI_Status status;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    char filename[32];
    sprintf(filename,"log%d",rank);
    FILE *log = fopen(filename,"w");
    while(true) {
      time_t start = clock();
      Real elapsed;

      MPI_Recv(&mjob,1,jobType,masterProc,TAG_MJOB,MPI_COMM_WORLD,&status);

      elapsed = ((Real)(clock()-start))/(Real)CLOCKS_PER_SEC;
      fprintf(log,"recv: %g\n",elapsed);
      start = clock();

      if(mjob.m < 0) break;
      Complex eps1(mjob.eps1r,mjob.eps1i);
      Complex eps2(mjob.eps2r,mjob.eps2i);

      //      SHSVisualize vis(zMap);
      Real eh = SphereHalfSpaceEH(mjob.a,mjob.gap,mjob.omega,eps1,eps2,mjob.m,mjob.nmax,mjob.eta);

      elapsed = ((Real)(clock()-start))/(Real)CLOCKS_PER_SEC;
      fprintf(log,"eh: %g\n",elapsed);
      start = clock();

      fflush(log);

      //int size = zMap.size();
      int size = 0;
      Array<Real> buf(1);
      buf[0] = eh;
      //memcpy(buf+1,vis.I.p,size*sizeof(Real));
      MPI_Send(buf,1+size,MPI_REALREAL,masterProc,TAG_EHI,MPI_COMM_WORLD);    

      elapsed = ((Real)(clock()-start))/(Real)CLOCKS_PER_SEC;
      fprintf(log,"send: %g\n",elapsed);
    }
    fclose(log);
    fprintf(stdout,"servant %d exiting...\n",rank);
  }
  
};

class Master {
public:
  
  MPI_Datatype &jobType;
  vector<Real> &as;
  vector<Real> &gaps;
  vector<Real> &freqs;
  Real kappa0;
  Real kappa1;
  Real kappa2;
  Real eta;
  //  ShapeMesh &mesh;
  //  map<Point*,int,ltpoint> &zMap;
  Dielectric *diel1;
  Dielectric *diel2;

  queue<int> procQ;
  queue<int> fjobQ;
  int na;
  int nf;
  int ngaps;
  int nfjobs;
  int nprocs;
  int nrequests;
  Array<FJob*> fjobs;
  Array<int> nmjobs;
  Array<int> procJ;
  Array<int> procM;
  Array<Real> procResult;
  Array<MPI_Request> procRequest;
  int bufsize;

  Master(MPI_Datatype &jobType_, vector<Real> &as_, vector<Real> &gaps_, vector<Real> &freqs_, Real kappa0_, Real kappa1_, Real kappa2_, Real eta_, 
         //ShapeMesh &mesh_, 
         //map<Point*,int,ltpoint> &zMap_, 
         Dielectric *diel1_, Dielectric *diel2_) :
    jobType(jobType_), as(as_), gaps(gaps_), freqs(freqs_), kappa0(kappa0_), kappa1(kappa1_), kappa2(kappa2_), eta(eta_),/* mesh(mesh_), zMap(zMap_),*/ diel1(diel1_), diel2(diel2_)
  {
  }
  
  void start() {
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Status status;

    na = as.size();
    nf = freqs.size();
    ngaps = gaps.size();
    nfjobs = na*nf*ngaps;
    fjobs.resize(nfjobs);
    nmjobs.resize(nfjobs);
    procRequest.resize(nprocs);
    //bufsize = (1+zMap.size());
    bufsize = 1;
    procResult.resize(nprocs*bufsize);
    procJ.resize(nprocs);
    procM.resize(nprocs);
    nrequests = 0;

    FILE *out = fopen("output","w");

    int j = 0;
    for(int l=0; l<na; l++) {
      Real a = as[l];
      for(int k=0; k<ngaps; k++) {
        Real gap = gaps[k];
        for(int i=0; i<nf; i++) {
          Real feV = freqs[i];
          Complex eps1 = diel1->eps(feV);
          Complex eps2 = diel2->eps(feV);
          FJob *fjob = new FJob(a,gap,feV,eps1,eps2,kappa0,kappa1,kappa2,eta/*mesh,zMap*/);
          fjobQ.push(j);
          fjobs[j] = fjob;
          nmjobs[j] = 0;
          j++;
        }
      }
    }

    for(int i=0; i<nprocs; i++) {
      procQ.push(i);
    }

    while(nfjobs) {
      while(!procQ.empty() && !fjobQ.empty()) {
        int proc = procQ.front(); procQ.pop();
        FJob *fjob;

        int j;
        do {
          j = fjobQ.front(); fjobQ.pop();
          fjob = fjobs[j];
          if(fjob && fjob->isDone()) {
            fjob = NULL;
          }
        } while(fjob == NULL && !fjobQ.empty());

        if(fjob == NULL) break;

        MJob mjob;
        fjob->getMJob(&mjob);
        MPI_Ssend(&mjob,1,jobType,proc,TAG_MJOB,MPI_COMM_WORLD);
        MPI_Irecv(procResult+proc*bufsize,bufsize,MPI_REALREAL,proc,TAG_EHI,MPI_COMM_WORLD,procRequest+proc);
        procJ[proc] = j;
        procM[proc] = mjob.m;
        nmjobs[j]++;
        nrequests = min(nprocs,nrequests+1);
        if(!fjob->isDone()) {
          fjobQ.push(j);
        }
      }

      int procDone;
      MPI_Waitany(nrequests,procRequest,&procDone,&status);

      if(procDone != MPI_UNDEFINED) {
        int j = procJ[procDone];
        FJob *fjob = fjobs[j];
        fjob->applyResult(procM[procDone],procResult[procDone*bufsize],procResult+(procDone*bufsize)+1);
        nmjobs[j]--;
        if(fjob->isDone() && (nmjobs[j]==0)) {
          Real a = fjob->a;
          Real gap = fjob->gap;
          Real feV = fjob->feV;
          Real p = fjob->getS();
          fprintf(out,"%.15g %.15g %.15g %.15g\n",a,gap,feV,p);
          fflush(out);

          /*
          Array<Real> &I = fjob->getI();
          char filename[128];
          sprintf(filename,"I_%.4g_%.4g_%.4g",a,gap,feV);
          ofstream out(filename);
          out << "I = [\n";
          for(vector< Array<Triangle>* >::const_iterator ti = mesh.refinedTriangle.begin(); ti != mesh.refinedTriangle.end(); ++ti) {
            Array<Triangle> *triangle = *ti;
            for(int j = 0; j<triangle->n; j++) {
              Triangle *t = (*triangle) + j;
              Real i1 = I[zMap[t->p1]];
              Real i2 = I[zMap[t->p2]];
              Real i3 = I[zMap[t->p3]];
              out << "[" << *(t->p1) << " " << *(t->p2) << " " << *(t->p3) << " [" << i1 << " " << i2 << " " << i3 << "]];\n";
            }
          }
          out << "]\n";
          out.close();
          */


          delete fjob;
          fjobs[j] = NULL;
          nfjobs--;
        }
        procQ.push(procDone);
      }
    }

    for(int i=0; i<nprocs; i++) {
      MJob mjob;
      mjob.m = -1;
      MPI_Ssend(&mjob,1,jobType,i,TAG_MJOB,MPI_COMM_WORLD);
    }
    
    fprintf(stdout,"master %d exiting...\n",rank);
    fclose(out);
  }
};
  
void *masterThreadCB(void *data) {
  
  Master *master = (Master*)data;
  master->start();

  return NULL;
}

void SHSMPI(vector<Real> &as, vector<Real> &gaps, vector<Real> &freqs, Real kappa0, Real kappa1, Real kappa2, Real eta, /*const char *shapeFile*/ Dielectric *diel1, Dielectric *diel2)
{
  int rank;
  int provided;
  MPI_Init_thread(NULL,NULL,MPI_THREAD_MULTIPLE,&provided);
  if(provided != MPI_THREAD_MULTIPLE) {
    MPI_Finalize();
    abort();
  }
  
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  initReferenceCountingPointerLock();

  MPI_Datatype jobtype, primtypes[2]; 
  int blockcounts[2];
  MPI_Aint offsets[2], extent;

  offsets[0] = 0;
  primtypes[0] = MPI_INT;
  blockcounts[0] = 2;
  MPI_Type_extent(MPI_INT, &extent);
  offsets[1] = 4 * extent;

  primtypes[1] = MPI_REALREAL;
  blockcounts[1] = 8;
  MPI_Type_struct(2, blockcounts, offsets, primtypes, &jobtype);
  MPI_Type_commit(&jobtype);


  Dielectric *dielectric;
  /*
  ShapeParser parser;
  Shape shape;
  map<Point*,int,ltpoint> zMap;
  parser.parse(shapeFile,shape,&dielectric);
  ShapeMesh mesh(shape);

  // index the mesh points
  for(vector< Array<Point>* >::const_iterator i = mesh.refinedPoint.begin(); i != mesh.refinedPoint.end(); ++i) {
    Array<Point> *point = *i;
    for(int j = 0; j<point->n; j++) {
      Point *z = (*point) + j;
      zMap[z] = 0;
    }
  }
  
  int q = 0;
  for(map<Point*,int,ltpoint>::iterator zi = zMap.begin(); zi != zMap.end(); ++zi) {
    zMap[zi->first] = q++;
  }

  */

  if(rank == 0) {
    Master master(jobtype,as,gaps,freqs,kappa0,kappa1,kappa2,eta,/*mesh,zMap,*/diel1,diel2);
    pthread_t masterThread;
    pthread_create(&masterThread, NULL, masterThreadCB, &master);
    Servant servant(0,jobtype/*,zMap*/);
    servant.start();
    pthread_join(masterThread,NULL);
  } else {
    Servant servant(0,jobtype/*,zMap*/);
    servant.start();
  }

  MPI_Type_free(&jobtype);
  MPI_Finalize();
}
