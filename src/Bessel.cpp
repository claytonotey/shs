#include "MathUtils.h"
#include "Bessel.h"

extern "C" {

  FORTRAN_RET_T
  FORTRAN(zbesj, ZBESJ) (double &zr,
                         double &zi,
                         double &nu,
                         int &kode,
                         int &n,
                         double *cyr,
                         double *cyi,
                         int &nz,
                         int &ierr);


  FORTRAN_RET_T
  FORTRAN(zbesh, ZBESH) (double &zr,
                         double &zi,
                         double &nu,
                         int &kode,
                         int &m,
                         int &n,
                         double *cyr,
                         double *cyi,
                         int &nz,
                         int &ierr);

  FORTRAN_RET_T
  FORTRAN(zbesk, ZBESK) (double &zr,
                         double &zi,
                         double &nu,
                         int &kode,
                         int &n,
                         double *cyr,
                         double *cyi,
                         int &nz,
                         int &ierr);
}

void BesselJ(const Complex &z,
             double nu,
             int n,
             Complex *j)
{
  double *jr = new double[n];
  double *ji = new double[n];
  int nz, ierr;

  double zr = real(z);
  double zi = imag(z);
  int kode = 1;
  FORTRAN(zbesj,ZBESJ)(zr,zi,nu,kode,n,jr,ji,nz,ierr);
  if(ierr) {
    fprintf(stderr,"zbesj returned with error %d\n",ierr);
    abort();
  }

  for(int i=0; i<n; i++) {
    j[i] = Complex(jr[i],ji[i]);
  }

  delete [] jr;
  delete [] ji;
}


void BesselK(const Complex &z,
             double nu,
             int n,
             Complex *j)
{
  double *jr = new double[n];
  double *ji = new double[n];
  int nz, ierr;

  double zr = real(z);
  double zi = imag(z);
  int kode = 1;
  FORTRAN(zbesk,ZBESK)(zr,zi,nu,kode,n,jr,ji,nz,ierr);
  if(ierr) {
    fprintf(stderr,"zbesj returned with error %d\n",ierr);
    abort();
  }

  for(int i=0; i<n; i++) {
    j[i] = Complex(jr[i],ji[i]);
  }

  delete [] jr;
  delete [] ji;
}

void SphericalBesselJ(const Complex &z,
                      Real nu,
                      int n,
                      Complex *j,
                      Complex *dj)
{
  Complex *J = new Complex[n+1];
  BesselJ(z,nu+0.5,n+1,J);

  Complex c = sqrt((Real)PI/((Real)2.*z));
  for(int i=0; i<n; i++) {
    j[i] = c * J[i];
    dj[i] = c * (-J[i+1] + (nu+(Real)(i+1))/z * J[i]); 
  }

  delete [] J;
}
                         
void BesselH(const Complex &z,
             double nu,
             int n,
             Complex *h)
{
  double *hr = new double[n];
  double *hi = new double[n];
  int nz, ierr;

  double zr = real(z);
  double zi = imag(z);

  int kode = 1;
  int m = 1;

  FORTRAN(zbesh,ZBESH)(zr,zi,nu,kode,m,n,hr,hi,nz,ierr);
  if(ierr) {
    fprintf(stderr,"zbesh returned with error %d\n",ierr);
    abort();
  }
    
  for(int i=0; i<n; i++) {
    h[i] = Complex(hr[i],hi[i]);
  }

  delete [] hr;
  delete [] hi;
}

void SphericalBesselH(const Complex &z,
                      Real nu,
                      int n,
                      Complex *h,
                      Complex *dh)
{
  Complex *H = new Complex[n+1];
  BesselH(z,nu+0.5,n+1,H);

  Complex c = sqrt((Real)PI/((Real)2.*z));
  for(int i=0; i<n; i++) {
    h[i] = c * H[i];
    dh[i] = c * (-H[i+1] + (nu+(Real)(i+1))/z * H[i]); 
  }

  delete [] H;
}
