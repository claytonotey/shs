#include "quadgk715.h"

static const Real xgk[15] = 
{
  0.0042723144395936940576063989283284, 
  0.025446043828620756865888097308925, 
  0.067567788320115451661251881887438, 
  0.12923440720030276995800022632466,
  0.20695638226615442611944217787823,
  0.29707742431130140792205907018797,
  0.3961075224960507457083735971537,
  0.5,
  0.6038924775039492542916264028463,
  0.7029225756886985365667896985542,
  0.79304361773384557388055782212177,
  0.87076559279969723004199977367534,
  0.93243221167988454833874811811256,
  0.97455395617137918762296067143325,
  0.99572768556040625043124236981384 
};
 
/* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule. 
   xgk[0], xgk[2], ... abscissae to optimally extend the 7-point gauss rule */
   
static const Real wgk15[15] =    /* weights of the 15-point kronrod rule */
{
  0.022935322010529224963732008058970,
  0.063092092629978553290700663189204,
  0.104790010322250183839876322541518,
  0.140653259715525918745189590510238,
  0.169004726639267902826583426598550,
  0.190350578064785409913256402421014,
  0.204432940075298892414161999234649,
  0.209482141084727828012999174891714,
  0.204432940075298892414161999234649,
  0.190350578064785409913256402421014,
  0.169004726639267902826583426598550,
  0.140653259715525918745189590510238,
  0.104790010322250183839876322541518,
  0.063092092629978553290700663189204,
  0.022935322010529224963732008058970
};

static const Real wg7[15] =     /* weights of the 7-point gauss rule */
{
  0,
  0.129484966168869693270611432679082,
  0,
  0.279705391489276667901467771423780,
  0,
  0.381830050505118944950369775488975,
  0,
  0.417959183673469387755102040816327,
  0,
  0.381830050505118944950369775488975,
  0,
  0.279705391489276667901467771423780,
  0,
  0.129484966168869693270611432679082,
  0,
};

#define MAXINTERVALS 4096
#define MAXPOINTS 61440
#define MAXITERATIONS 32

inline Real xx(Real u, Real A, Real B)
{
  return A + (B-A) * u * u;
}

inline Complex yy(Complex y, Real u)
{
  return 2. * y * u;
}

Complex quadgk715(Function *F, Real A, Real B, Real reltol, Real abstol, int N, bool bCV)
{
  Real a[MAXINTERVALS];
  Real b[MAXINTERVALS];
  Real u[MAXPOINTS];
  Real x[MAXPOINTS];
  Real yr[MAXPOINTS];
  Real yi[MAXPOINTS];
  Complex qi[MAXINTERVALS];
  Complex dqi[MAXINTERVALS];

  int stack[MAXINTERVALS];
  int nstack = 0;

  int top = 0;
  int next[MAXINTERVALS];
  int prev[MAXINTERVALS];
  int nlist = 0;

  Real D = B-A;
  Real D2 = D*D;
  Real tol;

  prev[top] = -1;
  for(int i=0;i<N;i++) {    
    a[i] = ((Real)i/(Real)N);
    b[i] = ((Real)(i+1)/(Real)N);
    if(i<N-1) {
      next[i] = i+1;
      prev[i+1] = i;
    } else {
      next[i] = -1;      
    }
  }
  nlist = N;

  Complex qok(0.);
  Complex dqok(0.);

  bool bMaxIntervals = false;

  int iterations = 0;
  while(true) {    
    int nx = nlist * 15;
    nx = 0;
    for(int i=top; i>=0; i = next[i]) {
      Real ai = a[i];
      Real di = b[i] - ai;      
      for(int j=0; j<15; j++) {
        u[nx] = ai + di * xgk[j];
        if(bCV)
          x[nx] = xx(u[nx],A,B);
        else
          x[nx] = u[nx];
        nx++;
      }
    }

    F->f(nx,x,yr,yi);

    nx = 0;
    Complex q = qok;
    
    for(int i=top; i>=0; i = next[i]) {
      Real di = 0.5 * (b[i] - a[i]);

      Complex q15(0.);
      Complex q7(0.);
      for(int j=0; j<15; j++) {
        Complex y(yr[nx],yi?yi[nx]:0.);
        if(bCV)
          y = yy(y,u[nx]);
        q15 += wgk15[j] * y;
        if(j%2==1) q7 += wg7[j] * y;
        nx++;
      }
      
      q += q15 * di;
      qi[i] = q15;
      dqi[i] = (q15 - q7);
    }

    tol = max(abstol,reltol*norm(q));

    nx = 0;
    Real dq = 0.;

    for(int i=top; i>=0;) {
      Real di = 0.5 * (b[i] - a[i]);
      Real ndqi = norm(dqi[i]);
      if(ndqi < tol) {
        qok += qi[i] * di;
        dqok += dqi[i] * di;
        if(bMaxIntervals) {
          i = next[i];
          continue;
        }
        stack[nstack++] = i;
        int p = prev[i];
        int n = next[i];
        if(p>=0) next[p] = n;
        if(n>=0) prev[n] = p;
        if(i==top) top = n;
        i = n;
        nlist--;
      } else {
        Real mi = 0.5 * (b[i] + a[i]);
        dq += ndqi * di * di;
        if(bMaxIntervals) {
          i = next[i];
          continue;
        }
        int ii;
        if(nstack) {
          nlist++;
          ii = stack[--nstack];
        } else {
          ii = nlist++;
          if(nlist >= MAXINTERVALS) {
            bMaxIntervals = true;            
            i = next[i];
            continue;
          }
        }
        int n = next[i];
        if(n>=0) prev[n] = ii;
        next[ii] = n;
        prev[ii] = i;
        next[i] = ii;
        a[ii] = mi;
        b[ii] = b[i];
        b[i] = mi;
        i = n;
      }
    }
    dq += norm(dqok);

    if(iterations++ > MAXITERATIONS || bMaxIntervals || dq < tol) {
      return D * q;
    }
  }
}
