#include "real.h"
#include <complex.h>
#include "fftwe.h"
#include <math.h>
#include <omp.h>

#define pi 3.1415926535897932384626

int enforce_hermitian(int n,fftwe_complex *array);

/* Formerly poisson. Now: given an NxNxN array rhohat with spacing dk,
   take the derivative as such: (d1,d2) = (0,-1) : -ikx/k^2.  
	(d1,d2) = (0,0): kx^2/k^2. (d1,d2) = (0,1) = kx*ky/k^2, etc.
   Careful with signs. Div Psi = -delta, but this routine only does
   the inverse divergence = -ikx/k^2, so actually delta*(ikx/k^2) = Psi.
   the array k[3] will contain the 3 components of k.
   2004-12-14. Fixed a 2-hour bug: the indices i1 and i2 hadn't been
   wrapping around.  Interestingly, this meant that as long as I was
   using the discrete case there would be no ill effects. */
int greens(int n,real dk, fftwe_complex *rhohat, int d1, int d2) {

  real  dx;
  long nby2p1,nby2p1xn;

  nby2p1 = n/2+1;
  nby2p1xn = nby2p1*n;
  dx = 2*pi/(dk*n);
  
#pragma omp parallel
{
  long i1,i2,i3,index;
  long i1shift,i2shift;
  real k[3],k2;
  int continuous = 1;

#pragma omp for
  for (i1=0;i1<n;i1++) {
  i1shift = (i1 >= n/2) ? (n) : 0;
  k[0] = dk * (i1 - i1shift); // k_z
  for (i2=0;i2<n;i2++) {
  i2shift = (i2 >= n/2) ? (n) : 0;
  k[1] = dk * (i2 - i2shift); // k_y
  for (i3=0;i3<nby2p1;i3++) {
    k[2] = dk * i3; 	      // k_x
    index = i1*nby2p1xn + i2*nby2p1 + i3;
    if (continuous == 1) { // --- continuous ---
      k2 = -(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
      if (i1 || i2 || i3) {
	rhohat[index] *= (I * k[d1]/k2);
	if (d2 > -1) rhohat[index] *= (I * k[d2]);
      } else {
	rhohat[index] = 0.;
      }
    } else { // --- discrete ---
      k2 = 2.*(cos(k[0]*dx) + cos(k[1]*dx) + cos(k[2]*dx) - 3.) / (dx*dx);
      if (i1 || i2 || i3) {
	rhohat[index] *= (I*sin(k[d1]*dx)/dx / k2);
	if (d2 > -1) rhohat[index] *= (I*sin(k[d2]*dx)/dx);
      } else {
	rhohat[index] = 0.;
      }
    }
  }
  }
  }
} //end omp region
 
// --- is there a better way to enforce hermitianness?
  enforce_hermitian(n,rhohat);

  return 0;
}
