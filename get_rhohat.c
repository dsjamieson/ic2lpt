#include "real.h"
#include <complex.h>
#include "fftwe.h"
#include <omp.h>


#define index(i,j,k) (i)*nby2p1xn + (j)*nby2p1 + (k)

int populate_rhohat(int n,int k,int j,int i,real dk,fftwe_complex *delta);
int enforce_hermitian(int n,fftwe_complex *array);

int get_rhohat(real dk,int n,fftwe_complex *rhohat) {
  // corresponds to ranphfld.
  /* 2004-02-21 ES: reorder assignment so that a given random seed
     can be trusted to produce the same large-scale power regardless of
     grid size.  See first few pages of 2nd logbook for diagram.
     i3, i2, i1: indices for looping
     The trick is in ordering i1, i2, and i3 so that they're in the order
     of plane, row, column.  
     To make a truncated grid, just change the limit on i3 here. 
   	nk: number of elements in k, p
	k, p: table for power spectrum
	n:  number of grid elements (in one dimension)
	dk: grid spacing in fourier space (= 2pi/boxsize) 
     if you want rhohat with "my" discrete fourier convention, then multiply
     the output of this routine by L^3/2. */
  int i1,i2,i3;
  long nby2p1,nby2p1xn;

  nby2p1 = n/2 + 1;
  nby2p1xn = nby2p1 * n;

#pragma omp parallel private(i1,i2,i3)
 {
#pragma omp for
  for (i1=0;i1<nby2p1;i1++) {
    for (i2=0;i2<n;i2++) {
      for (i3=0;i3<n;i3++) {
	rhohat[i3*n*nby2p1 + i2*nby2p1 + i1] = 0. + I*0.; } } }
 }

//#pragma omp for
  for(i3=0;i3<nby2p1;i3++) { // i3: indicator of how far out you've gone
    for (i2=0;i2<i3;i2++) { // i2: each level (except last)
      for (i1=0;i1<i3+1;i1++) { // i1: horizontal segment of cells, going right
        populate_rhohat(n,i2,i3,i1,dk,rhohat);
      }
      for (i1=0;i1<i3;i1++) { // i1: vertical segment of cells, going down
        populate_rhohat(n,i2,i1,i3,dk,rhohat);
      }
    }
    // populate i^th level, which has to be done completely.
    // i3 now takes the role of the last level, previously like i2
    for (i2=0;i2<i3+1;i2++) {
      for (i1=0;i1<i3+1;i1++) {
        populate_rhohat(n,i3,i2,i1,dk,rhohat);
      }
    }
  }

  // enforce complex conjugate constraints
  enforce_hermitian(n,rhohat);

  // temp: zero-out the lowest frequency modes
/*
  for (i1=-1;i1<=1;i1++) {
  for (i2=-1;i2<=1;i2++) {
  for (i3=0;i3<=1;i3++) { // i3 > 0 b/c of CC redundancy
  rhohat[index((i1+n)%n,(i2+n)%n,(i3+n)%n)] = 0.;
  } } }
*/

  // temp: zero-out everything
/*
  for (i1=0;i1<n;i1++) {
  for (i2=0;i2<n;i2++) {
  for (i3=0;i3<n/2+1;i3++) { // i3 > 0 b/c of CC redundancy
  rhohat[index(i1,i2,i3)] = 0.;
  } } }
  rhohat[index(1,0,0)] = -1000.;
  rhohat[index(-1+n,0,0)] = -1000.; // i think unnecessary, gets set ...
  rhohat[index(0,1,0)] = -1000.;
  rhohat[index(0,-1+n,0)] = -1000.; // ... in enforce_hermtian
  rhohat[index(0,0,1)] = -1000.;
*/



  return 0;
}
