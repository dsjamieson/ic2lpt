#include "real.h"
#include <complex.h>
#include "fftwe.h"
#include <omp.h>

#define index(i,j,k) (i)*nby2p1xn + (j)*nby2p1 + (k)
#define sqrt2 1.414213562373
/* enforce that 3-D array satisfies the hermitian constraints so that
   its Fourier transform is real. 
   Consider element [i,j,k] in the N x N x N/2+1 array.  Note that
   all elements with k > 0 are taken care of already by the FFTW data
   arrangement, with the exception of k = n/2 (Nyquist).  So we just have
   to modify the k=0 and k=n/2 planes.
   2004-10-02. Fixed indexing bug (changed k,i,j to i,j,k)
   2004-10-05. changed creal to cabs for those 8 elements constrained to be
   real. Is this correct?  Seems so, since the expectation value for 
   deltahat(0) should still equal P(0), right?
   2004-11-08. It seems it's best for delta(0) to be drawn from a Gaussian
   so that it's real, while still satisfying <delta> = 0 and <delta^2> = P(k).
   Since in populate_rhohat, expectation=sqrt(power/2.), we have to multiply
   creal by sqrt(2).  Equivalently, the mean of a Gaussian^2 is 1 but the 
   mean of a Rayleigh^2 is 2.
   2005-1-11. Note this routine won't work as expected for the wrong_xi demo,
   because this routine depends on the fact that you've already assigned
   amplitudes to every mode, including the DC offset, that are a factor
   of sqrt(2.) too low.
*/

int enforce_hermitian(int n,fftwe_complex *array) {
  long nby2,nby2p1,nby2p1xn;
  long i,j,k;

  nby2 = n/2;
  nby2p1 = nby2 + 1;
  nby2p1xn = nby2p1 * n;

#pragma omp parallel private(k,j,i)
 {
#pragma omp for 
  for (k=0;k<nby2+1;k+=nby2) { // do k=0 and k=Nyquist planes
  array[index(0,0,k)] = sqrt2*creal(array[index(0,0,k)]);
  array[index(0,nby2,k)] = sqrt2*creal(array[index(0,nby2,k)]);
  array[index(nby2,0,k)] = sqrt2*creal(array[index(nby2,0,k)]);
  array[index(nby2,nby2,k)] = sqrt2*creal(array[index(nby2,nby2,k)]);
  for (i=1;i<nby2;i++) {
    array[index(0,n-i,k)] = conj(array[index(0,i,k)]);
    array[index(nby2,n-i,k)] = conj(array[index(nby2,i,k)]);
    array[index(n-i,0,k)] = conj(array[index(i,0,k)]);
    array[index(n-i,nby2,k)] = conj(array[index(i,nby2,k)]);
  }
  for (j=1;j<nby2;j++) {
    for (i=1;i<nby2;i++) array[index(n-j,n-i,k)] = conj(array[index(j,i,k)]);
    for (i=nby2+1;i<n;i++) array[index(n-j,n-i,k)] = conj(array[index(j,i,k)]);
  }
  }
 } //end parallel
  

  return 0;
}
