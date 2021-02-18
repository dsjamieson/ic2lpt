#include "real.h"
#include <complex.h>
#include "fftwe.h"
#include <omp.h>


#define index(i, j, k) ((k) + ((nby2p1) * ((j) + ((n) * (i)))))

int populate_rhohat(int n, int k, int j, int i, real dk, fftwe_complex *delta); 
int enforce_hermitian(int n, fftwe_complex *array); 

int get_rhohat(real dk, int n, fftwe_complex *rhohat) {
	long nby2p1 = n / 2 + 1; 
	#pragma omp parallel for
	for (long i1 = 0; i1 < n; i1++) {
		for (long i2 = 0; i2 < n; i2++) {
			for (long i3 = 0; i3 < nby2p1; i3++) {
				rhohat[index(i1, i2, i3)] = 0. + I * 0.;	
			}
		}
	}
	for(long i3 = 0; i3 < nby2p1; i3++) { 
		for (long i2 = 0; i2 < i3; i2++) {
			for (long i1 = 0; i1 < i3+1; i1++) {
				populate_rhohat(n, i2, i3, i1, dk, rhohat); 
			}
			for (long i1 = 0; i1 < i3; i1++) {
				populate_rhohat(n, i2, i1, i3, dk, rhohat); 
			}
		}
		for (long i2 = 0; i2 < i3 + 1; i2++) {
			for (long i1 = 0; i1 < i3 + 1; i1++) {
				populate_rhohat(n, i3, i2, i1, dk, rhohat); 
			}
		}
	}
	enforce_hermitian(n, rhohat); 
	return 0; 
}
