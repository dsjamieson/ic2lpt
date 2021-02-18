#include "real.h"
#include <complex.h>
#include "fftwe.h"
#include <omp.h>

#define index(i, j, k) (i)*nby2p1xn + (j)*nby2p1 + (k)
#define sqrt2 1.414213562373

int enforce_hermitian(int n, fftwe_complex *array) {
	long nby2 = n/2;
	long nby2p1 = nby2 + 1;
	long nby2p1xn = nby2p1 * n;
	#pragma omp parallel
	{
		#pragma omp for 
		for (long k = 0; k < nby2 + 1; k += nby2) { 
			array[index(0, 0, k)] = sqrt2 * creal(array[index(0, 0, k)]);
			array[index(0, nby2, k)] = sqrt2 * creal(array[index(0, nby2, k)]);
			array[index(nby2, 0, k)] = sqrt2 * creal(array[index(nby2, 0, k)]);
			array[index(nby2, nby2, k)] = sqrt2 * creal(array[index(nby2, nby2, k)]);
			for (long i = 1; i < nby2;i ++) {
				array[index(0, n-i, k)] = conj(array[index(0, i, k)]);
				array[index(nby2, n-i, k)] = conj(array[index(nby2, i, k)]);
				array[index(n-i, 0, k)] = conj(array[index(i, 0, k)]);
				array[index(n-i, nby2, k)] = conj(array[index(i, nby2, k)]);
			}
			for (long j = 1; j < nby2; j++) {
				for (long i = 1; i < nby2; i++)
					array[index(n-j, n-i, k)] = conj(array[index(j, i, k)]);
				for (long i = nby2 + 1; i < n; i++)
					array[index(n-j, n-i, k)] = conj(array[index(j, i, k)]);
			}
		}
	} 
	return 0;
}
