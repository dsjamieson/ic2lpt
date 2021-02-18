#include "real.h"

#define index(a, b, c) (c) + ((twonby2p1) * ((b) + ((ng) * (a))))

int get_l2abc(int ng, real * x, real * y, real * z) {
	long twonby2p1 = 2 * (ng / 2 + 1);
	#pragma omp parallel 
	{
		long ind;
		#pragma omp for
		for (long i = 0; i < ng; i++) {
			for (long j = 0; j < ng; j++) {
				for (long k = 0; k < ng; k++) {
					ind = index(i, j, k);
					x[ind] = x[ind] * y[ind] + y[ind] * z[ind] + z[ind] * x[ind];
				}
			}
		}
	}
	return 0;
}

int get_l2def(int ng, real * x, real * y) {
	long twonby2p1 = 2 * (ng / 2 + 1);
	#pragma omp paralel
	{
		long ind;
		#pragma omp for
		for (long i = 0; i < ng; i++) {
			for (long j = 0; j < ng; j++) {
				for (long k = 0; k < ng; k++) {
					ind = index(i, j, k);
					x[ind] -= y[ind] * y[ind];
				}
			}
		}
	}
	return 0;
}
