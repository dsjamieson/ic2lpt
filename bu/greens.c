#include "real.h"
#include <complex.h>
#include "fftwe.h"
#include <math.h>
#include <omp.h>
#include <stdlib.h>

#define pi 3.1415926535897932384626
#define index(i, j, k) ((k) + ((nby2p1) * (j + ((n) * (i)))))

int enforce_hermitian(int n, fftwe_complex * array); 

int greens(int n, real dk, fftwe_complex * rhohat, int d1, int d2) {
	long nby2 = n / 2;
	long nby2p1 = nby2 + 1; 
	real dx = 2 * pi / (dk * n); 
	#pragma omp parallel
	{
		long ind; 
		long k[3];
		real k2; 
		int continuous = 1; 
		#pragma omp for
		for (long i1 = 0; i1 < n; i1++) {
			k[0] = (i1 > nby2) ? i1 - n : i1; 
			for (long i2 = 0; i2 < n; i2++) {
				k[1] = (i2 >= nby2) ? i2 - n : i2; 
				for (long i3 = 0; i3 < nby2p1; i3++) {
					k[2] = i3;
					ind = index(i1, i2, i3); 
					if (i1 || i2 || i3) {
						/*
						if (d1 != d2) {
							if (k[d1] == 0 || k[d1] == nby2) {
								rhohat[ind] = 0.;
								continue;
							}
						}
						*/
						k2 = -dk * ((double) (k[0] * k[0] + k[1] * k[1] + k[2] * k[2])); 
						rhohat[ind] *= I * ((double) k[d1]) / k2; 
						if (d2 > -1) {
							/*
							if (d1 != d2) {
								if (k[d2] == 0 || k[d2] == nby2) {
									rhohat[ind] = 0.;
									continue;
								}
							}
							*/
							rhohat[ind] *= I * dk * (double) k[d2]; 
						}
					}
					else
						rhohat[ind] = 0.;
				}
			}
		}
	}
	//rhohat[0] = 0.;
	enforce_hermitian(n, rhohat); 
	return 0; 
}
