#include <stdio.h>
#include <stdlib.h>

#include "real.h"
#include <complex.h>
#include "fftwe.h"
#include <math.h>

#define index(a, b, c) (c) + ((nby2p1) * ((b) + ((n * (a)))))

void write_power (real dk, int n, fftwe_complex * complex_delta, real g1, real bin, real start, char * fname) {
	long ind;
	long * kmod;
	double ks, power;
	double Pk[1000], rk[1000], jpower[1000];
	long nyquist = n / 2;
	long nby2p1 = nyquist + 1;
	int np = (int) floor((0.5 * dk * (n - 1)) / bin) + 3;
	for (long i = 0; i <= np; i++) {
		rk[i] = 0.;
		Pk[i] = 0.;
		jpower[i] = 0.;
	}
	/* this will help extracting information from data */
	kmod = (long *) calloc (n, sizeof (long));
	if (kmod == NULL) {
		fprintf (stderr, "not enough memory\n");
		exit (2);
	}
	for (long i = 0; i < n; i++) {
		if (i <= nyquist)
			kmod[i] = i;
		else
			kmod[i] = -(n - i);
	}
	for (long i = 0; i < n; i++) {
		for (long j = 0; j < n; j++) {
			for (long k = 0; k < nby2p1; k++) {
				if (k == 0 || k == nyquist) {
					if (j > nyquist)
						continue;
					else if (j == 0 || j == nyquist) {
						if (i > nyquist || i == 0)
							continue;
					}
				}
				ks = sqrt(kmod[i] * kmod[i] + kmod[j] * kmod[j] + kmod[k] * kmod[k]);
				if (ks <= (double) nyquist) {
					ind = (long) (ks - start) / bin + 0.5 + 1;
					if (ind <= np) { 
						power = pow(cabs(complex_delta[index(i, j, k)]) * g1, 2);
						if (power != 0.) {
							rk[ind] += dk * ks;
							Pk[ind] += power;
							jpower[ind]++;
						}
					}
				}
			}
		}
	}
	char outdata[256];
	sprintf(outdata,"pk_%s.txt", fname);
	FILE *fpout;
	if((fpout = fopen(outdata, "w")) == NULL) {
		printf("I cannot open %s\n", outdata);
		exit(1);
	}
	fprintf(fpout,"# Pk_in\n#dk = %g\t ng = %d\n", dk, n);
	for (long i = 1; i <= np; i++) {
		rk[i] /= jpower[i];
		Pk[i] /= jpower[i];
		if (jpower[i] != 0)
			fprintf(fpout,"%g %g %g %d\n", rk[i], Pk[i], Pk[i] / sqrt(jpower[i]), jpower[i]);	
	}
	fclose(fpout);	
	return;
}
