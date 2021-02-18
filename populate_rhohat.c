#include <stdlib.h>
#include "real.h"
#include <complex.h>
#include "fftwe.h"
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define pi 3.1415926535897932384626

#define index(i, j, k) ((k) + ((nby2p1) * ((j) + ((n) * (i)))))

extern const gsl_rng_type  *  T;
extern gsl_rng  *  rng;

real power_interp(real kk);

int populate_rhohat(int n, int i, int j, int k, real dk, fftwe_complex  * delta, int seed) {
	long ind;
	real kk, power, expectation, gauss1, gauss2;
	long nby2p1 = n / 2 + 1;
 	kk = dk * sqrt(i * i + j * j + k * k);
	power = power_interp(kk);
	expectation = sqrt(power / 2.);
	ind = index(i, j, k);
	gauss1 = gsl_ran_gaussian(rng, 1.);
	gauss2 = gsl_ran_gaussian(rng, 1.);
	delta[ind] = expectation * gauss1 + I * expectation * gauss2;
	if (i) {
		ind = index(i - n, j, k);
		gauss1 = gsl_ran_gaussian(rng, 1.);
		gauss2 = gsl_ran_gaussian(rng, 1.);
		delta[ind] = expectation * gauss1 + I * expectation * gauss2;
	}
	if (j) { 
		ind = index(i, n - j, k);
		gauss1 = gsl_ran_gaussian(rng, 1.);
		gauss2 = gsl_ran_gaussian(rng, 1.);
		delta[ind] = expectation * gauss1 + I * expectation * gauss2;
	}
	if (i && j) {
		ind = index(i - n, j - n, k);
		gauss1 = gsl_ran_gaussian(rng, 1.);
		gauss2 = gsl_ran_gaussian(rng, 1.);
		delta[ind] = expectation * gauss1 + I * expectation * gauss2;
	}
	return 0;
}
