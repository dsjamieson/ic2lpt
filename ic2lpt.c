#include "real.h"
#include <complex.h>
#include "fftwe.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "ic2lpt.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define MAXLINELEN 1000

void readIntParam(FILE * t_file, const char t_param_name[], int * t_value);
void readDoubleParam(FILE * t_file, const char t_param_name[], double * t_value);
void readStringParam(FILE * t_file, const char t_param_name[], char * t_value);

int nk;
long ranc_counter;
real * ktable, * ptable;
const gsl_rng_type * T;
gsl_rng * rng;
real dk, omega_m, omega_l, hubble;

int main(int argc, char *argv[]) {
	if (argc != 3) {
		printf("Usage: ic2lpt param_filename nthreads\n");
		return 0;
	}
	int seed, np, ng, nthreads, no_files;
	real boxsize, dnow, zthen, hthen, dthen, fthen, d2then, f2then;
	char pk_file[100], outprefix[100];
	long np3, ng3x, ng3k, il;
	real athen, g1, g2, gdot1, gdot2, sqrt_vol;
	fftwe_complex * bf0k, * bf1k, * bf2k, * bf3k, * bf4k;
	real * bf0x, * bf1x, * bf2x, * bf3x;
	fftwe_plan bf0k_plan, bf1k_plan, bf2k_plan, bf1x_plan, bf2x_plan;

   {
        FILE * f = fopen(argv[1],"r");
		if (f == NULL) { 
			fprintf(stderr, "Error, could not open parameter file %s for reading\n", argv[1]);
			exit(-1); 
		}
        char param_name[100];
        readIntParam(f, "NumPart1D", &np);
        readIntParam(f, "NumMesh1D", &ng);
        readIntParam(f, "Seed", &seed);
        readIntParam(f, "NumFiles", &no_files);
        readDoubleParam(f, "BoxLength", &boxsize);
        readDoubleParam(f, "HubbleParam", &hubble);
        readDoubleParam(f, "OmegaM", &omega_m);
        readDoubleParam(f, "D0", &dnow);
        readDoubleParam(f, "Zi", &zthen);
        readDoubleParam(f, "Hi", &hthen);
        readDoubleParam(f, "D1i", &dthen);
        readDoubleParam(f, "F1i", &fthen);
        readDoubleParam(f, "D2i", &d2then);
        readDoubleParam(f, "F2i", &f2then);
        readStringParam(f, "PowerFile", pk_file);
        readStringParam(f, "OutBase", outprefix);
    }
	nthreads = atoi(argv[2]);

	 omp_set_num_threads(nthreads);
	 if(fftw_init_threads() == 0){
			 printf("error with fftw_init_threads!\n");
			 exit(1);
	 }
	 fftw_plan_with_nthreads(nthreads);

	dk = 2. * pi / boxsize;
	{
		char line[MAXLINELEN];
		float ktemp, ttemp;
		FILE * f = fopen(pk_file,"r");
		if (f == NULL) { 
			perror(pk_file);
			exit(-1); 
		}
		nk = 0;
		while (fgets(line, MAXLINELEN, f) != NULL)
			nk++;
		rewind(f);
		ktable = (real *) malloc(nk * sizeof(real));
		if(ktable == NULL){
			printf("not enough memory!\n");
			exit(1);
		}
		ptable = (real *) malloc(nk * sizeof(real));
		if(ptable == NULL){
			printf("not enough memory!\n");
			exit(1);
		}
		for (int i = 0; i < nk; i++) {
			fgets(line, MAXLINELEN, f); // use fgets instead of fscanf to ignore 
			sscanf(line,"%f%f", &ktemp, &ttemp);	//	last three columns for TF file
			ktable[i] = ktemp;
			ptable[i] = ttemp;
		}
		fclose(f);
	}

	ng3x = ((long) ng) * ng * ng;
	ng3k = ((long) ng) * ng * (ng / 2 + 1);
	np3	= ((long) np) * np * np;

	xx = (real *) malloc(np3 * sizeof(real));
	if(xx==NULL){
			printf("not enough memory for particles %g!\n",(1.*np3*sizeof(real))/1.E9);
			exit(1);
	} 
	yy = (real *) malloc(np3 * sizeof(real));
	if(yy==NULL){
			printf("not enough memory for particles %g!\n",(1.*np3*sizeof(real))/1.E9);
			exit(1);
	}
	zz = (real *)malloc(np3*sizeof(real));
	if(zz==NULL){
			printf("not enough memory for particles %g!\n",(1.*np3*sizeof(real))/1.E9);
			exit(1);
	}
	vx = (real *) malloc(np3 * sizeof(real));
	if(vx==NULL){
			printf("not enough memory for particles %g!\n",(1.*np3*sizeof(real))/1.E9);
			exit(1);
	}
	vy = (real *) malloc(np3 * sizeof(real));
	if(vy==NULL){
			printf("not enough memory for particles %g!\n",(1.*np3*sizeof(real))/1.E9);
			exit(1);
	}
	vz = (real *) malloc(np3 * sizeof(real));
	if(vz==NULL){
			printf("not enough memory for particles %g!\n",(1.*np3*sizeof(real))/1.E9);
			exit(1);
	}
	bf0k = fftwe_malloc(ng3k * sizeof(fftwe_complex));
	if(bf0k==NULL){
			printf("not enough memory for bf0k!\n");
			exit(1);
	}
	bf1k = fftwe_malloc(ng3k * sizeof(fftwe_complex));
	if(bf1k==NULL){
			printf("not enough memory for bf1k!\n");
			exit(1);
	}
	bf2k = fftwe_malloc(ng3k * sizeof(fftwe_complex));
	if(bf2k==NULL){
			printf("not enough memory for bf2k!\n");
			exit(1);
	}
	bf3k = fftwe_malloc(ng3k * sizeof(fftwe_complex));
	if(bf3k == NULL){
			printf("not enough memoryfor bf3k!\n");
			exit(1);
	}
	bf0x = (real *) bf0k;
	bf1x = (real *) bf1k;
	bf2x = (real *) bf2k;
	bf3x = (real *) bf3k;
	bf0k_plan = fftwe_plan_dft_c2r_3d(ng, ng, ng, bf0k, bf0x, FFTW_ESTIMATE);
	bf1k_plan = fftwe_plan_dft_c2r_3d(ng, ng, ng, bf1k, bf1x, FFTW_ESTIMATE);
	bf2k_plan = fftwe_plan_dft_c2r_3d(ng, ng, ng, bf2k, bf2x, FFTW_ESTIMATE);
	bf1x_plan = fftwe_plan_dft_r2c_3d(ng, ng, ng, bf1x, bf1k, FFTW_ESTIMATE);
	bf2x_plan = fftwe_plan_dft_r2c_3d(ng, ng, ng, bf2x, bf2k, FFTW_ESTIMATE);

	sqrt_vol = pow(boxsize, 1.5);
	athen = 1. / (zthen + 1.);
	g1 = dthen / dnow;
	gdot1 = dthen * fthen * hthen * 100. * athen / dnow;
	g2 = d2then / dnow / dnow;
	gdot2 = d2then * f2then * hthen * 100. * athen / dnow / dnow;
	printf("g1, g2: %f %f %f %f\n", g1, g2, gdot1, gdot2);

	printf("Generating modes\n"); fflush(stdout);
	if (1) {
		gsl_rng_env_setup();
		T = gsl_rng_default;
		rng = gsl_rng_alloc(T);
		gsl_rng_set(rng,seed);
		ranc_ctr(seed);
		get_rhohat(dk, ng, bf1k);
	}
	else {
		FILE * modes_file = fopen("test_modes.dat", "rb");
		fread(bf1k, sizeof(double), 2 * ng3k, modes_file);
		fclose(modes_file);
		#pragma omp parallel for
		for (long i = 0; i < ng3k; i++)
			bf1k[i] /= sqrt_vol;
	}
	write_power(dk, ng, bf1k, g1, dk, dk, outprefix);

	printf("Computing 2LPT potential\n"); fflush(stdout);
	memcpy(bf3k, bf1k, ng3k * sizeof(fftwe_complex));
	memcpy(bf2k, bf1k, ng3k * sizeof(fftwe_complex));
	memcpy(bf0k, bf1k, ng3k * sizeof(fftwe_complex));
	greens(ng, dk, bf0k, 0, 0);
	greens(ng, dk, bf1k, 1, 1);
	greens(ng, dk, bf2k, 2, 2);
	fftwe_execute(bf0k_plan);
	fftwe_execute(bf1k_plan);
	fftwe_execute(bf2k_plan);
	get_l2abc(ng, bf2x, bf1x, bf0x);
	{
		int index1[3] = {0, 1, 2};
		int index2[3] = {1, 2, 0};
		for (int i = 0; i < 3; i++) {
			memcpy(bf1k, bf3k, ng3k * sizeof(fftwe_complex));
			greens(ng, dk, bf1k, index1[i], index2[i]);
			fftwe_execute(bf1k_plan);
			get_l2def(ng, bf2x, bf1x);
		}
	}
	fftwe_execute(bf2x_plan);

	{
		real norm_delta = -1. / sqrt_vol;
		real norm_l2 = 1. / (((real) ng3x) * boxsize * boxsize * boxsize);
		#pragma omp parallel for
		for (unsigned int ii = 0; ii < ng3k; ii++) {
			bf3k[ii] *= norm_delta;
			bf2k[ii] *= norm_l2;
		}
	}

	for (int i = 0; i < 3; i++) { 
		printf("Displacement %d \n", i); fflush(stdout);
			memcpy(bf0k, bf2k, ng3k * sizeof(fftwe_complex)); 
			memcpy(bf1k, bf3k, ng3k * sizeof(fftwe_complex)); 
			greens(ng, dk, bf0k, i, -1);
			greens(ng, dk, bf1k, i, -1);
			fftwe_execute(bf0k_plan);
			fftwe_execute(bf1k_plan);
			displace_2lpt(np, ng, boxsize, bf1x, g1, gdot1, bf0x, g2, gdot2, i); 
	}

	/*
	{
		real h_sqrta = hubble * sqrt(athen);
		long skip = np3 / 10;
		for (long i = skip / 2; i < np3; i += skip)
			printf("%06d pos: %+.4e %+.4e %+.4e | vel: %+.4e %+.4e %+.4e\n", i, xx[i], yy[i], zz[i], vx[i] / h_sqrta, vy[i] / h_sqrta, vz[i] / h_sqrta);
	}
	*/

	printf("Writing data ...\n"); fflush(stdout);
	omega_l = 0; 
	write_gadget(np,boxsize, zthen, hubble, omega_m, omega_l, outprefix, no_files);

	free(yy);
	free(zz);
	free(vx);
	free(vy);
	free(vz);
	fftwe_destroy_plan(bf0k_plan);
	fftwe_destroy_plan(bf1k_plan);
	fftwe_destroy_plan(bf2k_plan);
	fftwe_destroy_plan(bf2x_plan);
	fftw_cleanup_threads();
	fftwe_free(bf0k);
	fftwe_free(bf1k);
	fftwe_free(bf2k);
	fftwe_free(bf3k);
	free(xx);
	free(ktable);
	free(ptable);
	printf("done!\n");
	return 0;
}

void readIntParam(FILE * t_file, const char t_param_name[], int * t_value) {
    rewind(t_file);
	int value;
    char param_name[100];
	char line[MAXLINELEN];
    while (fgets(line, MAXLINELEN, t_file) != NULL) {
        sscanf(line,"%s", &param_name);
        if (!strcmp(param_name, t_param_name)) {
            sscanf(line,"%s %d", &param_name, &value);
            break;
        }
    }
    rewind(t_file);
    if (value == 0) {
        fprintf(stderr, "Error, no value found for %s in parameter file\n", t_param_name);
        exit(-1);
    }
	(*t_value) = value;
    return;
}

void readDoubleParam(FILE * t_file, const char t_param_name[], double * t_value) {
    rewind(t_file);
 	double value = 0;
    char param_name[100];
    char line[MAXLINELEN];
    while (fgets(line, MAXLINELEN, t_file) != NULL) {
        sscanf(line,"%s", &param_name);
        if (!strcmp(param_name, t_param_name)) {
            sscanf(line,"%s %lf", &param_name, &value);
            break;
        }
    }
    rewind(t_file);
    if (value == 0.) {
        fprintf(stderr, "Error, no value found for %s in parameter file\n", t_param_name);
        exit(-1);
    }
	(*t_value) = value;
    return;
}

void readStringParam(FILE * t_file, const char t_param_name[], char * t_value) {
    rewind(t_file);
    sprintf(t_value, "");
    char param_name[100];
    char line[MAXLINELEN];
    while (fgets(line, MAXLINELEN, t_file) != NULL) {
        sscanf(line,"%s", &param_name);
        if (!strcmp(param_name, t_param_name)) {
            sscanf(line,"%s %s", &param_name, &t_value[0]);
            break;
        }
    }
    rewind(t_file);
    if (!strcmp(t_value, "")) {
        fprintf(stderr, "Error, no value found for %s in parameter file\n", t_param_name);
        exit(-1);
    }
    return;
}
