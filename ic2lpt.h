#include "fftwe.h"
#include "real.h"

#define C_LIGHT 2997.92458
#define pi 3.1415926535897932384626

extern real dk;
extern real omega_m;
extern real omega_l;
extern real hubble;

void write_power (real dk, int NGRID, fftwe_complex* complex_delta, real g1, real BIN, real START, char *fname);
int enforce_hermitian(int n,fftwe_complex *array);
real power_interp(real kk);
double ranc_ctr(int iseed);
int get_rhohat(real dk,int n,fftwe_complex *rhohat);
int greens(int n,real dk, fftwe_complex *rhohat, int d1, int d2);
int get_l2abc(int ng, real *a, real *b, real *c);
int get_l2def(int ng, real *a, real *d);
int displace_2lpt(int np, int ng, real boxsize, real *psi1, real g1, real gdot1,real *psi2, real g2, real gdot2, int direction);
void write_gadget(int np,double b,double z_i, double h,double omega_m, double omega_lambda, char*outprefix,int no_files);
real power_interp(real kk);
int bsrch(int narr,real *arr,real val);
