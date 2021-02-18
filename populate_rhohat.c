#include <stdlib.h>
#include "real.h"
#include <complex.h>
#include "fftwe.h"
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define VERBOSE 0

#define pi 3.1415926535897932384626

int glass = 0;

extern const gsl_rng_type * T;
extern gsl_rng * rng;

//double gauss(){
//  double u1,u2,v;

// v=1;
// while (v>=1){
//     u1=drand48();
//     u2=drand48();
//     v=pow(2*u1-1,2)+pow(2*u2-1,2);
// }

// return (2*u1-1)*sqrt(-2*log(v)/v);

//}


//double cgasdev_ctr();
//double ranc_ctr();
real power_interp(real kk);


real test_var(real g1, real rk, real bin, int nmodes, int N){
  int i,j;
  real kk;
  real Pk[10000];
  real Pk_mean=0;
  real Pk_var=0;
  real gauss1,gauss2;
  real sigma;
  char string[32];
     
  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();
  
  T = gsl_rng_default;//gsl_rng_ranlxd1;
  r = gsl_rng_alloc (T);
     
  for (j=0;j<N;j++){
    Pk[j]=0;
    //        ranc_ctr(j);
    //    srand48(j);
    gsl_rng_set(r,j);
    for (i=0;i<nmodes;i++){
      //      kk=/*gsl_rng_uniform(r)drand48()*/ranc_ctr(0)*bin-0.5*bin+rk;
      kk=rk;
      gauss1 = gsl_ran_gaussian(r,1.);//gauss();//gsl_ran_gaussian(r,1.);
      gauss2 = gsl_ran_gaussian(r,1.);//gauss();//cgasdev_ctr();
      sigma=g1*sqrt(power_interp(kk)/2.);
      Pk[j]+=(pow(gauss1*sigma,2)+pow(gauss2*sigma,2));
    }
    Pk[j]/=nmodes;

    //   only look at the first 6 digits
    //    sprintf(string,"%g",Pk[j]);
    //    Pk[j]=atof(string);
    
    

    //    printf("%d: %g %s\n",j,Pk[j],string);
  }




  for (j=0;j<N;j++){
    Pk_mean+=Pk[j];
  }
  Pk_mean/=N;

  printf("Pk_mean %g\n",Pk_mean);

  for (j=0;j<N;j++){
    Pk_var+=((Pk[j]-Pk_mean)*(Pk[j]-Pk_mean));
  }

  Pk_var/=(N-1.);

  printf("Pk_var %g\n",Pk_var);

  gsl_rng_free (r);
  
  return Pk_var;

}


int populate_rhohat(int n,int k,int j,int i,real dk,fftwe_complex *delta,int seed) {
/* k: level, j: row, i: column.
   This routine will populate all 4 points in delta that correspond
   to the positive frequency (i,j,k).
   It would be 8, but remember half of the array is redundant.  Thus,
   there is no such considering negative "i."
   Complicated feature: see those if (j)'s? They used to be 
   if (j && j!=nby2) because you wouldn't want to overwrite the
   Nyquist frequencies you already calculated, right?  Turns out,
   in order to be consistent with other size grids, you should overwrite 
   them, because then the random seeds get called the right number of times.
   Not exactly true: when j==n/2, populate power twice to account for
   both positive and negative frequencies. Hence the +='s.
   dk is 2pi/boxsize.  */
/* 2004-10-05.  Should you += or just =?  Remember you're taking an infinite
   series and truncating at a Nyquist frequency.  This corresponds to 
   convolving the correlation function with a sinc.  So do you truncate
   in Fourier space with a box of length n, so that you only get one Nyquist
   frequency, or do you truncate with a box of length n+1, so that you alias
   the two contributions from opposite ends of the power spectrum? It seems
   more proper to truncate with box of size n.  Hence, just ='s.
*/
  long np1,np1xn,indpp,indpm,indmp,indmm;
  real kk,power,expectation,gauss1,gauss2;
  real kx,ky,kz,k_ny,wk,wi,wj,w,argu;



#if VERBOSE
  static int counter=0;
  counter++;
  printf("\n");
#endif

  np1 = n/2+1;
  np1xn = np1 * n;
  kk = dk * sqrt((real)(k*k + j*j + i*i));
  kx = dk * i;
  ky = dk * j;
  kz = dk * k;
  k_ny=n/2*dk;

  // find primary index in delta array (the positive frequency)
  indpp = k*np1xn + j*np1 + i;
  // next index: same k, negative "j", same i
  indpm = k*np1xn + (n-j)*np1 + i;
  // next index: negative k, same j, same i
  indmp = (n-k)*np1xn + j*np1 + i;
  // next index: negative k, negative j, same i
  indmm = (n-k)*np1xn + (n-j)*np1 + i;
#if VERBOSE
//  printf("indices %d %d %d %d %d\n",indpp,indpm,indmp,indmm);
#endif


  power = power_interp(kk);
  
  if (glass){
    argu=(pi*kx)/(2.*k_ny);
    if(i)
      wi=sin(argu)/argu;
    else
      wi=1;
    argu=(pi*ky)/(2.*k_ny);
    if(j)
      wj=sin(argu)/argu;
    else
      wj=1;
    argu=(pi*kz)/(2.*k_ny);
    if(k)
      wk=sin(argu)/argu;
    else
      wk=1;

    w=wi*wk*wj;
    w*=w;
    w*=w;

    power/=w;
  }

  expectation = sqrt(power/2.);
/* My discrete fourier convention is to include a factor of L^3 so that
   the discrete case corresponds exactly to the continuous case in the 
   limit N -> infinity, L -> infinity, L/N -> 0.  Therefore, there should
   be a factor of L^3/2 multiplied here (sqrt because the power spectrum
   is the fourier transform of the correlation function).  However, for
   computational efficiency, DO NOT do that here, but wait till the calling
   routine inverse fourier transforms delta hat, which introduces a factor
   of 1/L^3.  The net result is that in ic, I divide by L^3/2. These
   factors of L are especially important when mixing continuous and 
   discrete formulae, which I do to compute sigma8. */


  // always fill the positive frequency
  //  gauss1 = cgasdev_ctr();
  //  gauss2 = cgasdev_ctr();

  gauss1 = gsl_ran_gaussian(rng,1.);
  gauss2 = gsl_ran_gaussian(rng,1.);
  delta[indpp] = expectation*gauss1 + I*expectation*gauss2;
#if VERBOSE
  printf("pp %3d %3d %3d %7.2f %7.2f %7.2f\t%4d\n",
  		  k,j,i,expectation,gauss1,gauss2,counter);
#endif

  if (j) { // don't fill frequencies where n-j = n
    //    gauss1 = cgasdev_ctr();
    //    gauss2 = cgasdev_ctr();

    gauss1 = gsl_ran_gaussian(rng,1.);
    gauss2 = gsl_ran_gaussian(rng,1.);
    delta[indpm] = expectation*gauss1 + I*expectation*gauss2;
#if VERBOSE
    printf("pm %3d %3d %3d %7.2f %7.2f %7.2f\t%4d\n",
	   k,n-j,i,expectation,gauss1,gauss2,counter);
#endif
  }

  if (k) {
    //    gauss1 = cgasdev_ctr();
    //    gauss2 = cgasdev_ctr();
    gauss1 = gsl_ran_gaussian(rng,1.);
    gauss2 = gsl_ran_gaussian(rng,1.);
    delta[indmp] = expectation*gauss1 + I*expectation*gauss2;
#if VERBOSE
    printf("mp %3d %3d %3d %7.2f %7.2f %7.2f\t%4d\n",
	   n-k,j,i,expectation,gauss1,gauss2,counter);
#endif
  }
  
  if (j && k) {
    //    gauss1 = cgasdev_ctr();
    //    gauss2 = cgasdev_ctr();
    gauss1 = gsl_ran_gaussian(rng,1.);
    gauss2 = gsl_ran_gaussian(rng,1.);
    delta[indmm] = expectation*gauss1 + I*expectation*gauss2;
#if VERBOSE
    printf("mm %3d %3d %3d %7.2f %7.2f %7.2f\t%4d\n",
	   n-k,n-j,i,expectation,gauss1,gauss2,counter);
#endif
  }


  return 0;
}
