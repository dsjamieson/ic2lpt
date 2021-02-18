#include "ic2lpt.h"
#include <math.h>

int bsrch(int narr,real *arr,real val);


real power_interp(real kk) {
  int i0,i1;
  extern int nk;
  extern real *ktable, *ptable;
  real slope,power;

      
  // linearly interpolate in log space
  if (kk == 0) { // handle k=0 differently
    if (ktable[0] == 0.) 
      power = ptable[0];//pow(10.,ptable[0]); 
    else power = 0.;
  } else {
    i1= bsrch(nk,ktable,kk);
    i0 = i1-1;
    //    if ((ktable[0] == 0.) && (i0 == 0)) { // handle kk in first bin when k[0]=0.
      slope = (ptable[i1] - ptable[i0]) / (ktable[i1] - ktable[i0]);
      power = ptable[i0] + slope * (kk - ktable[i0]);
      //    } else {
      //      slope = (ptable[i1] - ptable[i0]) / log10(ktable[i1]/ktable[i0]);
      //      power = ptable[i0] + slope * log10(kk/ktable[i0]);
      //    }
      //      power = exp(power);
  }



  return power;
}
