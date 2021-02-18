#include "ic2lpt.h"
#include <math.h>

int bsrch(int narr, real *arr, real val);

real power_interp(real kk) {
  int i0, i1;
  extern int nk;
  extern real * ktable, * ptable;
  real slope,power; 
  if (kk == 0) {
    if (ktable[0] == 0.) 
      power = ptable[0]; 
    else 
		power = 0.;
  }
	else {
		i1= bsrch(nk, ktable, kk);
		i0 = i1-1;
		slope = (ptable[i1] - ptable[i0]) / (ktable[i1] - ktable[i0]);
		power = ptable[i0] + slope * (kk - ktable[i0]);
	}
	return power;
}
