#include "real.h"

int bsrch(int narr, real * arr, real val) {
	int low, high, mid;
  	low = 0;
	high = narr - 1;
	while ( high - low > 1) {
    	mid = (low + high) / 2;
    	if (*(arr + mid) > val)
			high = mid;
    	if (*(arr + mid) <= val)
			low = mid;
  }
  return high;
}
