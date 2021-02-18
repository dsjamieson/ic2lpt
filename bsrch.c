#include "real.h"

int bsrch(int narr,real *arr,real val)
// taken from Joe
/*  --- does a binary search for val in the array pointed
    by arr. array should be in ascending order */
{
  int low,high,mid;
  
  low=0;
  high=narr-1;
  
//printf("in bsrch searching for %f %d\n",val,narr);
  while((high-low)>1){
    mid=(low+high)/2;
    if(*(arr+mid)>val)high=mid;
    if(*(arr+mid)<=val)low=mid;
  }
  //  printf("end bsrch located %f at %d\n",val,high);
  return (high);
}


int bsrch_double(int narr, double *arr, double val)
// taken from Joe
/*  --- does a binary search for val in the array pointed
    by arr. array should be in ascending order */
{
  int low,high,mid;
  
  low=0;
  high=narr-1;
  
//printf("in bsrch searching for %f %d\n",val,narr);
  while((high-low)>1){
    mid=(low+high)/2;
    if(*(arr+mid)>val)high=mid;
    if(*(arr+mid)<=val)low=mid;
  }
  //  printf("end bsrch located %f at %d\n",val,high);
  return (high);
}


