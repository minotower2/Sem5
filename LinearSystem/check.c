#include "check.h"
#include <math.h>


double r1(double *a, double *a_rev, int n) {
  if (n <= 11000) return norm(a, a_rev, n);
  else return 0;
}

double r2(double *a, double *a_rev, int n){
  if (n <= 11000) return norm(a_rev, a, n);
  else return 0;
}

double norm(double *a, double *b, int n) {
  double sum = 0, res = 0, temp;
  int i, j, k, prod1, prod2;

  for (i = 0; i < n; i++) {
    prod1 = i*n;
    for (j = 0; j < i; j++) {
      prod2 = j*n;
      temp = b[prod1+j];
      b[prod1+j] = b[prod2+i];
      b[prod2+i] = temp;
    }
  }
  
  for (j = 0; j < n; j++) {
    prod2 = j*n;
    sum = 0;
    for (i = 0; i < n; i++) {
      prod1 = i*n;
      temp = 0;
      for (k = 0; k < n; k++) {
        temp += a[prod1 + k] * b[prod2 + k];
      }
      sum += (i != j ? fabs(temp) : fabs(temp-1));
    }
    //if (sum > res) res = sum;
    res = sum * (sum > res) + res * (res > sum);
  }
 
  for (i = 0; i < n; i++) {
    prod1 = i*n;
    for (j = 0; j < i; j++) {
      prod2 = j*n;
      temp = b[prod1+j];
      b[prod1+j] = b[prod2+i];
      b[prod2+i] = temp;
    }
  }
  return res;
}

