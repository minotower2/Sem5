#include "check.h"
#include <math.h>
#include <stdio.h>


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
  int i, j, k;
  for (j = 0; j < (n < 10 ? n : n/10); j++) {
    sum = 0;
    for (i = 0; i < n; i++) {
      temp = 0;
      for (k = 0; k < n; k++) {
        temp += a[i*n + k] * b[k*n + j];
      }
      sum += (i == j ? fabs(temp-1) : fabs(temp));
    }
    if (sum > res) res = sum;
  }
  return res;
}

