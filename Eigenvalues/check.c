#include "check.h"
#include <math.h>


double r1(double *a, double *x, int n, double norm) {
  double temp = 0;
  for (int i = 0; i < n; i++) temp+=(a[i*n + i] - x[i]);
  return fabs(temp)/norm;
}

double r2(double *a, double *x, int n, double norm) {
  int i, j;
  double temp1 = 0, temp2 = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      temp1 += a[i*n + j] * a[j*n + i];      
    }
  }
  for (i = 0; i < n; i++) temp2 += x[i] * x[i];
  return fabs(sqrt(temp1) - sqrt(temp2)) / norm;
}

