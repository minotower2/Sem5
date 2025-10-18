#include "check.h"

double dabs(double x) {return x > 0 ? x : -x;}

double r1(double *a, double *a_rev, int n) {
  if (n <= 11000) return norm(a, a_rev, n, n);
  else return 0;
}

double r2(double *a, double *a_rev, int n){
  if (n <= 11000) return norm(a_rev, a, n, n);
  else return 0;
}

double norm(double *a, double *b, int n, int m) { 
  double sum = 0, res = 0, temp;
  int i, j, k;
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      temp = 0;
      for (k = 0; k < m; k++) {
        temp += a[i*m + k] * b[k*n + j];
      }
      sum += (i == j ? dabs(temp-1) : dabs(temp));
    }
    if (sum > res) res = sum;
    sum = 0;
  }
  return res;
}

