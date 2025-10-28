#include "solve.h"
#include <math.h>
#include <stdio.h>

int equiv_double (double a, double b, double eps, double norm){
  return (fabs(b-a) < eps* norm) ? 1 : 0;
}

double norm_mat(double *a, int n) {
  double temp = 0, res = 0;
  int i, j;
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      temp += fabs(a[i*n + j]);
    }
    if (temp > res) res = temp;
    temp = 0;
  }
  return res;
}

int solve(double *a, int n, double *x, double eps, double norm) {
  int its = 0, i, j, sign, prod;
  double sigm = sigma(a, n), xs, ys, c, s, temp, aij;
  while (sigm > eps) {
    for (i = 0; i < n; i++) {
      prod = i*n;
      for (j = i+1; j < n; j++) {
        if(sigm < eps) break;
        aij = a[prod+j];
        its++;
        xs = -2 * aij;
        ys = a[prod+i] - a[j*n+j];
        if (equiv_double(ys, 0, eps, norm)) {
          c = 1/sqrt(2);
          s = 1/sqrt(2);
        }
        else {
          temp = (sqrt(xs*xs + ys*ys));
          c = sqrt(0.5*(1 + (fabs(ys))/temp));
          //sign = (xs*ys > 0 ? 1 : -1);
          sign = (xs*ys > 0) - (xs*ys <= 0);
          s = (sign * fabs(xs)) / (2*c * temp);
        }
        productMatrix(a, n, i, j, s, c);
        productMatrixHorizontal(a, n, i, j, s, c);
        sigm -= 2* aij*aij;
      }
    }
  }
  for (i = 0; i < n; i++) x[i] = a[i*n + i];
  return its;
}

void productVector(double *y, int i, int j, double s, double c) {
  double temp1, temp2;
  temp1 = y[i];
  temp2 = y[j];
  y[i] = temp1*c - temp2*s;
  y[j] = temp1*s + temp2*c;
}


void productMatrix(double *a, int n, int i, int j, double s, double c) {
  double temp1, temp2;
  int prod1 = i*n, prod2 = j*n;
  for (int k = 0; k < n; k++) {
    temp1 = a[prod1 + k];
    temp2 = a[prod2 + k];
    a[prod1 + k] = temp1*c - temp2*s;
    a[prod2 + k] = temp1*s + temp2*c;
  }
}

void productMatrixHorizontal(double *a, int n, int i, int j, double s, double c) {
  double temp1, temp2;
  for (int k = 0; k < n; k++) {
    temp1 = a[k*n+i];
    temp2 = a[k*n+j];
    a[k*n+i] = c*temp1 - s*temp2;
    a[k*n+j] = s*temp1 + c*temp2;
  }
}
double sigma(double *a, int n) {
  int i, j;
  double temp = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i != j) temp += a[i*n+j] * a[i*n+j];
    }
  }
  return temp;
}

