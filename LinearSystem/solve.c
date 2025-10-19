#include "solve.h"
#include "returns.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>
#define EPS 1e-16

int equiv_double (double a, double b, double norm){
  return (d_abs(b-a) < EPS * norm) ? 1 : 0;
}
double d_abs(double x) {return (x >= 0) ? x : -x;}

double norm_mat(double *a, int n) {
  double temp = 0, res = 0;
  int i, j;
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      temp += d_abs(a[i*n + j]);
    }
    if (temp > res) res = temp;
    temp = 0;
  }
  return res;
}

int solve(double *a, double *x, double *a_rev, int n, double norm){
  int k, j, i, sgn;
  double s, t, mod, modx;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i == j) a_rev[i*n+j] = 1;
      else a_rev[i*n+j] = 0;
    }
  }

  for (k = 0; k < n-1; k++) {
    modx = 0;
    s = 0;
    for (j = k+1; j < n; j++) {
      t = a[j*n + k];
      s += t*t;
    }
    t = a[k*n + k];
    sgn = t >= 0 ? 1 : -1;
    mod = sqrt(t*t + s);
    if (equiv_double(mod, 0, norm)) return DEV_BY_ZERO;
    modx = sgn*(a[k*n + k] - mod);
    x[k] = modx;
    modx = modx*modx+s;
    for (j = k+1; j < n; j++) {
      x[j] = sgn*a[j*n + k];
    }
    modx = sqrt(modx);
    if (equiv_double(modx, 0, norm)) continue;
    for (j = k; j < n; j++) {
      x[j] = x[j] / modx;
    }
    //printf("Vector x:\n");
    //for (j = k; j < n; j++) printf("%lf ", x[j]);
    //printf("\n");

    productMatrix(x, a, a, k, n);
    productVector(x, a_rev, a_rev, k, n);

    a[k*n+k] = mod;
    //printf("Matrix a:\n");
    //print_matrix(a, n, 5);
    //printf("Matrix a_rev:\n");
    //print_matrix(a_rev, n, 5);
  }

  for (i = n-1; i >= 0; i--) {
    s = a[i*n + i];
    if (equiv_double(s, 0, norm)) return DEV_BY_ZERO;
    for (j = 0; j < n; j++) {
      //a[i*n + j] /= s;
      a_rev[i*n + j] /= s;
    }
    for (j = i-1; j >= 0; j--) {
      s = a[j*n + i];
      if (!equiv_double(s, 0, norm)) {
        for (k = 0; k < n; k++) {
          //a[j*n + k] -= s * a[i*n + k];
          a_rev[j*n + k] -= s * a_rev[i*n + k];
        }
      }
    }
    //printf("Matrix a:\n");
    //print_matrix(a, n, 5);
    //printf("Matrix a_rev:\n");
    //print_matrix(a_rev, n, 5);
  }

  return SUCCESS;
}

void productVector(double *x, double *a, double *res, int start, int end) {
  int i, j;
  double scalarProduct;
  for (i = 0; i < end; i++) {
    scalarProduct = 0;
    for (j = start; j < end; j++) scalarProduct += x[j] * a[j*end+i];
    for (j = start; j < end; j++) res[j*end + i] = (a[j*end + i] - 2*scalarProduct*x[j]);
  }
}

void productMatrix(double *x, double *a, double *res, int start, int end) {
  int i, j;
  double scalarProduct;
  for (i = start+1; i < end; i++) {
    scalarProduct = 0;
    for (j = start; j < end; j++) scalarProduct += x[j] * a[j*end+i];
    for (j = start; j < end; j++) res[j*end + i] = (a[j*end + i] - 2*scalarProduct*x[j]);
  }
}
