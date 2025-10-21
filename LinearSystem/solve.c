#include "solve.h"
#include "returns.h"
#include "matrix.h"
#include <stdio.h>
#include <time.h>
#include <math.h>
#define EPS 1e-15

int equiv_double (double a, double b, double norm){
  return (fabs(b-a) < EPS * norm) ? 1 : 0;
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

int solve(double *a, double *x, double *a_rev, int n, double norm){
  int k, j, i;
  double s, t, mod, modx;
  double time;

  time = clock();
  for (k = 0; k < n-1; k++) {
    modx = 0;
    s = 0;
    for (j = k+1; j < n; j++) {
      t = a[k*n + j];
      s += t*t;
    }
    t = a[k*n + k];
    mod = sqrt(t*t + s);
    if (equiv_double(mod, 0, norm)) return DEV_BY_ZERO;
    modx = (t - mod);
    x[k] = modx;
    modx = modx*modx+s;
    if (equiv_double(modx, 0, norm)) continue;
    modx = 1. / sqrt(modx);
    x[k] *= modx;
    for (j = k+1; j < n; j++) {
      x[j] = a[k*n+ j] * (modx);
    }
    //printf("Vector x:\n");
    //for (j = k; j < n; j++) printf("%lf ", x[j]);
    //printf("\n");

    productOptimized(x, a, k, n);
    productHonest(x, a_rev, k, n);

    a[k*n+k] = mod;

    //printf("Matrix a:\n");
    //print_matrix(a, n, 5);
    //printf("Matrix a_rev:\n");
    //print_matrix(a_rev, n, 5);
  }

  time = clock() - time;
  time /= CLOCKS_PER_SEC;
  printf("Decomposition time: %lf\n", time);

  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      t = a_rev[i*n+j];
      a_rev[i*n+j] = a_rev[j*n+i];
      a_rev[j*n+i] = t;
    }
  }

  time = clock();
  for (i = n-1; i >= 0; i--) {
    s = a[i*n + i];
    if (equiv_double(s, 0, norm)) return DEV_BY_ZERO;
    s = 1. / s;
    for (j = 0; j < n; j++) {
      //a[i*n + j] /= s;
      a_rev[i*n + j] *= s;
    }
    for (j = i-1; j >= 0; j--) {
      s = a[i*n + j];
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
  time = clock() - time;
  time = time/CLOCKS_PER_SEC;
  printf("Inversion time: %lf\n", time);

  return SUCCESS;
}

void productHonest(double *x, double *a, int k, int n) {
  int i, j;
  double scalarProduct;
  for (i = 0; i < n; i++) {
    scalarProduct = 0;
    for (j = k; j < n; j++) scalarProduct += x[j] * a[i*n + j];
    scalarProduct *= 2;
    for (j = k; j < n; j++) a[i*n + j] -= scalarProduct * x[j];
  }
}

void productOptimized(double *x, double *a, int start, int end) {
  int i, j;
  double scalarProduct;
  for (i = start+1; i < end; i++) {
    scalarProduct = 0;
    for (j = start; j < end; j++) scalarProduct += x[j] * a[i*end+j];
    scalarProduct *= 2;
    for (j = start; j < end; j++) a[i*end + j] -= scalarProduct*x[j];
  }
}
