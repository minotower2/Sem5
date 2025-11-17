#include "solve.h"
#include "returns.h"
#include "matrix.h"
#include <math.h>
#include "synchronize.h"
#include "args_struct.h"
#include <stdio.h>
#include <time.h>
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

int solve(double *a, double *x, double *a_rev, int n, double norm, pthread_t * threads, int total_threads, ARGS* args){
  int k, j, i, prod1, prod2;
  double s, t, mod, modx;
  for (i = 0; i < total_threads; i++) {
    args[i].x = x;
    args[i].n = n;
    args[i].total_threads = total_threads;
    args[i].thread_num = i;
  }

  for (k = 0; k < n-1; k++) {
    for (i = 0; i < total_threads; i++) {
      args[i].k = k;
    }
    prod1 = k*n;
    modx = 0;
    s = 0;
    for (j = k+1; j < n; j++) {
      t = a[prod1 + j];
      s += t*t;
    }
    t = a[prod1 + k];
    mod = sqrt(t*t + s);
    if (equiv_double(mod, 0, norm)) return DEV_BY_ZERO;
    modx = (t - mod);
    x[k] = modx;
    modx = modx*modx+s;
    if (equiv_double(modx, 0, norm)) continue;
    modx = 1. / sqrt(modx);
    x[k] *= modx;
    for (j = k+1; j < n; j++) {
      x[j] = a[prod1 + j] * (modx);
    }
    //printf("Vector x:\n");
    //for (j = k; j < n; j++) printf("%lf ", x[j]);
    //printf("\n");

    for (i = 0; i < total_threads; i++) {
      args[i].a = a_rev;
      if (pthread_create(threads+i, 0, productHonestThreaded, args+i)) {
        printf("Can't create thread %d\n", i);
        return -1;
      }
    }
    for (i = 0; i < total_threads; i++) {
      if (pthread_join (threads[i], 0)){
        printf("Can't wait thread %d\n", i);
        return -1;
      }
    }

    for (i = 0; i < total_threads; i++) {
      args[i].a = a;
      if (pthread_create(threads+i, 0, productOptimizedThreaded, args+i)) {
        printf("Can't create thread %d\n", i);
        return -1;
      }
    }
    for (i = 0; i < total_threads; i++) {
      if (pthread_join (threads[i], 0)){
        printf("Can't wait thread %d\n", i);
        return -1;
      }
    }

    //for (i = 0; i < total_threads; i++) productOptimized(x, a, k, n, i, total_threads);
    //for (i = 0; i < total_threads; i++) productHonest(x, a_rev, k, n, i, total_threads);

    a[prod1 + k] = mod;

    //printf("Matrix a:\n");
    //print_matrix(a, n, 5);
    //printf("Matrix a_rev:\n");
    //print_matrix(a_rev, n, 5);
  }

  for (i = 0; i < n; i++) {
    prod1 = i*n;
    for (j = 0; j < i; j++) {
      prod2 = j*n;
      t = a_rev[prod1 + j];
      a_rev[prod1 + j] = a_rev[prod2 + i];
      a_rev[prod2 + i] = t;
    }
  }

  if (gauss(a, n, a_rev, norm)) return DEV_BY_ZERO;
    //printf("Matrix a:\n");
    //print_matrix(a, n, 5);
    //printf("Matrix a_rev:\n");
    //print_matrix(a_rev, n, 5);
  return SUCCESS;
}

void productHonest(double *x, double *a, int k, int n, int thread_num, int total_threads) {
  int i, j, prod;
  double scalarProduct;
  int first_row = n*thread_num;
  int last_row = n*thread_num + n;
  first_row /= total_threads;
  last_row = last_row / total_threads;
  for (i = first_row; i < last_row; i++) {
    prod = i*n;
    scalarProduct = 0;
    for (j = k; j < n; j++) scalarProduct += x[j] * a[prod + j];
    scalarProduct *= 2;
    for (j = k; j < n; j++) a[prod + j] -= scalarProduct * x[j];
  }
  synchronize(total_threads);
}

void productOptimized(double *x, double *a, int k, int n, int thread_num, int total_threads) {
  int i, j, prod;
  double scalarProduct;
  int c = (n-k) / total_threads;
  int first_row = k+1 + thread_num * c;
  int last_row = thread_num+1 != total_threads ?  k+1 + thread_num * c + c : n; 
  for (i = first_row; i < last_row; i++) {
    prod = i*n;
    scalarProduct = 0;
    for (j = k; j < n; j++) scalarProduct += x[j] * a[prod + j];
    scalarProduct *= 2;
    for (j = k; j < n; j++) a[prod + j] -= scalarProduct * x[j];
  }
  synchronize(total_threads);
}



void * productHonestThreaded(void *pa)
{
  ARGS *pargs = (ARGS*)pa;
  //printf ("Thread %d started\n", pargs->thread_num);
  productHonest(pargs->x, pargs->a, pargs->k, pargs->n, pargs->thread_num, pargs->total_threads);
  //printf("Thread %d finished, time = %f\n", pargs->thread_num, t/CLOCKS_PER_SEC);
  return 0;
}

void * productOptimizedThreaded(void *pa)
{
  ARGS *pargs = (ARGS*)pa;
  //printf ("Thread %d started\n", pargs->thread_num);
  productOptimized(pargs->x, pargs->a, pargs->k, pargs->n, pargs->thread_num, pargs->total_threads);
  //printf("Thread %d finished, time = %f\n", pargs->thread_num, t/CLOCKS_PER_SEC);
  return 0;
}


int gauss(double *a, int n, double *a_rev, double norm) {
  int i, prod1, j, prod2, k;
  double s;
  
  for (i = n-1; i >= 0; i--) {
    prod1 = i*n;
    s = a[prod1 + i];
    if (equiv_double(s, 0, norm)) return DEV_BY_ZERO;
    s = 1. / s;
    for (j = 0; j < n; j++) {
      //a[i*n + j] /= s;
      a_rev[prod1 + j] *= s;
    }
    for (j = i-1; j >= 0; j--) {
      prod2 = j*n;
      s = a[prod1 + j];
      if (!equiv_double(s, 0, norm)) {
        for (k = 0; k < n; k++) {
          //a[j*n + k] -= s * a[i*n + k];
          a_rev[prod2 + k] -= s * a_rev[prod1 + k];
        }
      }
    }
  }
  return 0;
}
