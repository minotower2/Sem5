#include <pthread.h>
#include "args_struct.h"

int equiv_double (double a, double b, double norm);
double norm_mat(double *a, int n);
int solve(double *a, double *x, double *a_rev, int n, double norm, pthread_t * threads, int total_threads, ARGS* args);
void inv(double *a, double *a_rev, int n, double norm);
void productHonest(double *x, double *a, int start, int end, int, int);
void productOptimized(double *x, double *y, int start, int end, int, int);
void * productHonestThreaded(void *pa);
void * productOptimizedThreaded(void *pa);
double inverse_sqrt_simple(double x, int iterations);
int gauss(double *a, int n, double *a_rev, double norm);
