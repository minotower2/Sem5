#include "returns.h"
#include "matrix.h"
#include <stdio.h>

int read_matrix(double *a, int n, const char * name) {
  FILE *fp;
  int i, j;
  if (!(fp = fopen(name, "r"))) return ERROR_OPEN;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (fscanf(fp, "%lf", a + i*n+j) != 1) {
        fclose(fp);
        return ERROR_READ;
      }
    }
  }
  fclose(fp);
  return SUCCESS;
}
int read_matrix_t(double *a, int n, const char * name) {
  FILE *fp;
  int i, j;
  if (!(fp = fopen(name, "r"))) return ERROR_OPEN;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (fscanf(fp, "%lf", a + j*n+i) != 1) {
        fclose(fp);
        return ERROR_READ;
      }
    }
  }
  fclose(fp);
  return SUCCESS;
}

void print_matrix(double *a, int n, int p){
  int np = (n > p ? p : n);
  int i, j;
  for (i = 0; i < np; i++) {
    for (j = 0; j < np; j++) {
      printf(" %10.3e", a[i*n + j]);
    }
    printf("\n");
  }
}

void print_matrix_t(double *a, int n, int p){
  int np = (n > p ? p : n);
  int i, j;
  for (i = 0; i < np; i++) {
    for (j = 0; j < np; j++) {
      printf(" %10.3e", a[j*n + i]);
    }
    printf("\n");
  }
}
void init_matrix(double *a, int n, int k) {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      a[i*n + j] = f(k, n, i+1, j+1);
    }
  }
}

double f(int k, int n, int i, int j) {
  switch (k) {
    case 1: return n - max(i, j) +  1;
    case 2: return max(i, j);
    case 3: return (i - j >= 0 ? i - j : j - i);
    case 4: return 1./(i + j - 1);
  }
  return -1e300;
}

int max(int n, int m) {return (n > m ? n : m);}


