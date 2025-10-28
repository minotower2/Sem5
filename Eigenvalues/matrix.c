#include "returns.h"
#include "matrix.h"
#include <stdio.h>

int read_matrix(double *a, int n, int m, const char * name) {
  FILE *fp;
  int l, lm = n * m;
  if (!(fp = fopen(name, "r"))) return ERROR_OPEN;
  for (l = 0; l < lm; l++){
    if (fscanf(fp, "%lf", a + l) != 1) {
      fclose(fp);
      return ERROR_READ;
    }
  }
  fclose(fp);
  return SUCCESS;
}

void print_matrix(double *a, int n, int m, int p){
  int np = (n>p ? p : n);
  int mp = (m>p ? p : m);
  int i, j;
  for (i = 0; i < np; i++) {
    for (j = 0; j < mp; j++) {
      printf(" %10.3e", a[i*m + j]);
    }
    printf("\n");
  }
}

void init_matrix(double *a, int n, int k) {
  int i, j;
  for (i = 0; i< n; i++) {
    for (j = 0; j < n; j++) {
      a[i*n + j] = f(k, n, i+1, j+1);
    }
  }
}

double f(int k, int n, int i, int j) {
  switch (k) {
    case 1: return n - max(i, j) + 1;
    case 2: return 2 * (i == j) - 1 * (mod(i-j) == 1);
    case 3: return (i == j && j < n) + i * (j == n && j != i) + j * (i == n);
    case 4: return 1./(i + j - 1);
  }
  return -1e300;
}

int max(int n, int m) {return (n>m ? n : m);}

int mod(int n) {return n> 0 ? n : -n;}
