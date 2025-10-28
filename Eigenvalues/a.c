#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "returns.h"
#include "matrix.h"
#include "check.h"
#include "solve.h"

int main(int argc, char ** argv) {
  double *a, *x;
  double t1 = 0, t2=0, res1 = 0, res2 = 0, eps = 0, norm = 0;  
  int n, m, k, res = 0, result = 0, its = 0;
  char *name = 0;

  if (!(
    (argc == 5 || argc == 6) && 
    sscanf(argv[1], "%d", &n) == 1 && n > 0 && 
    sscanf(argv[2], "%d", &m) == 1 && m > 0 && 
    sscanf(argv[3], "%lf", &eps) == 1 && 
    sscanf(argv[4], "%d", &k) == 1 && k >= 0 && k < 5)) 
  {
    printf("Usage %s size printed eps k (filename)\n", argv[0]);
    return 1;
  } 

  if (k == 0 && argc == 6) name = argv[5]; 
  else if ((argc == 6 && k > 0) || (argc == 5 && k == 0)) {
    printf("Usage %s size printed eps k (filename)\n", argv[0]); 
    return 1;
  }

  a = (double*)malloc(n*n*sizeof(double));
  if (!a) {printf("Not enough memory!\n"); return 2;} 
  x = (double*)malloc(n*sizeof(double));
  if (!x) {printf("Not enough memory!\n"); free(a); return 2;} 


  if (name) res = read_matrix(a, n, n, name);
  else init_matrix(a, n, k);
  if (res!= SUCCESS) {
    switch (res) {
      case ERROR_OPEN:
        printf("Couldn't open %s\n", name);
        break;
      case ERROR_READ: 
        printf("Couldn't read %s\n", name);
        break;
      default:
        printf("Unknown error %d in %s, %s\n", res, name ? "file" : "init", name ? name : " ");
        break;
     }      
    free(a);
    free(x);
    return 3;
  }
  printf("Initial matrix: \n");
  print_matrix(a, n, n, m);
  for (int i = 0; i < n; i++) x[i] = 0;

  norm = norm_mat(a, n);

  t1 = clock();
  its = solve(a, n, x, eps, norm);
  t1 = (clock() - t1)/CLOCKS_PER_SEC;

  if (name) res = read_matrix(a, n, n, name);
  else init_matrix(a, n, k);

  t2 = clock();
  res1 = r1(a, x, n, norm);
  res2 = r2(a, x, n, norm);
  t2 = (clock() - t2)/CLOCKS_PER_SEC;
  
  if (result == 0) { 
    printf("Eigenvalues: \n");
    print_matrix(x, 1, n, m);
    printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n", argv[0], res1, res2, its, its / n, t1, t2);  }
  else {
    printf("Division by 0\n");
  }
  free(a);
  free(x);
  return 0;
}
