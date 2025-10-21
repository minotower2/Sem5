#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "returns.h"
#include "solve.h"
#include "matrix.h"
#include "check.h"

int main(int argc, char ** argv) {
  double *a, *x, *a_rev;
  double t1, t2, res1 = 0, res2 = 0, norm;  
  int n, p, k, res = 0, result, task = 1;
  char *name = 0;

  if (!(
    (argc == 4 || argc == 5) && 
    sscanf(argv[1], "%d", &n) == 1 && n > 0 && 
    sscanf(argv[2], "%d", &p) == 1 && p >= 0 && 
    sscanf(argv[3], "%d", &k) == 1 && k >= 0 && k < 5)) 
  {
    printf("Usage %s size_of_matrix ammount_of_shown_elements formula (filename)\n", argv[0]);
    return 1;
  } 

  if (k == 0 && argc == 5) name = argv[4]; 
  else if ((argc == 5 && k > 0) || (argc == 4 && k == 0)) {
    printf("Usage %s size_of_matrix ammount_of_shown_elements formula (filename)\n", argv[0]); 
    return 1;
  }

  a = (double*)malloc(n*n*sizeof(double));
  if (!a) {printf("Not enough memory!\n"); return 2;} 
  x = (double*)malloc(n*sizeof(double));
  if (!x) {printf("Not enough memory!\n"); free(a); return 2;}
  a_rev = (double*)malloc(n*n*sizeof(double));
  if (!a_rev) {printf("Not enough memory!\n"); free(a); free(x); return 2;}

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j) a_rev[i*n+j] = 1;
      else a_rev[i*n+j] = 0;
    }
  }

  if (name) res = read_matrix_t(a, n, name);
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
    free(a_rev);
    return 3;
  }
  norm = norm_mat(a,n);
  printf("Initial matrix: \n");
  print_matrix_t(a, n,  p);
  t1 = clock();
  result = solve(a, x, a_rev, n, norm);
  t1 = (clock() - t1) / CLOCKS_PER_SEC;

  if (name) res = read_matrix(a, n, name);
  else init_matrix(a, n, k);

  if (result == SUCCESS) { 
    printf("Resulting matrix: \n");
    print_matrix(a_rev, n, p);
    t2 = clock();
    res1 = r1(a, a_rev, n);
    res2 = r2(a, a_rev, n);
    t2 = (clock() - t2) / CLOCKS_PER_SEC;
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d\n", argv[0], task, res1, res2, t1, t2, k, n);
  }
  else {
    printf("Division by 0\n");
  }
  free(a);
  free(x);
  free(a_rev);
  return 0;
}
