int equiv_double (double a, double b, double norm);
double d_abs(double x);
double norm_mat(double *a, int n);
int solve(double *a, double *x, double *a_rev, int n, double norm);
int sluL(double *a, int k, double *x, int n);
int sluU(double *a, double *r, double *x, int n, double norm);
void inv(double *a, double *a_rev, int n, double norm);
void productHonest(double *x, double *a, int start, int end);
void productOptimized(double *x, double *y, int start, int end);
