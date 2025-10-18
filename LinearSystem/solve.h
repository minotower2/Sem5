int equiv_double (double a, double b, double norm);
double d_abs(double x);
double norm_mat(double *a, int n);
int solve(double *a, double *x, double *a_rev, int n, double norm);
int sluL(double *a, int k, double *x, int n);
int sluU(double *a, double *r, double *x, int n, double norm);
void inv(double *a, double *a_rev, int n, double norm);
void productMatrix(double *x, double *a, double *res, int start, int end);
void productVector(double *x, double *y, double *res, int start, int end);
