void productVector(double *y, int i, int j, double s, double c);
int equiv_double (double a, double b, double eps, double norm);
int equiv_double (double a, double b, double eps, double norm);
void productMatrix(double *a, int n, int i, int j, double s, double c);
void productMatrixHorizontal(double *a, int n, int i, int j, double s, double c);
double norm_mat(double *a, int n);
int solve(double *a, int n, double *x, double eps, double norm);
double sigma(double *a, int n);
