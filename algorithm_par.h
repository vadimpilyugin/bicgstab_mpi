#include "create_csr.h"

void SpMV_par(CSRMatrix A, Vector X, Vector Y);
void axpby_par(Vector X, Vector Y, double a, double b);
double dot_par(Vector X, Vector Y);
void assign_vector_par(Vector X, Vector Y);
void fill_with_constant_par(Vector X, const int c);
int bicgstab_solve_par(CSRMatrix A, Vector BB, Vector XX, int N, double tol, int maxit, int info);