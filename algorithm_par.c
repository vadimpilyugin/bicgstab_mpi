#include "math.h"
#include "stdio.h"
#include "assert.h"
#include "time.h"

#include "core.h"
#include "algorithm_par.h"
#include "mpi.h"

extern MPI_Comm comm_grid;
extern int i_am_the_master;

/* Parallel operations */

static double max(double a1, double a2) {
  return a1 > a2 ? a1 : a2;
}

void SpMV_par(CSRMatrix A, Vector X, Vector Y) {
  assert(A.N == X.ni);
  assert(A.N == Y.ni);
  assert(X.nh == Y.nh);

  #pragma omp parallel for
  for (int i = 0; i < A.N; i++) {
    Y.data[i] = 0.0;
    int j;

    for (j = A.IA[i]; j < A.IA[i+1]; j++) {
      Y.data[i] += X.data[A.JA[j]] * A.A[j];
    }
  }

  do_exchange(Y);
}

void axpby_par(Vector X, Vector Y, double a, double b) {
  assert(X.size == Y.size);
  assert(X.ni == Y.ni);
  assert(X.nh == Y.nh);

  #pragma omp parallel for
  for (int i = 0; i < (X.ni + X.nh); i++) {
    X.data[i] = a*X.data[i] + b*Y.data[i];
  }
  // if X and Y are halo-consistent then we can skip exchange
  // and just calculate those values locally
  // X will stay halo-consistent
}

double dot_par(Vector X, Vector Y) {
  assert(X.size == Y.size);
  assert(X.ni == Y.ni);
  assert(X.nh == Y.nh);

  double s = 0;
  #pragma omp parallel for reduction(+:s)
  for (int i = 0; i < X.ni; i++) {
    s += X.data[i] * Y.data[i];
  }

  double s_all = 0;
  MPI_Allreduce(&s, &s_all, 1, MPI_DOUBLE, MPI_SUM, comm_grid);
  return s_all;
}
void assign_vector_par(Vector X, Vector Y) {
  assert(X.size == Y.size);
  assert(X.ni == Y.ni);
  assert(X.nh == Y.nh);

  #pragma omp parallel for
  for (int i = 0; i < (X.ni + X.nh); i++)
    X.data[i] = Y.data[i];
}
void fill_with_constant_par(Vector X, const int c) {
  #pragma omp parallel for
  for (int i = 0; i < (X.ni + X.nh); i++) {
    X.data[i] = c;
  }
}

int bicgstab_solve_par(const CSRMatrix A, const Vector BB, Vector XX, int N, 
  double tol, int maxit, int info) {

  double Rhoi_1 = 1.0;
  double alphai = 1.0;
  double wi = 1.0;
  double betai_1 = 1.0;
  double Rhoi_2 = 1.0;
  double alphai_1 = 1.0;
  double wi_1 = 1.0;
  double RhoMin = 1E-60;
  double mineps = 1E-15;

  /*
    Before cycle:
    axpby - 0
    spmv - 0
    dot - 1

    First iteration:
    axpby - 4
    spmv - 4
    dot - 4

    Other iterations:
    axpby - 6
    spmv - 4
    dot - 5

    TOTAL_GFLOP = 1*DOT_GFLOP + 4*(AXPBY_GFLOP + SPMV_GFLOP + DOT_GFLOP) + 
      i * (6*AXPBY_GFLOP + 4*SPMV_GFLOP + 5*DOT_GFLOP)
  */

  CSRMatrix DD = csr_diag_inverse(A);
  Vector RR   = make_vector(N);
  Vector RR2  = make_vector(N);
  Vector PP   = make_vector(N);
  Vector PP2  = make_vector(N);
  Vector VV   = make_vector(N);
  Vector SS   = make_vector(N);
  Vector SS2  = make_vector(N);
  Vector TT   = make_vector(N);

  fill_with_constant_par(XX, 0);
  assign_vector_par(RR, BB);
  assign_vector_par(RR2, BB);
  assign_vector_par(PP, RR);

  double initres = sqrt(dot_par(RR,RR));
  double eps = max(mineps, tol*initres);
  double res = initres;
  int i;

  for (i = 0; i < maxit; i++) {
    if (info && i_am_the_master) {
      printf("It %d: res = %e tol=%e\n", i, res, res / initres); 
    }
    if (res < eps)
      break; 
    if (res > initres / mineps)
      return -1; 
    if (i == 0)
      Rhoi_1 = initres * initres;
    else
      Rhoi_1 = dot_par(RR2, RR);
    if (fabs(Rhoi_1) < RhoMin)
      return -1;
    if (i > 0) {
      betai_1 = (Rhoi_1 * alphai_1) / (Rhoi_2 * wi_1); 
      axpby_par(PP, RR, betai_1, 1.0);
      axpby_par(PP, VV, 1.0, -wi_1 * betai_1); 
    }
    SpMV_par(DD, PP, PP2); 
    SpMV_par(A, PP2, VV);  
    alphai = dot_par(RR2, VV); 
    if (fabs(alphai) < RhoMin)
      return -3;
    alphai = Rhoi_1 / alphai;
    assign_vector_par(SS, RR);
    axpby_par(SS, VV, 1.0, -alphai);
    SpMV_par(DD, SS, SS2); 
    SpMV_par(A, SS2, TT);
    wi = dot_par(TT, TT);
    if (fabs(wi) < RhoMin)
      return -4; 
    wi = dot_par(TT, SS) / wi;  
    if (fabs(wi) < RhoMin)
      return -5;
    axpby_par(XX, PP2, 1.0, alphai); 
    axpby_par(XX, SS2, 1.0, wi);
    assign_vector_par(RR, SS);
    axpby_par(RR, TT, 1.0, -wi); 
    alphai_1 = alphai;  
    Rhoi_2 = Rhoi_1;  
    wi_1 = wi; 
    res = sqrt(dot_par(RR,RR)); 
  }
  if (info && i_am_the_master)
    printf("Solver_BiCGSTAB_par: outres: %g\n", res);
  free_vector(&RR);
  free_vector(&RR2);
  free_vector(&PP);
  free_vector(&PP2);
  free_vector(&VV);
  free_vector(&SS);
  free_vector(&SS2);
  free_vector(&TT);
  csr_free(&DD);
  return i; 
}