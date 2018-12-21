#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "algorithm_par.h"
#include "division.h"
#include "create_csr.h"
#include "core.h"
#include "mpi.h"
#include "omp.h"

#define FAIL_NPROC 2

int myrank, nproc, i_am_the_master;
MPI_Comm comm_grid;
int old_nproc;

const char *WRONG_NPROC = "Number of processes is not a power of 2";

void my_assert (int expr, const char *s) {
  if (!expr) {
    fprintf(stderr, "%s\n", s);
    exit(1);
  }
}

void test_create_csr() {

  int max_values[N_DIMENSIONS] = {2,2,2};

  old_nproc = nproc;
  int P[N_DIMENSIONS];
  int grid_nproc = make_grid(max_values, P);

  if (myrank < grid_nproc) {
    nproc = grid_nproc;

    int n_local;
    int *Rows = NULL;
    int p[N_DIMENSIONS];
    int ps[N_DIMENSIONS];
    int n_start[N_DIMENSIONS];
    grid_params(myrank, max_values, P, p, ps, n_start, &n_local);
    CSRMatrix m = csr_matrix(max_values, n_start, ps, n_local, &Rows);
    my_assert(Rows != NULL, "Rows are NULL!");

    int IA[] = {0, 4, 8, 12, 16, 20, 24, 28, 32};
    int JA[] = {0, 1, 2, 4, 0, 1, 3, 5, 0, 2, 3, 6, 1, 2, 3, 7, 0, 4, 5, 6, 1, 4, 5, 7, 2, 4, 6, 7, 3, 5, 6, 7};
    double A[] = {2.2102758805035565, 0.9092974268256817, 0.1411200080598672, -0.9589242746631385, 0.9092974268256817, 2.77772913022837, -0.9589242746631385, 0.6569865987187891, 0.1411200080598672, 0.9159193906506048, -0.27941549819892586, 0.4121184852417566, -0.9589242746631385, -0.27941549819892586, 2.462162977354045, -0.9999902065507035, -0.9589242746631385, 2.753229151313533, -0.5440211108893698, -0.9999902065507035, 0.6569865987187891, -0.5440211108893698, 1.7832922210782798, 0.4201670368266409, 0.4121184852417566, -0.9999902065507035, 2.6429876522360636, 0.9906073556948704, -0.9999902065507035, 0.4201670368266409, 0.9906073556948704, 2.6518410589794366};

    if (i_am_the_master)
      print_csr(m);

    double eps = 1e-7;

    for (int i = 0; i < m.N; i++) {
      int glob_row = Rows[i];
      int j1 = m.IA[i];
      int j2 = IA[glob_row];
      int len1 = m.IA[i+1] - m.IA[i];
      int len2 = IA[glob_row+1] - IA[glob_row];
      my_assert(len1 == len2, "Lengths of rows are not equal!");
      for (int k = 0; k < len1; k++) {
        my_assert(m.JA[j1] == JA[j2], "JA are not equal!");
        my_assert(fabs(m.A[j1] - A[j2]) < eps, "A's are not equal!");
        j1++;
        j2++;
      }
    }
    csr_free(&m);
  }

  nproc = old_nproc;

  if (i_am_the_master)
    printf("----------------\n");
  MPI_Barrier(MPI_COMM_WORLD);
}

void test_csr_inverse() {

  int IA[] = {0, 3, 6, 9};
  int JA[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  double AA[] = {1, 2, 3, 15, 5, 6, 7, 13, 9};
  int N = 3;

  CSRMatrix A;
  A.IA = IA;
  A.JA = JA;
  A.A = AA;
  A.N = N;
  CSRMatrix D = csr_diag_inverse(A);

  int ID[] = {0, 1, 2, 3};
  int JD[] = {0, 1, 2};
  double DD[] = {1.0, 0.2, 0.11111111};
  int ND = 3;
  double eps = 1e-7;
  int i;

  if (i_am_the_master)
    print_csr(D);

  for (i = 0; i < sizeof(ID)/sizeof(int); i++) {
    my_assert(ID[i] == D.IA[i], "ID are not equal!");
  }

  for (i = 0; i < sizeof(JD)/sizeof(int); i++) {
    my_assert(JD[i] == D.JA[i], "JD are not equal!");
  }

  for (i = 0; i < sizeof(DD)/sizeof(double); i++) {
    if (i_am_the_master)
      printf("%lf vs %lf\n", D.A[i], DD[i]);
    my_assert(fabs(D.A[i] - DD[i]) < eps, "DD's are not equal!");
  }

  my_assert(D.N == ND, "Dimension is incorrect!");

  csr_free(&D);

  if (i_am_the_master)
    printf("----------------\n");
  MPI_Barrier(MPI_COMM_WORLD);
}

void test_big_slau() {
  int max_values[] = {100, 100, 100};

  // volume of the cube
  int N = 1;
  for (int i = 0; i < N_DIMENSIONS; i++) {
    N *= max_values[i];
  }
  if (i_am_the_master)
    printf("N = %d\n", N);

  old_nproc = nproc;
  int P[N_DIMENSIONS];
  int grid_nproc = make_grid(max_values, P);

  if (myrank < grid_nproc) {
    nproc = grid_nproc;

    int n_local;
    int *Rows;
    int p[N_DIMENSIONS];
    int ps[N_DIMENSIONS];
    int n_start[N_DIMENSIONS];
    grid_params(myrank, max_values, P, p, ps, n_start, &n_local);
    CSRMatrix A = csr_matrix(max_values, n_start, ps, n_local, &Rows);

    int *Part = get_part(max_values, P);

    initialize(A, Rows, Part, N);
    Vector XX = make_vector(N);
    Vector BB = make_vector(N);
    for (int i = 0; i < BB.ni; i++) {
      BB.data[i] = sin(Rows[i] + 1);
    }
    do_exchange(BB);

    double eps = 1e-4;
    double tol = eps*eps;
    int maxit = 50;
    int info = 1;
    int i = bicgstab_solve_par(A, BB, XX, N, tol, maxit, info);
    Vector Y = make_vector(N);
    SpMV_par(A, XX, Y);
    for (i = 0; i < Y.ni; i++) {
      // printf("%lf vs %lf\n", Y.data[i], BB.data[i]);
      my_assert(fabs(Y.data[i] - BB.data[i]) < eps, "Solution is incorrect!");
    }

    finalize();
    free_vector(&BB);
    free_vector(&XX);
    free_vector(&Y);
    csr_free(&A);
  }

  nproc = old_nproc;

  if (i_am_the_master)
    printf("----------------\n");
  MPI_Barrier(MPI_COMM_WORLD);
}

void test_dot() {
  int max_values[N_DIMENSIONS] = {10, 10, 10};
  // volume of the cube
  int N = 1;
  for (int i = 0; i < N_DIMENSIONS; i++) {
    N *= max_values[i];
  }

  int P[N_DIMENSIONS];
  old_nproc = nproc;
  int grid_nproc = make_grid(max_values, P);

  if (myrank < grid_nproc) {
    int n_local;
    int p[N_DIMENSIONS];
    int ps[N_DIMENSIONS];
    int n_start[N_DIMENSIONS];
    grid_params(myrank, max_values, P, p, ps, n_start, &n_local);

    int *Rows = NULL;
    CSRMatrix A = csr_matrix(max_values, n_start, ps, n_local, &Rows);

    int *Part = get_part(max_values, P);

    initialize(A, Rows, Part, N);

    Vector XX = make_vector(N);
    Vector YY = make_vector(N);

    for (int i = 0; i < ni; i++) {
      XX.data[i] = sin(Rows[i] + 1);
      YY.data[i] = cos(Rows[i] + 1);
    }

    double d = dot_par(XX, YY);
    double true_ans = 0.452019;
    double eps = 1e-6;
    if (i_am_the_master) {
      printf("dot(XX, YY) for X = sin(i), Y = cos(i), N = %d: %lf\n", N, d);
      printf("True ans: %f\n", true_ans);
    }

    my_assert(fabs(d - true_ans) < eps, "Not equal!");

    free_vector(&XX);
    free_vector(&YY);
    csr_free(&A);

    finalize();
  }

  nproc = old_nproc;

  if (i_am_the_master)
    printf("----------------\n");
  MPI_Barrier(MPI_COMM_WORLD);
}

void test_axpby() {
  int max_values[N_DIMENSIONS] = {5,10,2};
  // volume of the cube
  int N = 1;
  for (int i = 0; i < N_DIMENSIONS; i++) {
    N *= max_values[i];
  }

  int P[N_DIMENSIONS];
  old_nproc = nproc;
  int grid_nproc = make_grid(max_values, P);

  if (myrank < grid_nproc) {
    int n_local;
    int p[N_DIMENSIONS];
    int ps[N_DIMENSIONS];
    int n_start[N_DIMENSIONS];
    grid_params(myrank, max_values, P, p, ps, n_start, &n_local);

    int *Rows = NULL;
    CSRMatrix A = csr_matrix(max_values, n_start, ps, n_local, &Rows);

    int *Part = get_part(max_values, P);

    initialize(A, Rows, Part, N);

    Vector XX = make_vector(N);
    Vector YY = make_vector(N);
    for (int i = 0; i < ni; i++) {
      XX.data[i] = sin(Rows[i] + 1);
      YY.data[i] = cos(Rows[i] + 1);
    }
    do_exchange(XX);
    do_exchange(YY);

    axpby_par(XX, YY, 1.0, -1.0);
    double s = 0;
    for (int i = 0; i < ni; i++) {
      s += XX.data[i];
    }

    double s_all = 0;
    MPI_Allreduce(&s, &s_all, 1, MPI_DOUBLE, MPI_SUM, comm_grid);

    double true_ans = 0.405118;
    double eps = 1e-6;
    if (i_am_the_master) {
      printf("axpby(XX, YY, a, b) = %lf\n", s_all);
      printf("True ans: %lf\n", true_ans);
    }

    my_assert(fabs(s_all - true_ans) < eps, "Not equal!");

    free_vector(&XX);
    free_vector(&YY);
    csr_free(&A);

    finalize();
  }

  nproc = old_nproc;

  if (i_am_the_master)
    printf("----------------\n");
  MPI_Barrier(MPI_COMM_WORLD);
}

void test_smpv() {
  int max_values[N_DIMENSIONS] = {10, 10, 10};
  // volume of the cube
  int N = 1;
  for (int i = 0; i < N_DIMENSIONS; i++) {
    N *= max_values[i];
  }

  int P[N_DIMENSIONS];
  old_nproc = nproc;
  int grid_nproc = make_grid(max_values, P);

  if (myrank < grid_nproc) {
    int n_local;
    int p[N_DIMENSIONS];
    int ps[N_DIMENSIONS];
    int n_start[N_DIMENSIONS];
    grid_params(myrank, max_values, P, p, ps, n_start, &n_local);

    int *Rows = NULL;
    CSRMatrix A = csr_matrix(max_values, n_start, ps, n_local, &Rows);

    int *Part = get_part(max_values, P);

    initialize(A, Rows, Part, N);

    Vector XX = make_vector(N);
    Vector YY = make_vector(N);
    for (int i = 0; i < ni; i++) {
      XX.data[i] = sin(Rows[i] + 1);
    }
    do_exchange(XX);
    SpMV_par(A, XX, YY);

    double s = 0;
    for (int i = 0; i < ni; i++) {
      s += YY.data[i];
    }
    double s_all = 0;
    MPI_Allreduce(&s, &s_all, 1, MPI_DOUBLE, MPI_SUM, comm_grid);
    double true_ans = 4.042218;
    double eps = 1e-6;
    if (i_am_the_master) {
      printf("spMV(A, XX, YY) = %lf\n", s_all);
      printf("True ans: %lf\n", true_ans);
    }

    my_assert(fabs(s_all - true_ans) < eps, "Not equal!");

    free_vector(&XX);
    free_vector(&YY);
    csr_free(&A);

    finalize();
  }

  nproc = old_nproc;

  if (i_am_the_master)
    printf("----------------\n");
  MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  i_am_the_master = myrank == 0;

  if (deg_2(nproc) == ERROR) {
    if (i_am_the_master) {
      fprintf(stderr, "%s\n", WRONG_NPROC);
      exit(FAIL_NPROC);
    }
  }

  if (i_am_the_master) {
    printf("Hello, world!\n");
    printf("Comm size: %d\n", nproc);
    printf("Comm rank: %d\n", myrank);
    printf("----------------\n");
  }

  omp_set_num_threads(1);

  test_create_csr();
  test_csr_inverse();
  test_big_slau();
  test_dot();
  test_axpby();
  test_smpv();
  MPI_Finalize();
  return 0;
}