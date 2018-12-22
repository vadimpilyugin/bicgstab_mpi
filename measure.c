#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "algorithm_par.h"
#include "division.h"
#include "create_csr.h"
#include "core.h"
#include "measure.h"
#include "mpi.h"
#include "omp.h"

#define MASTER 0

#define ARG_NX 1
#define ARG_NY 2
#define ARG_NZ 3
#define ARG_TOL 4
#define ARG_MAXIT 5
#define ARG_NT 6
#define ARG_N_DOT 7
#define ARG_N_SPMV 8
#define ARG_N_AXPBY 9
#define ARG_N_SOLVER 10
#define N_ARGS 11

#define FAIL_ARGS 1
#define FAIL_NPROC 2
#define EXIT_GRID 3

int myrank, nproc, i_am_the_master;
MPI_Comm comm_grid;

static int max_values[N_DIMENSIONS];
static double tol;
static int maxit;
int n_threads;

static int *Rows;
static int *Part;
int N;

int n_tests[N_OPS];
const char *op_str[N_OPS] = {
  "dot",
  "SpMV",
  "axpby",
  "solver",
};

static CSRMatrix A;
static Vector BB;
static Vector XX;
static Vector YY;

const char *USAGE_S = "Usage: main [NX] [NY] [NZ] [TOL] [MAXIT] [THREADS] [N_DOT] [N_SPMV] N_AXPBY] [N_SOLVER]";

static void usage() {
  if (i_am_the_master)
    fprintf(stderr, "%s\n", USAGE_S);
  MPI_Finalize();
  exit(FAIL_ARGS);
}

void init(int argc, char **argv) {
  // initialize MPI Library
  MPI_Init(&argc, &argv);

  // get nproc and myrank
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  i_am_the_master = myrank == 0;
  if (argc != N_ARGS) {
    for (int i = 0; i < argc; i++) {
      fprintf(stderr, "Arg %d: %s, N_ARGS=%d\n", i, argv[i], N_ARGS);
    }
    usage();
  }

  // if nproc is not a power of 2 then exit
  if (deg_2(nproc) == -1) {
    if (i_am_the_master) {
      fprintf(stderr, "Nproc=%d\n", nproc);
      fprintf(stderr, "%s\n", "Number of processes is not a power of 2");
      exit(FAIL_NPROC);
    }
  }

  // read command-line arguments
  max_values[X_DIM] = atoi(argv[ARG_NX]);
  max_values[Y_DIM] = atoi(argv[ARG_NY]);
  max_values[Z_DIM] = atoi(argv[ARG_NZ]);
  tol               = strtod(argv[ARG_TOL], NULL);
  maxit             = atoi(argv[ARG_MAXIT]);
  n_threads         = atoi(argv[ARG_NT]);
  n_tests[TEST_DOT]     = atoi(argv[ARG_N_DOT]);
  n_tests[TEST_SPMV]    = atoi(argv[ARG_N_SPMV]);
  n_tests[TEST_AXPBY]   = atoi(argv[ARG_N_AXPBY]);
  n_tests[TEST_SOLVER]  = atoi(argv[ARG_N_SOLVER]);

  // volume of the cube
  N = 1;
  for (int i = 0; i < N_DIMENSIONS; i++) {
    N *= max_values[i];
  }

  if (i_am_the_master) {
    printf("Testing BiCGSTAB solver for a 3D grid domain\n");
    printf("nx=%d ny=%d nz=%d tol=%lf maxit=%d n_threads=%d nproc=%d\n", 
      max_values[X_DIM], max_values[Y_DIM], max_values[Z_DIM], 
      tol, maxit, n_threads, nproc);
    printf("n_dot=%d n_spmv=%d n_axpby=%d n_solver=%d\n\n", 
      n_tests[TEST_DOT], n_tests[TEST_SPMV], n_tests[TEST_AXPBY], n_tests[TEST_SOLVER]);
    printf("N   = %d (Nx=%d, Ny=%d, Nz=%d\n", N, 
      max_values[X_DIM], max_values[Y_DIM], max_values[Z_DIM]);
    printf("Aij = sin(i+j+1), i != j\n");
    printf("Aii = 1.1*sum(fabs(Aij))\n");
    printf("Bi  = sin(i+1)\n");
    printf("tol = %.10e\n\n", tol);
  }

  // calculate process grid dimensions
  int P[N_DIMENSIONS];
  int grid_nproc = make_grid(max_values, P);
  
  // if process is outside of grid then it exits
  if (myrank >= grid_nproc) {
    printf("--- Exiting: rank=%d, max_rank=%d\n", myrank, grid_nproc);
    exit(EXIT_GRID);
  }

  nproc = grid_nproc;

  // get other grid parameters
  int n_local;  // small cube volume
  int p[N_DIMENSIONS]; // process coordinates
  int ps[N_DIMENSIONS]; // small cube dimensions
  int n_start[N_DIMENSIONS]; // coordinates of starting point
  grid_params(myrank, max_values, P, p, ps, n_start, &n_local);

  // calculate matrix from a given cube
  A = csr_matrix(max_values, n_start, ps, n_local, &Rows);

  // calculate part array based on two grids
  Part = get_part(max_values, P);

  // initialize the core with three pieces of data (N is the length of Part)
  initialize(A, Rows, Part, N);

  // prepare vectors for tests
  BB = make_vector(N); // rhs
  XX = make_vector(N); // solution
  YY = make_vector(N); // vec1
  for (int i = 0; i < ni; i++) {
    BB.data[i] = sin(Rows[i] + 1);
    XX.data[i] = sin(Rows[i] + 1);
    YY.data[i] = cos(Rows[i] + 1);
  }
  do_exchange(BB);
  do_exchange(XX);
  do_exchange(YY);
}

void fin() {
  finalize();
  free(Part);
  free(Rows);
  csr_free(&A);
  free_vector(&BB);
  free_vector(&XX);
  free_vector(&YY);
  MPI_Finalize();
}

double test_op(int op_num) {
  double start = MPI_Wtime();
  if (i_am_the_master)
    printf("Starting tests for %s (%d): loop %d times\n",
    op_str[op_num], op_num, n_tests[op_num]);

  for (int i = 0; i < n_tests[op_num]; i++) {
    if (op_num == TEST_AXPBY) {
      axpby_par(XX, YY, 1.1, 0.99);
    } else if (op_num == TEST_SPMV) {
      SpMV_par(A, XX, YY);
    } else if (op_num == TEST_DOT) {
      dot_par(XX, YY);
    } else if (op_num == TEST_SOLVER) {
      bicgstab_solve_par(A, BB, XX, N, tol, maxit, 0);
    }
  }
  return MPI_Wtime() - start;
}