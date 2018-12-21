#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "algorithm_par.h"
#include "division.h"
#include "create_csr.h"
#include "core.h"
#include "mpi.h"
#include "omp.h"

#define MASTER 0

#define ARG_NX 1
#define ARG_NY 2
#define ARG_NZ 3
#define ARG_TOL 4
#define ARG_MAXIT 5
#define ARG_NT 6
#define ARG_QA 7
#define ARG_N_DOT 8
#define ARG_N_SPMV 9
#define ARG_N_AXPBY 10
#define ARG_N_SOLVER 11
#define N_ARGS 12

#define FAIL_ARGS 1
#define FAIL_NPROC 2
#define EXIT_GRID 3

#define TEST_AXPBY 0
#define TEST_SPMV 1
#define TEST_DOT 2
#define TEST_SOLVER 3
#define N_OPS 4

int myrank, nproc, i_am_the_master;
MPI_Comm comm_grid;

static int max_values[N_DIMENSIONS];
static double tol;
static int maxit;
static int nt;
static int qa;

static int n_tests[N_OPS];

static int *Rows;
static int *Part;
static int N;

static CSRMatrix A;
static Vector BB;
static Vector XX;
static Vector YY;

const char *USAGE_S = "Usage: main [NX] [NY] [NZ] [TOL] [MAXIT] [THREADS] [QA] [N_DOT] [N_SPMV] N_AXPBY] [N_SOLVER]";

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
    usage();
  }

  // if nproc is not a power of 2 then exit
  if (deg_2(nproc) == ERROR) {
    if (i_am_the_master) {
      fprintf(stderr, "%s\n", "Number of processes is not a power of 2");
      exit(FAIL_NPROC);
    }
  }

  if (i_am_the_master) {
    printf("Hello, world!\n");
    printf("Comm size: %d\n", nproc);
    printf("Comm rank: %d\n", myrank);
  }

  // read command-line arguments
  max_values[X_DIM] = atoi(argv[ARG_NX]);
  max_values[Y_DIM] = atoi(argv[ARG_NY]);
  max_values[Z_DIM] = atoi(argv[ARG_NZ]);
  tol               = strtod(argv[ARG_TOL], NULL);
  maxit             = atoi(argv[ARG_MAXIT]);
  nt                = atoi(argv[ARG_NT]);
  qa                = atoi(argv[ARG_QA]);
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
    printf("Cube: %dx%dx%d\n", max_values[X_DIM], max_values[Y_DIM], max_values[Z_DIM]);
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
  
  if (i_am_the_master) {
    printf("Division: %dx%dx%d process grid, processes [0..%d]\n", 
      P[X_DIM], P[Y_DIM], P[Z_DIM], grid_nproc-1);
  }

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

  printf("Hi! I am process %d! My coords are [%d, %d, %d]! My cube is %dx%dx%d! I start from %d,%d,%d!\n", myrank, 
    p[X_DIM], p[Y_DIM], p[Z_DIM], ps[X_DIM], ps[Y_DIM], ps[Z_DIM], 
    n_start[X_DIM], n_start[Y_DIM], n_start[Z_DIM]);

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