#include "stdlib.h"
#include "assert.h"
#include "math.h"
#include "stdio.h"
#include "division.h"

#include "create_csr.h"
#include "core.h"

extern int ni;
extern int nh;

Vector make_vector(int size) {
  Vector V;
  V.size = size;
  V.ni = ni;
  V.nh = nh;
  V.data = (double *)malloc((ni + nh) * sizeof(double));
  return V;
}

void free_vector(Vector *X) {
  assert(X -> data != NULL);
  free(X -> data);
  X -> size = 0;
  X -> ni = 0;
  X -> nh = 0;
  X -> data = NULL;
}

static void fill_neighbours(int nx, int ny, int nz, int i,
  int neighbours[MAX_NEIGHBOURS], int max_values[N_DIMENSIONS]) {

  neighbours[0] = nz == 0 ? -1 : i - max_values[X_DIM] * max_values[Y_DIM];
  neighbours[1] = ny == 0 ? -1 : i - max_values[X_DIM];
  neighbours[2] = nx == 0 ? -1 : i - 1;
  neighbours[3] = i;
  neighbours[4] = nx == max_values[X_DIM] - 1 ? -1 : i + 1;
  neighbours[5] = ny == max_values[Y_DIM] - 1 ? -1 : i + max_values[X_DIM];
  neighbours[6] = nz == max_values[Z_DIM] - 1 ? -1 : i + max_values[X_DIM] * max_values[Y_DIM];
}

static void csr_init_malloc(CSRMatrix *A, int N, int *curr_cap) {
  A -> IA = (int *)malloc((N + 1) * sizeof(int));
  A -> JA = (int *)malloc(DEFAULT_SIZE * sizeof(int));
  A -> A = (double *)malloc(DEFAULT_SIZE * sizeof(double));
  *curr_cap = DEFAULT_SIZE;
}

static void csr_realloc(CSRMatrix *A, int *curr_cap) {
  *curr_cap *= *curr_cap;

  void *tmp = (int *)realloc(A -> JA, (*curr_cap) * sizeof(int));
  assert(tmp != NULL);
  A -> JA = tmp;

  tmp = (double *)realloc(A -> A, (*curr_cap) * sizeof(double));
  assert(tmp != NULL);
  A -> A = tmp;
}

void csr_free(CSRMatrix *A) {
  free(A -> IA);
  free(A -> JA);
  free(A -> A);

  A -> N = 0;
  A -> IA = NULL;
  A -> JA = NULL;
  A -> A = NULL;
}

CSRMatrix csr_matrix(int max_values[N_DIMENSIONS], int n_start[N_DIMENSIONS], 
  int ps[N_DIMENSIONS], int n_local, int **_Rows) {

  int *Rows = *_Rows;

  CSRMatrix A;
  A.N = n_local;
  int curr_cap;
  csr_init_malloc(&A, n_local, &curr_cap);
  int M = 0; // number of non-zero elements in the local matrix
  int neighbours[MAX_NEIGHBOURS]; // indices of neighbouring cells
  int i = 0; // row number in local numeration
  Rows = (int *)malloc(n_local * sizeof(int));

  for (int nz = n_start[Z_DIM]; nz < n_start[Z_DIM] + ps[Z_DIM]; nz++)
    for (int ny = n_start[Y_DIM]; ny < n_start[Y_DIM] + ps[Y_DIM]; ny++)
      for (int nx = n_start[X_DIM]; nx < n_start[X_DIM] + ps[X_DIM]; nx++) {
        A.IA[i] = M;

        int global_row = glob_row(max_values, nx,ny,nz);
        fill_neighbours(nx, ny, nz, global_row, neighbours, max_values);
        Rows[i] = global_row;

        int middle = -1; // index in JA corresponding to element A[glob_row, glob_row]
        double s = 0; // row sum except the diagonal
        for (int j = 0; j < MAX_NEIGHBOURS; j++) {
          if (neighbours[j] >= 0) {
            if (M == curr_cap) {
              csr_realloc(&A, &curr_cap);
            }
            A.JA[M] = neighbours[j];
            if (global_row == neighbours[j]) {
              A.A[M] = 0;
              middle = M;
            } else {
              double sin_ij = sin(global_row + neighbours[j] + 1);
              A.A[M] = sin_ij;
              s += fabs(sin_ij);
            }
            M += 1;
          }
        }
        assert(middle != -1);
        A.A[middle] = 1.1*s;
        i++;
      }

  A.IA[n_local] = M;
  *_Rows = Rows;

  return A;
}

static void csr_fixed_malloc(CSRMatrix *A, int N, int M) {
  A -> IA = (int *)malloc((N + 1) * sizeof(int));
  A -> JA = (int *)malloc(M * sizeof(int));
  A -> A = (double *)malloc(M * sizeof(double));
  A -> N = N;
  A -> IA[N] = M; // M non-zero elements in A
}

void print_csr(CSRMatrix A) {
  int i;
  printf("\nIA: [");
  for (i = 0; i <= A.N; i++) {
    printf("%d, ", A.IA[i]);
  }
  printf("]\n");
  printf("JA: [");
  for (i = 0; i < A.IA[A.N]; i++) {
    printf("%d, ", A.JA[i]);
  }
  printf("]\n");
  printf("A: [");
  for (i = 0; i < A.IA[A.N]; i++) {
    printf("%lf, ", A.A[i]);
  }
  printf("]\n\n");
}

CSRMatrix csr_diag_inverse(CSRMatrix A) {
  CSRMatrix D;
  csr_fixed_malloc(&D, A.N, A.N); // Matrix of size N with N non-zero elements
  int i;
  for (i = 0; i < A.N; i++) {
    D.IA[i] = i;
    D.JA[i] = i;
    int j;
    for (j = A.IA[i]; j < A.IA[i+1]; j++) {
      if (A.JA[j] == i) {
        D.A[i] = 1/A.A[j];
      }
    }
  }
  return D;
}
