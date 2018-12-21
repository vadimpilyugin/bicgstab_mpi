#pragma once

#define N_DIMENSIONS 3
#define MAX_NEIGHBOURS 7
#define DEFAULT_SIZE 8
#define X_DIM 0
#define Y_DIM 1
#define Z_DIM 2

typedef struct CSRMatrix {
  int *IA;
  int *JA;
  double *A;
  int N;
} CSRMatrix;

typedef struct Vector {
  double *data;
  int size;
  int ni;
  int nh;
} Vector;

void csr_free(CSRMatrix *A);
CSRMatrix csr_matrix(int max_values[N_DIMENSIONS], int n_start[N_DIMENSIONS], int ps[N_DIMENSIONS], int n_local, int **_Rows);
CSRMatrix csr_diag_inverse(CSRMatrix A);
void print_csr(CSRMatrix A);
Vector make_vector(int size);
void free_vector(Vector *X);
