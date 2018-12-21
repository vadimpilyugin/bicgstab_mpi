#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "mpi.h"
#include "division.h"

#define X_DIM 0
#define Y_DIM 1
#define Z_DIM 2

extern int nproc;
extern MPI_Comm comm_grid;

// int deg_2(int n) {
//   int i = 0;
//   while (n > 1) {
//     if (n % 2 != 0) {
//       return ERROR;
//     }
//     i += 1;
//     n /= 2;
//   }
//   return i;
// }

int max_deg(int n) {
  return log(n) / log(2);
}

int deg_2(int n) {
  int deg = max_deg(n);
  if ((1 << deg) != n) {
    return -1;
  }
  return deg;
}

int max_ind(int a[N_DIMENSIONS]) {
  int max_i = 0;
  int max = a[0];
  for (int i = 1; i < N_DIMENSIONS; i++) {
    if (a[i] > max) {
      max = a[i];
      max_i = i;
    }
  }
  return max_i;
}

int *g_ar;

int comp(const void * elem1, const void * elem2) 
{
    int f = *((int *) elem1);
    int s = *((int *) elem2);
    if (g_ar[f] > g_ar[s]) return -1;
    if (g_ar[f] < g_ar[s]) return 1;
    return 0;
}

void sort_ind(int a[N_DIMENSIONS], int ind[N_DIMENSIONS]) {
  for(int i = 0; i < N_DIMENSIONS; i++) {
    ind[i] = i;
  }
  g_ar = (int *)a;
  qsort((int *)ind, N_DIMENSIONS, sizeof(int), comp);
}

void divide_cube(int max_values[N_DIMENSIONS], int division[N_DIMENSIONS], int nproc) {
  int degree = deg_2(nproc);
  int i;
  for (i = 0; i < N_DIMENSIONS; i++) {
    division[i] = max_deg(max_values[i]);
  }
  int indices[N_DIMENSIONS];
  sort_ind(max_values, indices);
  int j = 0;
  for (i = 0; i < degree; i++) {
    int j_begin = j;
    // printf("Next j: %d\n", j);
    while (division[indices[j]] == 0) {
      j = (j + 1) % N_DIMENSIONS;
      if (j == j_begin) {
        break;
      }
    }
    if (division[indices[j]] == 0) {
      break;
    } else {
      division[indices[j]]--;
      j = (j + 1) % N_DIMENSIONS;
    }
  }
  for (i = 0; i < N_DIMENSIONS; i++) {
    division[i] = 1 << (max_deg(max_values[i]) - division[i]);
  }
}

int part_size(int N, int P, int p) {
  return N / P + (p < N % P);
}

int min(int x, int y) {
  if (x < y) {
    return x;
  } else {
    return y;
  }
}

int start_pos(int N, int P, int p) {
  return p * (N / P) + min(p, N % P);
}

int glob_row(int max_values[N_DIMENSIONS], int nx, int ny, int nz) {
  return nx + ny * max_values[X_DIM] + nz * max_values[X_DIM] * max_values[Y_DIM];
}

int * get_part(int max_values[N_DIMENSIONS], int P[N_DIMENSIONS]) {
  // volume of the cube
  int N = 1;
  for (int i = 0; i < N_DIMENSIONS; i++) {
    N *= max_values[i];
  }
  int *Part = (int *)malloc(N * sizeof(int));

  // number of processes
  int nproc = 1;
  for (int i = 0; i < N_DIMENSIONS; i++) {
    nproc *= P[i];
  }

  for (int rank = 0; rank < nproc; rank++) {
    // my coords in a process grid
    int p[N_DIMENSIONS];
    int mod = 1;
    int div = 1;
    for (int i = 0; i < N_DIMENSIONS; i++) {
      mod *= P[i];
      p[i] = (rank % mod) / div;
      div *= P[i];
    }

    // dimensions of the small cube
    int ps[N_DIMENSIONS];
    for (int i = 0; i < N_DIMENSIONS; i++) {
      ps[i] = part_size(max_values[i], P[i], p[i]);
    }

    // starting position
    int n_start[N_DIMENSIONS];
    for (int i = 0; i < N_DIMENSIONS; i++) {
      n_start[i] = start_pos(max_values[i], P[i], p[i]);
    }
    
    for (int nz = n_start[Z_DIM]; nz < n_start[Z_DIM] + ps[Z_DIM]; nz++)
      for (int ny = n_start[Y_DIM]; ny < n_start[Y_DIM] + ps[Y_DIM]; ny++)
        for (int nx = n_start[X_DIM]; nx < n_start[X_DIM] + ps[X_DIM]; nx++) {
          Part[glob_row(max_values, nx,ny,nz)] = rank;
        }
  }

  return Part;
}

void grid_params(int rank, int max_values[N_DIMENSIONS], int P[N_DIMENSIONS], 
  int p[N_DIMENSIONS], int ps[N_DIMENSIONS], int n_start[N_DIMENSIONS], int *_n_local) {

  // my coords in a process grid
  int mod = 1;
  int div = 1;
  for (int i = 0; i < N_DIMENSIONS; i++) {
    mod *= P[i];
    p[i] = (rank % mod) / div;
    div *= P[i];
  }

  // dimensions of the small cube
  int n_local = 1; // volume of the small cube
  for (int i = 0; i < N_DIMENSIONS; i++) {
    ps[i] = part_size(max_values[i], P[i], p[i]);
    n_local *= ps[i];
  }

  // starting position
  for (int i = 0; i < N_DIMENSIONS; i++) {
    n_start[i] = start_pos(max_values[i], P[i], p[i]);
  }

  *_n_local = n_local;
}

int make_grid(int max_values[N_DIMENSIONS], int P[N_DIMENSIONS]) {
  divide_cube(max_values, P, nproc);
  int grid_nproc = 1;
  for (int i = 0; i < N_DIMENSIONS; i++) {
    grid_nproc *= P[i];
  }

  int *process_ranks = (int*) malloc(grid_nproc * sizeof(int));
  for(int i = 0; i < grid_nproc; i++)
    process_ranks[i] = i;

  MPI_Group group_world;
  MPI_Group grid_group;
  MPI_Comm_group(MPI_COMM_WORLD, &group_world);
  MPI_Group_incl(group_world, grid_nproc, process_ranks, &grid_group);
  MPI_Comm_create(MPI_COMM_WORLD, grid_group, &comm_grid);
  free(process_ranks);

  return grid_nproc;
}

// int main() {
//   printf("Deg 2(1024) = %d\n", deg_2(1024));
//   printf("Deg 2(768) = %d\n", deg_2(768));
//   printf("Deg 2(1) = %d\n", deg_2(1));
//   printf("Deg 2(2) = %d\n", deg_2(2));
//   printf("Deg 2(3) = %d\n", deg_2(3));
//   printf("Max deg(50) = %d\n", max_deg(50));
//   printf("Max deg(20) = %d\n", max_deg(20));
//   printf("Max deg(10) = %d\n", max_deg(10));
//   printf("Max deg(1) = %d\n", max_deg(1));

//   int max_values[N_DIMENSIONS] = {10,50,100};
//   int division[N_DIMENSIONS];
//   int nproc = 256;
//   printf("Dividing cube (%dx%dx%d) into cubes:\n", max_values[0], max_values[1], max_values[2]);
//   printf("Nproc = %d\n", nproc);
//   divide_cube(max_values, division, nproc);
//   printf("Divided processes into grid: %dx%dx%d\n", division[0], division[1], division[2]);

//   printf("Sorting indices:\n");
//   sort_ind(max_values, division);
//   printf("Indices: %d,%d,%d\n", division[0], division[1], division[2]);
// }