#pragma once

#define N_DIMENSIONS 3
#define ERROR -1
#define OK 0

int deg_2(int n);
void divide_cube(int max_values[N_DIMENSIONS], int division[N_DIMENSIONS], int nproc);
int part_size(int N, int P, int p);
int start_pos(int N, int P, int p);
int glob_row(int max_values[N_DIMENSIONS], int nx, int ny, int nz);
int * get_part(int max_values[N_DIMENSIONS], int P[N_DIMENSIONS]);
void grid_params(int rank, int max_values[N_DIMENSIONS], int P[N_DIMENSIONS], 
  int p[N_DIMENSIONS], int ps[N_DIMENSIONS], int n_start[N_DIMENSIONS], int *_n_local);
int make_grid(int max_values[N_DIMENSIONS], int P[N_DIMENSIONS]);