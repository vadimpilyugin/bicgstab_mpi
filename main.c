#include <stdio.h>
#include "measure.h"
#include "omp.h"

#define LOOP_THREADS 12345

extern int i_am_the_master;
extern int nproc;
extern int N;
extern int n_threads;
extern const char *op_str[N_OPS];
const char *format = "%s\ttime=%6.3fs NTR=%d nproc=%d N=%d\n";

void tests() {
  int t_begin;
  int t_end;
  if (n_threads == LOOP_THREADS) {
    t_begin = 1;
    t_end = omp_get_num_procs();
  } else {
    t_begin = n_threads;
    t_end = n_threads;
  }
  // for each number of threads measure time taken by each operation
  // N, P are given from command line
  for(int ntr = t_begin; ntr <= t_end; ntr *= 2) {
    if (i_am_the_master)
      printf("\ntesting parallel ops for ntr=%d:\n", ntr);
    omp_set_num_threads(ntr);
    for (int op_num = 0; op_num < N_OPS; op_num++) {
      double t = test_op(op_num);
      if (i_am_the_master)
        printf(format, op_str[op_num], t, ntr, nproc, N);
    }
  }
}

int main(int argc, char **argv) {
  init(argc, argv);
  tests();
  fin();
  return 0;
}