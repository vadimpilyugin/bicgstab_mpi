#include <stdio.h>
#include "measure.h"
#include "omp.h"

extern int i_am_the_master;
extern int nproc;
extern const char *op_str[N_OPS];
const char *format = "%s\ttime=%6.3fs NTR=%d nproc=%d\n";

void tests() {
  // for each number of threads measure time taken by each operation
  // N, P are given from command line
  const int NTR = omp_get_num_procs();
  for(int ntr = 1; ntr <= NTR; ntr *= 2) {
    if (i_am_the_master)
      printf("\ntesting parallel ops for ntr=%d:\n", ntr);
    omp_set_num_threads(ntr);
    for (int op_num = 0; op_num < N_OPS; op_num++) {
      double t = test_op(op_num);
      if (i_am_the_master)
        printf(format, op_str[op_num], t, ntr, nproc);
    }
  }
}

int main(int argc, char **argv) {
  init(argc, argv);
  tests();
  fin();
  return 0;
}