#include <stdio.h>
#include "measure.h"
#include "omp.h"

extern int i_am_the_master;

int main(int argc, char **argv) {
  init(argc, argv);
  omp_set_num_threads(2);
  double t = test_op(TEST_DOT);
  if (i_am_the_master)
    printf("Dot: %lf seconds\n", t);
  t = test_op(TEST_SPMV);
  if (i_am_the_master)
    printf("SpMV: %lf seconds\n", t);
  t = test_op(TEST_SOLVER);
  if (i_am_the_master)
    printf("Solver: %lf seconds\n", t);
  fin();
  return 0;
}