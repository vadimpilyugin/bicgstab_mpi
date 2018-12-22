#pragma once

#define TEST_DOT 0
#define TEST_SPMV 1
#define TEST_AXPBY 2
#define TEST_SOLVER 3
#define N_OPS 4

void init(int argc, char **argv);
void fin();
double test_op(int op_num);