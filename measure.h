#pragma once

#define TEST_AXPBY 0
#define TEST_SPMV 1
#define TEST_DOT 2
#define TEST_SOLVER 3

void init(int argc, char **argv);
void fin();
double test_op(int op_num);