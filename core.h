#pragma once

#include "create_csr.h"

int ni;
int nh;

void do_exchange(Vector V);
void initialize(CSRMatrix A, int *Rows, int *Part, int N);
void finalize();