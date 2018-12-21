main: main.c measure.c division.c create_csr.c algorithm_par.c core.c
	mpicc -Wpedantic -Wall -o main main.c measure.c division.c create_csr.c algorithm_par.c core.c -lm -fopenmp -O3 -march=native -mtune=native 

test: test.c division.c create_csr.c algorithm_par.c core.c
	mpicc -Wpedantic -Wall -o test test.c division.c create_csr.c algorithm_par.c core.c -lm -fopenmp -O3 -march=native -mtune=native

.PHONY: clean
clean:
	rm -f *.o
	rm -f main
	rm -f main2
	rm -f test