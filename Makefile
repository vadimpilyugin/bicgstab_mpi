main: main.o algorithm_par.o create_csr.o
	gcc -o main main.o algorithm_par.o create_csr.o -Wall -std=gnu11 -lm -fopenmp -O3 -march=native -mtune=native 

main2: main2.c
	mpicc -Wpedantic -Wall -o main2 main2.c measure.c division.c create_csr.c algorithm_par.c core.c -lm -fopenmp -O3 -march=native -mtune=native 

division: division.c
	gcc -o division division.c -lm

main.o: main.c
	gcc -o main.o -c main.c -Wall -std=gnu11 -O3 -march=native -mtune=native 

create_csr.o: create_csr.c
	gcc -o create_csr.o -c create_csr.c -Wall -std=gnu11 -O3 -march=native -mtune=native 

algorithm.o: algorithm.c
	gcc -o algorithm.o -c algorithm.c -Wall -std=gnu11 -O3 -march=native -mtune=native 

algorithm_par.o: algorithm_par.c
	gcc -o algorithm_par.o -c algorithm_par.c -Wall -std=gnu11 -fopenmp -O3 -march=native -mtune=native 

test.o: test.c
	gcc -o test.o -c test.c -Wall -std=gnu11

test: test.c
	mpicc -Wpedantic -Wall -o test test.c division.c create_csr.c algorithm_par.c core.c -lm -fopenmp -O3 -march=native -mtune=native

.PHONY: clean
clean:
	rm -f *.o
	rm -f main
	rm -f main2
	rm -f test