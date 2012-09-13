all: DirectPut

a.out: main.c spi.h
	mpixlc_r -O3 main.c -lpthread

DirectPut: DirectPut.c DirectPut.h Makefile
	mpixlc_r -O3 -qsmp=omp DirectPut.c -o DirectPut -lpthread
