all: a.out DirectPut

a.out: main.c spi.h
	mpixlc_r -O3 main.c -lpthread

DirectPut: DirectPut.c spi.h msg_common.c
	mpixlc_r -O3 DirectPut.c msg_common.c -o DirectPut -lpthread
