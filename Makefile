a.out: main.c spi.h
	mpixlc_r -O3 main.c -lpthread
