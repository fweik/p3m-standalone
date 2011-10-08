CC=mpicc
CFLAGS=-Wall
LFLAGS=-lfftw3 -lm

all: p3mstandalone

p3mstandalone: charge-assign.o common.o error.o ewald.o interpol.o io.o main.o p3m-common.o p3m-ik.o realpart.o timings.o Makefile main.c main.h
	$(CC) $(LFLAGS) charge-assign.o common.o error.o ewald.o interpol.o io.o main.o p3m-common.o p3m-ik.o realpart.o main.c -o p3m

charge-assign.o: charge-assign.c charge-assign.h
	$(CC) -c $(CFLAGS) charge-assign.c
common.o: common.c common.h
	$(CC) -c $(CFLAGS) common.c
error.o: error.c error.h
	$(CC) -c $(CFLAGS) error.c
ewald.o: ewald.c ewald.h
	$(CC) -c $(CFLAGS) ewald.c
interpol.o: interpol.c interpol.h
	$(CC) -c $(CFLAGS) interpol.c
io.o: io.c io.h
	$(CC) -c $(CFLAGS) io.c
p3m-common.o: p3m-common.c p3m-common.h
	$(CC) -c $(CFLAGS) p3m-common.c
p3m-ik.o: p3m-ik.c p3m-ik.h
	$(CC) -c $(CFLAGS) p3m-ik.c
realpart.o: realpart.c realpart.h
	$(CC) -c $(CFLAGS) realpart.c

.PHONY = all
