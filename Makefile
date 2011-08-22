CC=gcc
CFLAGS=-Wall -fopenmp -g -lfftw3 -lm

all: p3mstandalone

p3mstandalone: charge-assign.o common.o dipol.o error.o ewald.o interpol.o io.o main.o p3m-common.o p3m-error.o p3m-ik.o realpart.o Makefile main.c main.h
	$(CC) charge-assign.o common.o dipol.o error.o ewald.o interpol.o io.o main.o p3m-common.o p3m-error.o p3m-ik.o realpart.o main.c -o p3m

charge-assign.o: charge-assign.c charge-assign.h
	$(CC) $(CFLAGS) charge-assign.c
common.o: common.c common.h
	$(CC) $(CFLAGS) common.c
dipol.o: dipol.c dipol.h
	$(CC) $(CFLAGS) dipol.c
error.o: error.c error.h
	$(CC) $(CFLAGS) error.c
ewald.o: ewald.c ewald.h
	$(CC) $(CFLAGS) ewald.c
interpol.o: interpol.c interlpol.h
	$(CC) $(CFLAGS) interpol.c
io.o: io.c io.h
	$(CC) $(CFLAGS) io.c
p3m-common.o: p3m-common.c p3m-common.h
	$(CC) $(CFLAGS) p3m-common.c
p3m-error.o: p3m-error.c p3m-error.h
	$(CC) $(CFLAGS) p3m-error.c
p3m-ik.o: p3m-ik.c p3m-ik.h
	$(CC) $(CFLAGS) p3m-ik.c
realpart.o: realpart.c realpart.h
	$(CC) $(CFLAGS) realpart.c

.PHONY = all