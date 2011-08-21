CC=gcc
CFLAGS=-Wall -fopenmp -g -lfftw3 -lm

all: p3mstandalone

p3mstandalone: charge-assign.o common.o dipol.o error.o ewald.o interpol.o io.o main.o p3m-common.o p3m-error.o p3m-ik.o realpart.o
	$(CC) charge-assign.o common.o dipol.o error.o ewald.o interpol.o io.o main.o p3m-common.o p3m-error.o p3m-ik.o realpart.o -o p3m

charge-assign.o: charge-assign.c
	$(CC) $(CFLAGS) charge-assign.c
common.o: common.c
	$(CC) $(CFLAGS) common.c
dipol.o: dipol.c
	$(CC) $(CFLAGS) dipol.c
error.o: error.c
	$(CC) $(CFLAGS) error.c
ewald.o: ewald.c
	$(CC) $(CFLAGS) ewald.c
interpol.o: interpol.c
	$(CC) $(CFLAGS) interpol.c
io.o: io.c
	$(CC) $(CFLAGS) io.c
main.o: main.cs 
	$(CC) $(CFLAGS) main.c
p3m-common.o: p3m-common.c
	$(CC) $(CFLAGS) p3m-common.c
p3m-error.o: p3m-error.c
	$(CC) $(CFLAGS) p3m-error.c
p3m-ik.o: p3m-ik.c
	$(CC) $(CFLAGS) p3m-ik.c
realpart.o: realpart.c
	$(CC) $(CFLAGS) realpart.c