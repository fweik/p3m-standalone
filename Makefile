CC=mpicc
CFLAGS=-Wall -g
LFLAGS=-lgsl -lgslcblas -lfftw3 -lm

all: p3mstandalone

p3mstandalone: window-functions.o greens.o charge-assign.o common.o error.o ewald.o interpol.o io.o p3m-common.o p3m-ik.o realpart.o timings.o p3m-ik-i.o p3m-ad.o p3m-ad-i.o  Makefile main.c
	$(CC) $(CFLAGS) $(LFLAGS) window-functions.o greens.o charge-assign.o common.o error.o ewald.o interpol.o io.o p3m-common.o p3m-ik.o  p3m-ik-i.o p3m-ad.o p3m-ad-i.o realpart.o main.c -o p3m

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o p3m

.PHONY = all
