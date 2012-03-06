CC=mpicc
CFLAGS=-Wall -g
LFLAGS=-lgsl -lgslcblas -lfftw3 -lm

all: p3mstandalone

p3mstandalone: generate_system.o visit_writer.o window-functions.o greens.o charge-assign.o common.o error.o ewald.o interpol.o io.o p3m-common.o p3m-ik.o realpart.o timings.o p3m-ik-i.o p3m-ad.o p3m-ad-i.o  Makefile main.c
	$(CC) $(CFLAGS) $(LFLAGS) generate_system.o visit_writer.o window-functions.o greens.o charge-assign.o common.o error.o ewald.o interpol.o io.o p3m-common.o p3m-ik.o  p3m-ik-i.o p3m-ad.o p3m-ad-i.o realpart.o main.c -o p3m

random_system: window-functions.o greens.o charge-assign.o common.o error.o ewald.o interpol.o io.o p3m-common.o p3m-ik.o realpart.o timings.o p3m-ik-i.o p3m-ad.o p3m-ad-i.o generate_system.o  Makefile main_random.c
	$(CC) $(CFLAGS) $(LFLAGS) generate_system.o window-functions.o greens.o charge-assign.o common.o error.o ewald.o interpol.o io.o p3m-common.o p3m-ik.o  p3m-ik-i.o p3m-ad.o p3m-ad-i.o realpart.o main_random.c -o random_system

madelung: window-functions.o greens.o charge-assign.o common.o error.o ewald.o interpol.o p3m-common.o p3m-ik.o realpart.o timings.o p3m-ik-i.o p3m-ad.o p3m-ad-i.o generate_system.o  Makefile main_madelung.c
	$(CC) $(CFLAGS) $(LFLAGS) generate_system.o window-functions.o greens.o charge-assign.o common.o error.o ewald.o interpol.o p3m-common.o p3m-ik.o  p3m-ik-i.o p3m-ad.o p3m-ad-i.o realpart.o main_madelung.c -o madelung

estimates: window-functions.o greens.o charge-assign.o common.o error.o ewald.o interpol.o io.o p3m-common.o p3m-ik.o realpart.o timings.o p3m-ik-i.o p3m-ad.o p3m-ad-i.o generate_system.o  Makefile main_estimates.c
	$(CC) $(CFLAGS) $(LFLAGS) generate_system.o window-functions.o greens.o charge-assign.o common.o error.o ewald.o interpol.o io.o p3m-common.o p3m-ik.o  p3m-ik-i.o p3m-ad.o p3m-ad-i.o realpart.o main_estimates.c -o estimates
makefile.dep : *.[ch]
	for i in *.[c]; do $(CC) -MM $(CFLAGS) "$${i}"; done > $@

include makefile.dep

visit_writer.o:
	gcc -I./tools -c tools/visit_writer.c

clean:
	rm -rf *.o p3m

.PHONY = all
