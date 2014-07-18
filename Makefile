CC=mpicc
#CFLAGS=-Wall -O5 -DNDEBUG -mavx
CFLAGS=-Wall -O5 -DNDEBUG -march=native
#CFLAGS=-Wall -O0 -march=native -g
CFLAGS+=-std=gnu99
#LFLAGS=-L/home/fweik/Base/lib -lgsl -lgslcblas -lfftw3 
LFLAGS=-L/scratch/fweik/Base/lib -lgsl -lgslcblas -lfftw3 
#Uncomment to add long double 
#LFLAGS+=-lfftw3l
LFLAGS+=-lm

CUDA_COMPILER=nvcc
CUDA_COMPILER_FLAGS=-arch=sm_30 -g -G
CUDA_COMPILER_LFLAGS=-lcufft

OBJECTS=sort.o generate_system.o visit_writer.o window-functions.o  charge-assign.o common.o error.o ewald.o interpol.o io.o p3m-common.o p3m-ik.o realpart.o timings.o p3m-ik-i.o p3m-ad.o p3m-ad-i.o p3m-ad-self-forces.o domain-decomposition.o statistics.o tuning.o p3m-ik-real.o parameters.o p3m-ad-real.o q_ik.o q_ad.o q_ik_i.o q_ad_i.o find_error.o q.o p3m-ik-real-ns.o

BINARIES=prof_ca time_assignment test_tuning p3m tuning_density

all: p3mstandalone

profile_charge_assignment: $(OBJECTS) Makefile profiling/prof_charge_assignment.c
	$(CC) $(CFLAGS) -o prof_ca profiling/prof_charge_assignment.c $(OBJECTS) $(LFLAGS)

make_reference: $(OBJECTS) Makefile make_reference.c
	$(CC) $(CFLAGS) -o make_reference make_reference.c $(OBJECTS) $(LFLAGS)

time_assignment: $(OBJECTS) Makefile profiling/time_assignment.c
	$(CC) $(CFLAGS) -I. -o time_assignment profiling/time_assignment.c $(OBJECTS) $(LFLAGS)

test_tuning: $(OBJECTS) Makefile tuning_test.c
	$(CC) $(CFLAGS) -o test_tuning tuning_test.c $(OBJECTS) $(LFLAGS)

tuning_density: $(OBJECTS) Makefile tuning_density.c
	$(CC) $(CFLAGS) -o tuning_density tuning_density.c $(OBJECTS) $(LFLAGS)

p3mstandalone: $(OBJECTS) Makefile main.c
	$(CC) $(CFLAGS) -o p3m main.c $(OBJECTS) $(LFLAGS)

test_vtf: $(OBJECTS) Makefile test_vtf.c
	$(CC) $(CFLAGS) -o test_vtf test_vtf.c $(OBJECTS) $(LFLAGS)

dipolar_system: $(OBJECTS) Makefile dipolar_system.c
	$(CC) $(CFLAGS) -o dipolar_system dipolar_system.c $(OBJECTS) $(LFLAGS)

test: $(OBJECTS) Makefile test.c
	$(CC) $(CFLAGS) -o test test.c $(OBJECTS) $(LFLAGS)

makefile.dep : *.[ch] Makefile
	for i in *.[c]; do $(CC) -MM $(CFLAGS) "$${i}"; done > $@

include makefile.dep

visit_writer.o:
	gcc -I./tools -c tools/visit_writer.c

clean:
	rm -rf *.o
	rm -rf $(BINARIES)

.PHONY = all
