CC=mpicc
CFLAGS=-Wall
CFLAGS+=-march=native -O5 -DNDEBUG
CFLAGS+=-std=c99
#LFLAGS=-L/home/fweik/Base/lib -lgsl -lgslcblas -lfftw3 
LFLAGS=-L/scratch/fweik/Base/lib -lgsl -lgslcblas -lfftw3 
#Uncomment to add long double 
#LFLAGS+=-lfftw3l
LFLAGS+=-lm

CUDA_COMPILER=nvcc
CUDA_COMPILER_FLAGS=-arch=sm_30 -g -G
CUDA_COMPILER_LFLAGS=-lcufft

OBJECTS=sort.o generate_system.o visit_writer.o window-functions.o  charge-assign.o common.o error.o ewald.o interpol.o io.o p3m-common.o p3m-ik.o realpart.o timings.o p3m-ik-i.o p3m-ad.o p3m-ad-i.o p3m-ad-self-forces.o domain-decomposition.o statistics.o tuning.o p3m-ik-real.o parameters.o p3m-ad-real.o q_ik.o q_ad.o q_ik_i.o q_ad_i.o find_error.o q.o

all: p3mstandalone

shift_test: $(OBJECTS) shift-test.c Makefile
	$(CC) $(CFLAGS) -o shift_test shift-test.c $(OBJECTS) $(LFLAGS)

profile_charge_assignment: $(OBJECTS) profiling/prof_charge_assignment.c
	$(CC) $(CFLAGS) -o prof_ca profiling/prof_charge_assignment.c $(OBJECTS) $(LFLAGS)

tuning4cuda: $(OBJECTS) Makefile tuning4cuda.ca
	$(CC) $(CFLAGS) -o tuning4cuda tuning4cuda.c $(OBJECTS) $(LFLAGS)

test_tuning: $(OBJECTS) Makefile tuning_test.c
	$(CC) $(CFLAGS) -o test_tuning tuning_test.c $(OBJECTS) $(LFLAGS)

tune_ik: $(OBJECTS) Makefile tune_ik
	$(CC) $(CFLAGS) -o tune_ik tune_ik.c $(OBJECTS) $(LFLAGS)

time_assignment: $(OBJECTS) Makefile time_assignment.c
	$(CC) $(CFLAGS) -o time_assignment time_assignment.c $(OBJECTS) $(LFLAGS)

test_dd: $(OBJECTS) Makefile test-dd.c domain-decomposition.o
	$(CC) $(CFLAGS) -o test-dd test-dd.c domain-decomposition.o $(OBJECTS) $(LFLAGS)

p3mstandalone: $(OBJECTS) Makefile main.c
	$(CC) $(CFLAGS) -o p3m main.c $(OBJECTS) $(LFLAGS)

test_vtf: $(OBJECTS) Makefile test_vtf.c
	$(CC) $(CFLAGS) -o test_vtf test_vtf.c $(OBJECTS) $(LFLAGS)

dipolar_system: $(OBJECTS) Makefile dipolar_system.c
	$(CC) $(CFLAGS) -o dipolar_system dipolar_system.c $(OBJECTS) $(LFLAGS)

test: $(OBJECTS) Makefile test.c
	$(CC) $(CFLAGS) -o test test.c $(OBJECTS) $(LFLAGS)

time_method: $(OBJECTS) Makefile time_method.c
	$(CC) $(CFLAGS) -o time_method time_method.c $(OBJECTS) $(LFLAGS)

makefile.dep : *.[ch] Makefile
	for i in *.[c]; do $(CC) -MM $(CFLAGS) "$${i}"; done > $@

include makefile.dep

visit_writer.o:
	gcc -I./tools -c tools/visit_writer.c

p3m_ik_cuda_i.o:
	$(CUDA_COMPILER) $(CUDA_COMPILER_FLAGS) -c -o p3m_ik_cuda_i.o cuda/p3m-ik.cu $(CUDA_COMPILER_LFLAGS)

clean:
	rm -rf *.o p3m

.PHONY = all
