
mpicc -Wall -g -O3 -lfftw3 -lm -o tuning_detail main_tuning_detail.c generate_system.c common.c realpart.c tuning.c p3m-common.c timings.c interpol.c ewald.c p3m-ik.c charge-assign.c error.c p3m-ad-i.c p3m-ad.c p3m-ik-i.c

