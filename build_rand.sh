#!/bin/bash

set -e

mpicc -O3 -Wall -DNDEBUG -o random_system -lfftw3 -lm main_random.c ewald.c error.c realpart.c common.c io.c p3m-common.c charge-assign.c p3m-ik.c timings.c p3m-ik-i.c greens.c p3m-ad.c generate_system.c p3m-ad-i.c





