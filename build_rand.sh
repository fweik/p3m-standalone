#!/bin/bash

set -e

mpicc -O3 -Wall -DNDEBUG -lm -lfftw3 -o random_system main_random.c ewald.c error.c realpart.c common.c io.c p3m-common.c charge-assign.c p3m-ik.c timings.c p3m-ik-i.c greens.c p3m-ad.c generate-system.c p3m-ad-i.c





