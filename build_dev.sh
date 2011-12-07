#!/bin/bash

set -e

mpicc -g -lm -lfftw3 main.c ewald.c error.c realpart.c common.c io.c p3m-common.c 
charge-assign.c p3m-ik.c timings.c p3m-ik-i.c greens.c p3m-ad.c p3m-ad-i.c





