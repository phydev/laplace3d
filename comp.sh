#!/bin/bash
mpifort -O3 -r8 src/random_m.F90 src/routines.F90 src/laplace_m.F90 src/main.F90 -o run_lapl
