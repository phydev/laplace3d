#!/bin/bash
mpifort -O3 -r8 random_m.F90 laplace_m.F90 main.F90 -o run_lapl
