#!/bin/bash

export LANG=C
export ARCH=irene_skl
export JOREK_HOST=irene_skl
export compilethreads=32

USE_INTELMPI=$(mpirun -V | grep -c -i "Intel.*MPI")
USE_OPENMPI=$(mpirun -V | grep -c -i "Open.*MPI")
MPIRUN="srun -n "
export MKL_DOMAIN_NUM_THREADS="MKL_BLAS=1"
export MKL_NUM_THREADS="1"
export MKL_DYNAMIC="FALSE"
export PRERUN=""

if [ -n "$OMP_NUM_THREADS" ]; then
  ((PPN=48/$OMP_NUM_THREADS))
else
  PPN=1
fi
export I_MPI_DEBUG=4
export KMP_AFFINITY=verbose,scatter
export MPIRUN
export BATCHCOMMAND="ccc_msub"
