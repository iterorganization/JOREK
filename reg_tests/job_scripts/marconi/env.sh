#!/bin/bash

export LANG=C
export ARCH=marconi
export JOREK_HOST=marconi
export compilethreads=32
USE_INTELMPI=$(mpirun -V | grep -c -i "Intel.*MPI")
USE_OPENMPI=$(mpirun -V | grep -c -i "Open.*MPI")
MPIRUN="mpirun -np"
export MKL_DOMAIN_NUM_THREADS="MKL_BLAS=1"
export MKL_NUM_THREADS="1"
export MKL_DYNAMIC="FALSE"
export PRERUN="uniq $PBS_NODEFILE > hosts.txt"
#-genv I_MPI_PIN=1 -genv I_MPI_PIN_DOMAIN=socket -genv I_MPI_PIN_ORDER=spread -genv OMP_PROC_BIND=1
if [ -n "$OMP_NUM_THREADS" ]; then
  ((PPN=36/$OMP_NUM_THREADS))
else
  PPN=1
fi
export I_MPI_DEBUG=4
export KMP_AFFINITY=verbose,scatter
if [ ${USE_INTELMPI} -gt 0 ]; then
  MPIRUN="mpiexec.hydra -genvall -ppn $PPN -f hosts.txt -np"
fi
if [ ${USE_OPENMPI} -gt 0 ]; then
  MPIRUN="mpirun -x KMP_AFFINITY=verbose -x OMP_PROC_BIND=1 -n"
fi
export MPIRUN
export BATCHCOMMAND="qsub"
