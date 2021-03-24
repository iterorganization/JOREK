#!/bin/bash

#module purge
#module load intel/12.1.0 mvapich2/1.8.1 mkl/10.3

export LANG=C
export JOREK_HOST=poincare
export compilethreads=8
export BATCHCOMMAND="llsubmit"
#export KMP_STACK_SIZE=16M
export KMP_AFFINITY="verbose,norespect,compact" # if this is not set, the OpenMP threads are confined on one core


export MKL_NUM_THREADS=1
export MKL_DYNAMIC=0
export MV2_ENABLE_AFFINITY=0
export MPIRUN="mpiexec -envlist LANG,OMP_NUM_THREADS,KMP_STACK_SIZE,KMP_AFFINITY,MKL_NUM_THREADS,MKL_DYNAMIC,MV2_ENABLE_AFFINITY -launcher-exec /opt/ibmll/LoadL/scheduler/full/bin/llspawn.stdio -ppn 1 -binding -verbose -f $LOADL_HOSTFILE -n"

