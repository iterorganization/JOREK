#!/bin/sh

export LANG=C
export ARCH=frioul
export JOREK_HOST=frioul_knl
export compilethreads=32
export FORT_BUFFERED=true
export I_MPI_DEBUG=5; 
export KMP_AFFINITY=verbose,scatter

export MKL_DOMAIN_NUM_THREADS="MKL_BLAS=1"
export MKL_NUM_THREADS="1"
export MKL_DYNAMIC="FALSE"
export MKL_ENABLE_INSTRUCTIONS=AVX512_MIC

export MPIRUN="srun -n"
export BATCHCOMMAND="sbatch"
