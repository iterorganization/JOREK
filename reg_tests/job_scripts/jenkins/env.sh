#!/bin/bash

export LANG=C
export JOREK_HOST=jenkins
export compilethreads=2
export BATCHCOMMAND="bash"
export PRERUN="export OMP_NUM_THREADS=2"
export MPIRUN="mpirun.mpich2 -np "
