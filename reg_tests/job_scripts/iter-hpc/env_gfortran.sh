#!/bin/bash
eval `tclsh /work/imas/opt/modules-tcl/modulecmd.tcl $(basename $SHELL) autoinit`

module purge

module use /work/imas/opt/EasyBuild/modules/all
module use /work/imas/etc/modules/all

# There is an issue with PaStiX/5.2.3, use an older version
module load MUMPS/5.1.2-foss-2018a-metis
module load PaStiX/5.2.3-foss-2018a
module load HDF5/1.10.1-foss-2018a
#export PASTIX_HOME=/work/imas/opt/EasyBuild/software/PaStiX/5.2.2.22-goolf-1.5.16
#export SCOTCH_HOME=/work/imas/opt/EasyBuild/software/SCOTCH/6.0.4-goolf-1.5.16


export LANG=C
export JOREK_HOST=iter-hpc
export compilethreads=4
export MAKEFLAGS="-j$compilethreads"
export PRERUN="export OMP_NUM_THREADS=4"
export MPIRUN="mpirun -np "
export BATCHCOMMAND="qsub"

export http_proxy=${JOREK_HTTP_PROXY}
export https_proxy=${JOREK_HTTP_PROXY}

export MPIRUN='mpirun --allow-run-as-root -n'
