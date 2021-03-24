#!/bin/bash
eval `tclsh /work/imas/opt/modules-tcl/modulecmd.tcl $(basename $SHELL) autoinit`

module purge

module use /work/imas/etc/modules/all

module load cURL/7.58.0-GCCcore-6.4.0
module load MUMPS/5.1.2-intel-2018a-metis
#module load PaStiX/5.2.3-intel-2018a
module load FFTW/3.3.7-intel-2018a
module load HDF5/1.10.1-intel-2018a
module load STRUMPACK/3.3.0-intel-2018a
module load ParMETIS/4.0.3-intel-2018a

module avail STRUMPACK
module show STRUMPACK/3.3.0-intel-2018a
module list
ls $EBROOTSTRUMPACK/lib

export METIS_HOME=$EBROOTPARMETIS
export STRUMPACK_HOME=$EBROOTSTRUMPACK
export LANG=C
export JOREK_HOST=iter-hpc
export compilethreads=4
export MAKEFLAGS="-j$compilethreads"
export PRERUN="export OMP_NUM_THREADS=4"
export MPIRUN="mpirun -np "
export BATCHCOMMAND="qsub"
export CXXFLAGS=-O0 # problem with stdio library on ITER http://gcc.1065356.n8.nabble.com/g-4-8-fails-with-Ox-option-td953876.html

export http_proxy=${JOREK_HTTP_PROXY}
export https_proxy=${JOREK_HTTP_PROXY}
