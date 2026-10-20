---
title: "TGCC-CEA"
nav_order: 15
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---

# Running JOREK on TGCC (CEA)

## Getting Access

1. Currently limited to CEA staff, under a single project (PI: Marina Becoulet).
2. French research teams can ask for resources through GENCI via DARI calls: <http://www.edari.fr>
3. Scientists and researchers from academia and industry can ask for resources through PRACE.

## README

- TGCC documentation is available here: <https://www-hpc.cea.fr/tgcc-public/en/html/tgcc-public.html>
- Information on the cluster: <https://www-hpc.cea.fr/en/TGCC.html>

## Login

```bash
ssh -Y <user>@irene-fr.ccc.cea.fr
```

## Software Available

As usual, `module avail`, `module load`, `module list`, and `module show` allow access to the available software.

## Compiling JOREK

Below are some examples of the `.bashrc` and `Makefile.inc` setup.

- Add the following lines to your `.bashrc` file:

```bash
module purge

module load intel/20.0.4
module load mkl/20.0.4
#module load mpi/openmpi
module load mpi/intelmpi
#module load mpi/intelmpi/23.2.0
module load scotch/6.1.0
module load metis
module load flavor/parmetis/r64 parmetis/4.0.3
module load parmetis
module load fftw3/mkl/20.0.4
module load flavor/hdf5/parallel
module load hdf5
module load hwloc/2.5.0
module load nedit
```

- Example `Makefile.inc`:

```makefile
# Makefile.inc file for compiling JOREK on the TGGC supercomputer from CEA

# --- Select physics model
MODEL             = model600

# --- Compiler and options
#FC                = mpif90
#CC                = mpicc
#CXX               = mpicxx -std=c++11

FC                = mpiifort -mt_mpi
CC                = mpiicc   -mt_mpi
CXX               = mpiicpc  -mt_mpi -std=c++11

FFLAGS            = -r8 -i4 -O2 -mavx2 -mp1
FFLAGS           += -DUSE_BLOCK -DUSE_FFTW  #-DWORLDWAR2
FFLAGS           += -fpp -align
#FFLAGS           += -axCOMMON-AVX512
FFLAGS           += -mcmodel=medium
#FFLAGS           += -xHost -O2
FFLAGS           += -qopenmp

#FFLAGS           += -DAUG_WALL

# --- DEBUGGING
#DEBUGFLAGS  = -g
#DEBUGFLAGS += -check bounds -traceback -fpe-all=0 -ftrapuv -check uninit
#DEBUGFLAGS += -traceback -fpe-all=0
#DEBUGFLAGS += -gen-interfaces -warn interfaces
#DEBUGFLAGS += -warn all,nounused -check all,noarg_temp_created -debug all -debug-parameters -fstack-security-check -ftrapuv
#DEBUGFLAGS += -warn all -check all
#DEBUGFLAGS += -check all -check noarg_temp_created -traceback -fpe0
#DEBUGFLAGS += -check stack

# --- Various switches
USE_PASTIX        = 0
USE_PASTIX_MURGE  = 0
USE_HDF5          = 1
USE_FFTW          = 1
USE_MUMPS         = 1
USE_STRUMPACK     = 1
USE_WSMP          = 0
USE_DIRECT_CONSTRUCTION = 0


LIBDIR = /ccc/work/cont003/gen15002/filalexa/lib_rome

# --- BLAS/LAPACK/SCALAPACK
LIBLAPACK = -L$(MKL_LIBDIR) -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
SCALAPACK = -L$(MKL_LIBDIR) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

# --- HWLOC library
HWLOC_HOME  = $(HWLOC_ROOT)/hwloc/install
HWLOC_LIB   = -L$(HWLOC_ROOT)/lib -lhwloc
HWLOC_INC   = -I$(HWLOC_ROOT)/include

# --- METIS, PARMETIS
METIS_HOME = $(METIS_ROOT)
METIS_INC  = $(METIS_ROOT)/include
METIS_LIB  = $(METIS_ROOT)/lib
PARMETIS_HOME = $(PARMETIS_ROOT)
PARMETIS_INC  = $(PARMETIS_ROOT)/include
PARMETIS_LIB  = $(PARMETIS_ROOT)/lib

# --- Scotch 6 library
SCOTCH_HOME  = $(SCOTCH_ROOT)
LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lptscotch -lscotch -lptscotcherr
INC_SCOTCH   = -I$(SCOTCH_HOME)/include

# --- PaStiX solver
PASTIX_HOME       = $(LIBDIR)/pastix_5.2.3/install
#PASTIX_HOME       = $(LIBDIR)/pastix/install
LIB_PASTIX        = -L$(PASTIX_HOME) -lpastix $(LIB_SCOTCH) $(HWLOC_LIB)
LIB_PASTIX_MURGE  = -L$(PASTIX_HOME) -lpastix_murge $(LIB_SCOTCH) $(HWLOC_LIB)
#LIB_PASTIX_BLAS   =
INC_PASTIX        = -I$(PASTIX_HOME)

#  --- STRUMPACK library
STRUMPACK_HOME = $(LIBDIR)/STRUMPACK-7.2.0/install
STRUMPACKINC = -I$(STRUMPACK_HOME)/include -I$(MKL_INC) -I$(METIS_INC) $(INC_SCOTCH) -I$(PARMETIS_INC)  -I$(MKL_INC)
STRUMPACKLIB = -L$(STRUMPACK_HOME)/lib64 -lstrumpack -L$(PARMETIS_LIB) -lparmetis -L$(METIS_LIB) -lmetis $(LIBLAPACK) $(SCALAPACK)

# --- FFTW Library
#LIBFFTW           = $(FFTW_LIB)/libfftw3.a
#INC_FFTW          = -I$(FFTW_INC)

# --- Usually not to be changed
FFLAGS_FIXEDFORM := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NOBOUNDS  := $(FFLAGS)               $(FFLAGS_OMP)
FFLAGS           := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NO_OMP     = $(FFLAGS) $(DEBUGFLAGS)

HDF5INCLUDE = $(HDF5_ROOT)/include
#HDF5LIB     = -L$(HDF5_ROOT)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -L$(SZIP_LIB) -lsz -lz
HDF5LIB     = -L$(HDF5_ROOT)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz

# --- MUMPS solver
MUMPS_HOME = $(LIBDIR)/MUMPS_5.3.5
MUMPS_LIB = $(MUMPS_HOME)/lib
MUMPS_INC = $(MUMPS_HOME)/include
LIB_MUMPS = -L$(MUMPS_LIB) -ldmumps -lzmumps -lmumps_common  -lpord -L$(PARMETIS_LIB) -L$(METIS_LIB) -lparmetis -lmetis  $(SCALAPACK)
INC_MUMPS = $(MUMPS_INC)

# --- Dierckx
LIBDIERCKX = $(LIBDIR)/Dierckx/dierckx.a

LIBS += $(STRUMPACK_LIB) $(LIB_SCOTCH) -lm -ldl -limf -lsvml -lrt -lirng -lc -lmpicxx -lgcc_s -lmpi -lintlc
INCLUDES += -I./tools
DEFINES += -DUSE_STRUMPACK
EXTRA_FLAGS += -std=c++11 -cxxlib
```

## File Systems

The following storage areas are available: `HOME`, `WORK`, and `SCRATCH`.

- To reach your `WORK` area, go to:

```text
/ccc/work/cont003/gen15002/
```

and create your own folder.

- To reach your `SCRATCH` area, go to:

```text
/ccc/scratch/cont003/gen15002/
```

and create your own folder.

- To store simulation files in compressed form, go to:

```text
/ccc/store/cont003/gen15002/
```

and create your own folder.

## Batch System

### Basic commands

```bash
sbatch <jobscript>                              # Submit a job script
sbatch -d afterok:<jobid> <jobscript>           # Submit a job script to start after <jobid> is finished

squeue -t running                               # List all currently running jobs
squeue -u <username>                            # List all own jobs

scancel <jobid>                                 # Delete a queued or running batch job
```

### Example job script

```bash
#!/bin/bash
#MSUB -r JOREK #name
#MSUB -A gen15002 #project
#MSUB -o output_%J.o
#MSUB -e output_%J.e
#MSUB -N 2 #number of nodes
#MSUB -n 8 #number of processes, 4*N
#MSUB -c 32 #number of cores per process to use
#MSUB -T 20000 #in seconds
#MSUB -q rome
#MSUB -Q normal
##MSUB -m work

set -x

ulimit -s unlimited
export OMP_STACKSIZE=256m

export OMP_NUM_THREADS=32
export I_MPI_PIN_MODE=lib

echo "We are in the directory :"
pwd
echo ""
echo "We are running on host :"
hostname
echo ""
echo "Here are some useful environment variables"
echo "------------------------------------------"
echo "To get the current job ID :               BRIDGE_MSUB_JOBID=${BRIDGE_MSUB_JOBID}"
echo "To get the time limit in secondes :       BRIDGE_MSUB_MAXTIME=${BRIDGE_MSUB_MAXTIME}"
echo "To get the number of processes used :     BRIDGE_MSUB_NPROC=${BRIDGE_MSUB_NPROC}"
echo "To get the number of nodes used :         BRIDGE_MSUB_NNODE=${BRIDGE_MSUB_NNODE}"

ccc_mprun ./jorek_model600 < in_n0 > logrun
```

## Accounting

The command `ccc_myproject` gives information about the accounting of your project(s).
