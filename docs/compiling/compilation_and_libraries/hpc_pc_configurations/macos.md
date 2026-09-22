---
title: "macOS"
nav_order: 16
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---

# Building JOREK on Apple computers with macOS using GCC

Here is a short description of how to build JOREK on Apple computers with macOS using GCC. This tutorial is based on [Homebrew](https://brew.sh/), but should work similarly with [MacPorts](https://www.macports.org/). The description is not optimized for performance, but just to get the code up and running.

On macOS there are some differences compared to Linux:

- the default compiler is Clang, not GCC
- Unix tools such as `readlink` are not compatible
- OpenMPI does not seem to work with GCC 7.3, so MPICH is used instead

## Required steps

- Install tools, compiler, MPICH, and HDF5 with MPI and Fortran support. It is important here to compile MPICH with `gcc-7`.

```bash
brew install coreutils findutils gcc
brew install --cc=gcc-7 mpich
brew install hdf5 --with-fortran --with-mpi
```

- Set the path to use GNU tools instead of the macOS ones.

```bash
export PATH="/usr/local/opt/coreutils/libexec/gnubin:$PATH"
```

- Compile SCOTCH ([see the main compilation guide](/JOREK/compiling/compilation_and_libraries/compiling.html)) with:

```makefile
EXE=
LIB= .a
OBJ= .o

MAKE= make
AR= ar
ARFLAGS= -ruv
CAT= cat
CCS= gcc-7
CCP= mpicc
CCD= gcc-7
CFLAGS= -O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_RENAME_PARSER -DSCOTCH_PTHREAD -DIDXSIZE64 -DINTSIZE32
CLIBFLAGS=
LDFLAGS= -g -lz -lm
CP= cp
LEX= flex -Pscotchyy -olex.yy.c
LN= ln
MKDIR= mkdir
MV= mv
RANLIB= ranlib
YACC= bison -pscotchyy -y -b y
```

- Compile PaStiX ([see the main compilation guide](/JOREK/compiling/compilation_and_libraries/compiling.html)) with:

```makefile
HOSTARCH    = i686_mac
VERSIONBIT  = _64bit
EXEEXT      =
OBJEXT      = .o
LIBEXT      = .a
CCPROG      = gcc-7 -Wall
CFPROG      = gfortran
CF90PROG    = gfortran
MCFPROG     = mpif90
CF90CCPOPT  = -fpp
# Compilation options for optimization (make expor)
CCFOPT      = -O3
# Compilation options for debug (make | make debug)
CCFDEB      = -g3

LKFOPT      =
MKPROG      = make
MPCCPROG    = mpicc -Wall
CPP         = cpp-7
ARFLAGS     = ruv
ARPROG      = ar
# Uncomment the correct line
EXTRALIB    = -lifcore -lm


VERSIONMPI  = _mpi
VERSIONSMP  = _smp
VERSIONSCH  = _static
VERSIONINT  = _int
VERSIONPRC  = _simple
VERSIONFLT  = _real
VERSIONORD  = _scotch

###################################################################
#                  SETTING INSTALL DIRECTORIES                    #
###################################################################
# ROOT       = /path/to/install/directory
# INCLUDEDIR = ${ROOT}/include
# LIBDIR     = ${ROOT}/lib
# BINDIR     = ${ROOT}/bin

###################################################################
#                          INTEGER TYPE                           #
###################################################################
# Uncomment the following lines for integer type support (Only 1)

#VERSIONINT  = _long
#CCTYPES     = -DFORCE_LONG -DLONG
#---------------------------
#VERSIONINT  = _int32
CCTYPES     = -DFORCE_INT32 -DINTSIZE32
#---------------------------
#VERSIONINT  = _int64
#CCTYPES     = -DFORCE_INT64 -DINTSSIZE64

###################################################################
#                           FLOAT TYPE                            #
###################################################################
CCTYPESFLT  =
# Uncomment the following lines for double precision support
VERSIONPRC  = _double
CCTYPESFLT := $(CCTYPESFLT) -DFORCE_DOUBLE -DPREC_DOUBLE

# Uncomment the following lines for float=complex support
#VERSIONFLT  = _complex
#CCTYPESFLT := $(CCTYPESFLT) -DFORCE_COMPLEX -DTYPE_COMPLEX


###################################################################
#                          MPI/THREADS                            #
###################################################################

# Uncomment the following lines for sequential (NOMPI) version
#VERSIONMPI  = _nompi
#CCTYPES    := $(CCTYPES) -DFORCE_NOMPI
#MPCCPROG    = $(CCPROG)
#MCFPROG     = $(CFPROG)

# Uncomment the following lines for non-threaded (NOSMP) version
#VERSIONSMP  = _nosmp
#CCTYPES    := $(CCTYPES) -DFORCE_NOSMP

# Uncomment the following line to enable a progression thread
#CCPASTIX   := $(CCPASTIX) -DTHREAD_COMM

# Uncomment the following line if your MPI doesn't support MPI_THREAD_MULTIPLE level
CCPASTIX   := $(CCPASTIX) -DPASTIX_FUNNELED

# Uncomment the following line if your MPI doesn't support MPI_Datatype correctly
CCPASTIX   := $(CCPASTIX) -DNO_MPI_TYPE

# Uncomment the following line if you want to use semaphore barrier
# instead of MPI barrier (with IPARM_AUTOSPLIT_COMM)
#CCPASTIX    := $(CCPASTIX) -DWITH_SEM_BARRIER

# Uncomment the following lines to enable StarPU.
#CCPASTIX   := ${CCPASTIX} `pkg-config libstarpu --cflags` -DWITH_STARPU
#EXTRALIB   := $(EXTRALIB) `pkg-config libstarpu --libs`

# Uncomment the following line to enable StarPU profiling
# ( IPARM_VERBOSE > API_VERBOSE_NO ).
#CCPASTIX   := ${CCPASTIX} -DSTARPU_PROFILING

# Uncomment the following line to disable CUDA (StarPU)
#CCPASTIX   := ${CCPASTIX} -DFORCE_NO_CUDA

###################################################################
#                          Options                                #
###################################################################

# Show memory usage statistics
CCPASTIX   := $(CCPASTIX) -DMEMORY_USAGE

# Show memory usage statistics in solver
CCPASTIX   := $(CCPASTIX) -DSTATS_SOPALIN

# Uncomment following line for dynamic thread scheduling support
#CCPASTIX   := $(CCPASTIX) -DPASTIX_DYNSCHED

# Uncomment the following lines for Out-of-core
#CCPASTIX   := $(CCPASTIX) -DOOC -DOOC_NOCOEFINIT -DOOC_DETECT_DEADLOCKS

###################################################################
#                      GRAPH PARTITIONING                         #
###################################################################

# Uncomment the following lines for using metis ordering
#VERSIONORD  = _metis
#METIS_HOME  = ${HOME}/metis-4.0
#CCPASTIX   := $(CCPASTIX) -DMETIS -I$(METIS_HOME)/Lib
#EXTRALIB   := $(EXTRALIB) -L$(METIS_HOME) -lmetis

# Scotch always needed to compile
SCOTCH_HOME ?= $(HOME)/src/scotch_5.1.12
SCOTCH_INC ?= $(SCOTCH_HOME)/include
SCOTCH_LIB ?= $(SCOTCH_HOME)/lib
# Uncomment on of this blocks
#scotch
#CCPASTIX   := $(CCPASTIX) -I$(SCOTCH_INC) -DWITH_SCOTCH
#EXTRALIB   := $(EXTRALIB) -L$(SCOTCH_LIB) -lscotch -lscotcherrexit
#ptscotch
CCPASTIX   := $(CCPASTIX) -I$(SCOTCH_INC) -DDISTRIBUTED -DWITH_SCOTCH
EXTRALIB   := $(EXTRALIB) -L$(SCOTCH_LIB) -lscotch -lptscotch -lscotcherrexit

###################################################################
#                Portable Hardware Locality                       #
###################################################################
# If HwLoc library is available, uncomment the following lines to bind threads correctly to CPUs
#HWLOC_HOME ?= /opt/hwloc/
#HWLOC_INC  ?= $(HWLOC_HOME)/include
#HWLOC_LIB  ?= $(HWLOC_HOME)/lib
#CCPASTIX   := $(CCPASTIX) -I$(HWLOC_INC) -DWITH_HWLOC
#EXTRALIB   := $(EXTRALIB) -L$(HWLOC_LIB) -lhwloc

###################################################################
#                             MARCEL                              #
###################################################################

# Uncomment following lines for marcel thread support
#VERSIONSMP := $(VERSIONSMP)_marcel
#CCPASTIX   := $(CCPASTIX) `pm2-config --cflags` -I${PM2_ROOT}/marcel/include/pthread
#EXTRALIB   := $(EXTRALIB) `pm2-config --libs`
# ---- Thread Posix ------
EXTRALIB   := $(EXTRALIB) -lpthread

# Uncomment following line for bubblesched framework support (need marcel support)
#VERSIONSCH  = _dyn
#CCPASTIX   := $(CCPASTIX) -DPASTIX_BUBBLESCHED

###################################################################
#                              BLAS                               #
###################################################################

# Choose BLAS library (only 1)
# Do not forget to set BLAS_HOME if it is not in your environment
#BLAS_HOME=${MKLROOT}/lib/intel64
#----  Blas    ----
BLASLIB =  -lblas
#---- Gotoblas ----
#BLASLIB =  -L$(BLAS_HOME) -lgoto
#----  MKL     ----
# Uncomment the correct line
#BLASLIB =  -L$(BLAS_HOME) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#BLASLIB =  -L$(BLAS_HOME) -lmkl_intel -lmkl_sequential -lmkl_core
#----  Acml    ----
#BLASLIB =  -L$(BLAS_HOME) -lacml

###################################################################
#                         PYTHON WRAPPER                          #
###################################################################
#MPI4PY_DIR    = /path/to/mpi4py
#MPI4PY_INC    = $(MPI4PY_DIR)/src/include/
#MPI4PY_LIBDIR = $(MPI4PY_DIR)/build/lib.linux-x86_64-2.7/
#PYTHON_INC    = /usr/include/python2.7/
#CCTYPES      := $(CCTYPES) -fPIC

###################################################################
#                          DO NOT TOUCH                           #
###################################################################

FOPT      := $(CCFOPT)
FDEB      := $(CCFDEB)
CCHEAD    := $(CCPROG) $(CCTYPES) $(CCFOPT)
CCFOPT    := $(CCFOPT) $(CCTYPES) $(CCPASTIX)
CCFDEB    := $(CCFDEB) $(CCTYPES) $(CCPASTIX)


###################################################################
#                        MURGE COMPATIBILITY                      #
###################################################################
# Uncomment if y
```

- Compile JOREK ([see the main compilation guide](/JOREK/compiling/compilation_and_libraries/compiling.html)) using:

```bash
CC=gcc-7 make -j8
```

See more below. Use the following `Makefile.inc`:

```makefile
# Configuration file for jorek
# model directory
MODEL = model199

# Fortran compiler
FC     = mpif90
CC     = mpicc
CXX    = mpicxx

FFLAGS_OMP       = -fopenmp
FFLAGS           = -g -O0
FFLAGS          := $(FFLAGS) -fdefault-real-8 -fdefault-double-8
FFLAGS          := $(FFLAGS) -ffree-line-length-none
FFLAGS          := $(FFLAGS) -x f95-cpp-input
## TODO find equivalents of this ifort options
#FFLAGS          := $(FFLAGS) -traceback -xT -align

FFLAGS_FIXEDFORM = $(FFLAGS) -ffixed-form  -fbounds-check $(FFLAGS_OMP)
FFLAGS_NOBOUNDS  = $(FFLAGS) -ffixed-form                 $(FFLAGS_OMP)
FFLAGS          := $(FFLAGS)               -fbounds-check $(FFLAGS_OMP)

#FFLAGS_NO_OMP   = -fpp -O2 -xT -align
FFLAGS          := $(FFLAGS) -DMURGE_INTERFACE_MAJOR_VERSION=1 -DMURGE_INTERFACE_MINOR_VERSION=0
FFLAGS := $(FFLAGS) -cpp  -DUSE_R3_INFO_MPI -DUSE_R3_INFO
CFLAGS :=  -DUSE_R3_INFO_MPI -DUSE_R3_INFO

# Solvers dependencies
USE_HIPS   = 0
USE_PASTIX = 1
USE_BLOCK  = 1
USE_PASTIX_MURGE = 0
USE_MUMPS  = 0
USE_WSMP   = 0
USE_FFTW   = 1
USE_HDF5   = 1

#Scotch library
SCOTCH_HOME  = $(HOME)/src/scotch_5.1.12
LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lscotch -lscotcherr
INC_SCOTCH   = -I$(SCOTCH_HOME)/include

# MUMPS
MUMPS_HOME ?=
LIB_MUMPS  = $(MUMPS_HOME)/lib/libdmumps.a
LIB_ZMUMPS = $(MUMPS_HOME)/lib/libzmumps.a
INC_MUMPS  = $(MUMPS_HOME)/include/
ORDLIB     = $(MUMPS_HOME)/lib/libpord.a

# PASTIX
PASTIX_HOME      ?= $(HOME)/src/pastix_release_4492/install
LIB_PASTIX       = `$(PASTIX_HOME)/pastix-conf --libs`
LIB_PASTIX_MURGE = `$(PASTIX_HOME)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS  = -lblas -llapack
INC_PASTIX       = `$(PASTIX_HOME)/pastix-conf --incs`
PASTIX_MEMORY_USAGE = 0

#HIPS
LIBHIPS   = -L$(HOME)/hips_trunk/trunk/LIB -lhips -lio $(LIBSCOTCH)
INCHIPS   = -I$(HOME)/hips_trunk/trunk/LIB

#WSMP
LIB_WSMP = -L$(WSMP_HOME) -lpwsmp64

# LIBFFTW
FFTWDIR ?= /usr/local
LIBFFTW   = -L$(FFTWDIR)/lib -lfftw3
INC_FFTW  = -I$(FFTWDIR)/include

# LIBDIERCKX
DIERCKX_HOME ?= /home/ITER/vanvugd/source/Dierckx
LIBDIERCKX = -L$(DIERCKX_HOME) -ldierckx

#LAPACK
#LIBLAPACK = -L$(MKL) -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#SCALAP    = -L$(MKL) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

HDF5_HOME=/usr/local
HDF5INCLUDE=$(HDF5_HOME)/include
HDF5LIB=$(HDF5_HOME)/lib/libhdf5hl_fortran.a \
        $(HDF5_HOME)/lib/libhdf5_hl.a        \
        $(HDF5_HOME)/lib/libhdf5_fortran.a   \
        $(HDF5_HOME)/lib/libhdf5.a           \
        -lz -lsz

        #$(EBROOTZLIB)/lib/libz.a
```

## Possible compilation error

When compiling, you may encounter an error like:

```bash
Undefined symbols for architecture x86_64:
  "_r3_info_begin_", referenced from:
      ___construct_matrix_mod_MOD_construct_matrix in construct_matrix_mod.o
      ___construct_matrix_murge_mod_MOD_construct_matrix_murge in construct_matrix_murge_mod.o
      ___solve_mat_n_MOD_solve_matrix_n in solve_mat_n.o
      _gmres_driver_ in gmres_driver.o
  "_r3_info_end_", referenced from:
      ___construct_matrix_mod_MOD_construct_matrix in construct_matrix_mod.o
      ___solve_mat_n_MOD_solve_matrix_n in solve_mat_n.o
      _gmres_driver_ in gmres_driver.o
  "_r3_info_init_", referenced from:
      _MAIN__ in jorek2_main.o
  "_r3_info_print_", referenced from:
      _MAIN__ in jorek2_main.o
  "_r3_info_summary_", referenced from:
      _MAIN__ in jorek2_main.o
ld: symbol(s) not found for architecture x86_64
collect2: error: ld returned 1 exit status
make: *** [jorek2_main] Error 1
```

In this case, add these two lines:

```bash
jorek2_main: .obj/r3_info.o
jorek2_main: .obj/r3_ctlk.o
```

to `.dep/jorek2_main.d`, comment out the line

```c
#include <sys/procfs.h>
```

in `tools/r3_ctlk.c`, and run `make` again.
