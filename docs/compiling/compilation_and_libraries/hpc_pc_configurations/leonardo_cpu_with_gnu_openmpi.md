---
title: "Leonardo with GNU OpenMPI"
nav_order: 2
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---
# Compiling JOREK on Leonardo-CPU with Intel and Oneapi
  ## bashrc
  ```bash
  module load gcc/12.2.0
  module load openmpi/4.1.6--gcc--12.2.0
  module load openblas/0.3.24--gcc--12.2.0
  module load netlib-scalapack/2.2.0--openmpi--4.1.6--gcc--12.2.0
  module load fftw/3.3.10--openmpi--4.1.6--gcc--12.2.0
  module load boost/1.83.0--openmpi--4.1.6--gcc--12.2.0
  module load hdf5/1.14.3--openmpi--4.1.6--gcc--12.2.0
  module load metis/5.1.0--gcc--12.2.0
  module load parmetis/4.0.3--openmpi--4.1.6--gcc--12.2.0
  # Other module loads
  module load python/3.10.8--gcc--8.5.0
  
  # PUBLIC OpenMPI compiled libraries with LEONARDO METIS/PARMETIS
  export BASE_LIBRARY_PATH="/leonardo/pub/userexternal/csommari"
  export MUMPS_PUBLIC_LEO_ROOT="$BASE_LIBRARY_PATH/jorek_libraries/openmpi416_gcc122/MUMPS_v541_openmpi416_gcc122_leonardo_metis510_parmetis403_13022025"
  export MUMPS_PUBLIC_LEO_LIB="$MUMPS_PUBLIC_LEO_ROOT/lib"
  export MUMPS_PUBLIC_LEO_INC="$MUMPS_PUBLIC_LEO_ROOT/include"
  export STRUMPACK_PUBLIC_LEO_ROOT="$BASE_LIBRARY_PATH/jorek_libraries/openmpi416_gcc122/STRUMPACK_v710_openmpi416_gcc122_leonardo_metis510_parmetis403_13022025"
  export STRUMPACK_PUBLIC_LEO_LIB="$STRUMPACK_PUBLIC_LEO_ROOT/lib/cmake/STRUMPACK"
  export STRUMPACK_PUBLIC_LEO_LIB64="$STRUMPACK_PUBLIC_LEO_ROOT/lib64"
  export STRUMPACK_PUBLIC_LEO_INC="$STRUMPACK_PUBLIC_LEO_ROOT/include"
```

## Makefile.inc 

```bash
# model directory
MODEL = model600

# Fortran compiler
FC     = mpif90
CC     = mpicc
CXX    = mpicxx

FFLAGS_OMP       = -fopenmp
FFLAGS           = -w -fallow-argument-mismatch -fdefault-real-8 -fdefault-double-8 -DFUNNELED
#FFLAGS          = -w -fallow-argument-mismatch -fdefault-real-8 -fdefault-double-8 -DUSE_BLOCK -DFUNNELED -DWORLDWAR2 -DUNIT_TESTS -DUNIT_TESTS_AFIELDS #-DCOMPARE_ELEMENT_MATRIX
FFLAGS          := $(FFLAGS) -O2 -msse2 -march=native

# --- DEBUGGING
#   --- Debug symbols for debuggers
DEBUGFLAGS  = -g -p
#DEBUGFLAGS += -warn all,nounused -check all,noarg_temp_created -debug all -debug-parameters -fstack-security-check -traceback

FFLAGS_FIXEDFORM := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NOBOUNDS  := $(FFLAGS)               $(FFLAGS_OMP)
FFLAGS           := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)

#FFLAGS_NO_OMP    = $(FFLAGS) $(DEBUGFLAGS)

# Solvers dependencies
USE_STRUMPACK = 1
USE_HIPS   = 0
USE_PASTIX = 0
USE_PASTIX_MURGE = 0
USE_PASTIX6 = 0
USE_MUMPS = 1
USE_WSMP   = 0
USE_FFTW = 1
USE_HDF5 = 1
USE_BOOST = 0
USE_STD_BESSELK = 1
USE_TASKLOOP = 0

#Scotch library
#SCOTCH_HOME  = $(SCOTCH_ROOT)
#LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lptesmumps -lptscotch -lesmumps -lscotch -lptscotcherrexit -lptscotcherr -lscotcherrexit -lscotcherr
#INC_SCOTCH   = -I$(SCOTCH_HOME)/include

LIBLAPACK   = -L$(NETLIB_SCALAPACK_LIB) -lscalapack $(OPENBLAS_LIB)/libopenblas.a

# METIS
#METIS_HOME = $(METIS_ROOT)
LIB_METIS = -L$(METIS_LIB) -lmetis
INC_METIS = -I$(METIS_INCLUDE)

# PARMETIS
#PARMETIS_HOME = $(PARMETIS_ROOT)
LIB_PARMETIS = -L$(PARMETIS_LIB) -lparmetis
INC_PARMETIS = -I$(PARMETIS_INCLUDE)

# STRUMPACK
STRUMPACKLIB = -L$(STRUMPACK_PUBLIC_LEO_LIB64) -lstrumpack $(LIB_PARMETIS) $(LIB_METIS) $(LIBLAPACK)
STRUMPACKINC = -I$(STRUMPACK_PUBLIC_LEO_INC) $(INC_PARMETIS) $(INC_METIS)

# MUMPS
#MUMPS_HOME = $(MUMPS_ROOT)
LIB_MUMPS  = -L$(MUMPS_PUBLIC_LEO_LIB) -lzmumps -ldmumps -lmumps_common $(LIB_PARMETIS) $(LIB_METIS)
INC_MUMPS  = $(MUMPS_PUBLIC_LEO_INC) $(INC_SCOTCH) $(INC_PARMETIS) $(INC_METIS)
ORDLIB     = -L$(MUMPS_PUBLIC_LEO_ROOT)/PORD/lib -lpord

# PASTIX
#LIB_PASTIX        = `$(PASTIX_INSTALL_ROOT)/pastix-conf --libs` -L$(MKLROOT) -lpciaccess
#LIB_PASTIX_MURGE  = `$(PASTIX_INSTALL_ROOT)/pastix-conf --libs_murge`
#LIB_PASTIX_BLAS   = `$(PASTIX_INSTALL_ROOT)/pastix-conf --blas`
#INC_PASTIX        = `$(PASTIX_INSTALL_ROOT)/pastix-conf --incs`

#PASTIX6_HOME      = -I$(PASTIX_ROOT)
#LIB_PASTIX6       = -L$(PASTIX6_HOME)/lib -lpastixf -lspmf -lpastix -lspm -lpastix_kernels $(LIB_SCOTCH)
#LIB_PASTIX6_BLAS  =
#INC_PASTIX6       = $(PASTIX6_HOME)/include $(INC_SCOTCH)

# BOOST
#BOOST_HOME = $(BOOST_ROOT)
#LIB_BOOST  = -L$(BOOST_HOME)/lib -lboost_math_tr1
#INC_BOOST  = -I$(BOOST_HOME)/include

#HIPS
#LIBHIPS   = -L$(HOME)/hips_trunk/trunk/LIB -lhips -lio $(LIBSCOTCH)
#INCHIPS   = -I$(HOME)/hips_trunk/trunk/LIB

#WSMP
#LIB_WSMP = -L$(WSMP_HOME) -lpwsmp64

# LIBFFTW
FFTW_HOME = $(FFTW_ROOT)
LIBFFTW   = -L$(FFTW_LIB) -lfftw3_mpi -lfftw3
INC_FFTW  = -I$(FFTW_INC)

# HDF5
HDF5INCLUDE = $(HDF5_INC)
HDF5LIB     = $(HDF5_LIB)/libhdf5_hl_fortran.a \
              $(HDF5_LIB)/libhdf5_fortran.a \
              $(HDF5_LIB)/libhdf5_hl_f90cstub.a \
              $(HDF5_LIB)/libhdf5_f90cstub.a \
              $(HDF5_LIB)/libhdf5_hl.a \
              $(HDF5_LIB)/libhdf5_tools.a \
              $(HDF5_LIB)/libhdf5.a $(ZLIB_LIB)/libz.a -ldl          
```
