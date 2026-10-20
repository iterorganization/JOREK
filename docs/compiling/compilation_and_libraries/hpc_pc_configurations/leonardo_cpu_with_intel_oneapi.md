---
title: "Leonardo with Intel oneAPI"
nav_order: 3
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---
# Compiling JOREK on Leonardo-CPU with Intel and Oneapi
## bashrc file
```bash
# Load modules for setting-up the environment
# INTEL ONEAPI
module load intel-oneapi-compilers
module load intel-oneapi-vtune
module load intel-oneapi-inspector
module load intel-oneapi-itac
module load intel-oneapi-mpi/2021.10.0
module load intel-oneapi-tbb/2021.10.0
module load intel-oneapi-mkl/2023.2.0
module load hdf5/1.14.3--intel-oneapi-mpi--2021.10.0--oneapi--2023.2.0
module load python/3.10.8--gcc--8.5.0

#INTEL compiled libraries
export INTEL_GKLIB_LIB_LOC="/leonardo_work/FUAL8_MHD/csommari/install/GKlib_oneapi202320_intelmpi2021_23012025/lib"
export INTEL_GKLIB_INC_LOC="/leonardo_work/FUAL8_MHD/csommari/install/GKlib_oneapi202320_intelmpi2021_23012025/include"
export INTEL_GKLIB_ROOT_LOC="/leonardo_work/FUAL8_MHD/csommari/install/GKlib_oneapi202320_intelmpi2021_23012025"
export INTEL_PARMETIS_LIB_LOC="/leonardo_work/FUAL8_MHD/csommari/install/ParMETIS_github_oneapi202320_intelmpi2021_gklib_metis_v521_i32_r64_23012025/lib"
export INTEL_METIS_LIB_LOC="/leonardo_work/FUAL8_MHD/csommari/install/METIS_v521_i32_r64_oneapi202320_intelmpi2021_gklib_23012025/lib"
export INTEL_PARMETIS_INC_LOC="/leonardo_work/FUAL8_MHD/csommari/install/ParMETIS_github_oneapi202320_intelmpi2021_gklib_metis_v521_i32_r64_23012025/include"
export INTEL_METIS_INC_LOC="/leonardo_work/FUAL8_MHD/csommari/install/METIS_v521_i32_r64_oneapi202320_intelmpi2021_gklib_23012025/include"
export INTEL_PARMETIS_ROOT_LOC="/leonardo_work/FUAL8_MHD/csommari/install/ParMETIS_github_oneapi202320_intelmpi2021_gklib_metis_v521_i32_r64_23012025"
export INTEL_METIS_ROOT_LOC="/leonardo_work/FUAL8_MHD/csommari/install/METIS_v521_i32_r64_oneapi202320_intelmpi2021_gklib_23012025"
export INTEL_MUMPS_ROOT_LOC="/leonardo_work/FUAL8_MHD/csommari/install/MUMPS_v541_oneapi202320_intelmpi2021_intelmkl202320_gklib_metis_v531_parmetis_i32_r64_23012025"
export INTEL_MUMPS_LIB_LOC="/leonardo_work/FUAL8_MHD/csommari/install/MUMPS_v541_oneapi202320_intelmpi2021_intelmkl202320_gklib_metis_v531_parmetis_i32_r64_23012025/lib"
export INTEL_MUMPS_INC_LOC="/leonardo_work/FUAL8_MHD/csommari/install/MUMPS_v541_oneapi202320_intelmpi2021_intelmkl202320_gklib_metis_v531_parmetis_i32_r64_23012025/include"
export INTEL_STRUMPACK_ROOT_LOC="/leonardo_work/FUAL8_MHD/csommari/install/STRUMPACK_v710_oneapi202320_intelmpi2021_intelmkl202320_gklib_metis_v531_parmetis_i32_r64_23012025"
export INTEL_STRUMPACK_INC_LOC="/leonardo_work/FUAL8_MHD/csommari/install/STRUMPACK_v710_oneapi202320_intelmpi2021_intelmkl202320_gklib_metis_v531_parmetis_i32_r64_23012025/include"
export INTEL_STRUMPACK_LIB_LOC="/leonardo_work/FUAL8_MHD/csommari/install/STRUMPACK_v710_oneapi202320_intelmpi2021_intelmkl202320_gklib_metis_v531_parmetis_i32_r64_23012025/lib/cmake/STRUMPACK"
export INTEL_STRUMPACK_LIB64_LOC="/leonardo_work/FUAL8_MHD/csommari/install/STRUMPACK_v710_oneapi202320_intelmpi2021_intelmkl202320_gklib_metis_v531_parmetis_i32_r64_23012025/lib64"
export INTEL_FFTW_ROOT_LOC="/leonardo_work/FUAL8_MHD/csommari/install/FFTW_v3310_oneapi202320_intelmpi2021_24012025"
export INTEL_FFTW_LIB_LOC="/leonardo_work/FUAL8_MHD/csommari/install/FFTW_v3310_oneapi202320_intelmpi2021_24012025/lib"
export INTEL_FFTW_INC_LOC="/leonardo_work/FUAL8_MHD/csommari/install/FFTW_v3310_oneapi202320_intelmpi2021_24012025/include"
```

## Makefile.inc file 
```bash
# model directory
MODEL = model600

# Fortran compiler
FC         = mpiifort -fc=ifx
CC         = mpiicx
CXX        = mpiicpx
FFLAGS_OMP = -qopenmp
FFLAGS     = -double-size 64 -integer-size 32 -DFUNNELED
#FFLAGS    = -double-size 64 -integer-size 32 -DUSE_BLOCK -DFUNNELED -DWORLDWAR2 -DUNIT_TESTS -DUNIT_TESTS_AFIELDS #-DCOMPARE_ELEMENT_MATRIX
FFLAGS    := $(FFLAGS) -O2 -fpp -mcmodel=medium -axCORE-AVX512 -mtune=icelake-server -align array64byte -qopt-zmm-usage=high
#FFLAGS    := $(FFLAGS) -Ofast -fpp -mcmodel=medium -axCORE-AVX512 -mtune=icelake-server -align array64byte -qopt-zmm-usage=high

# --- DEBUGGING
#   --- Debug symbols for debuggers
#DEBUGFLAGS      = -g
#DEBUGFLAGS      += -warn all,nounused
#DEBUGFLAGS      += -check all,noarg_temp_created
#DEBUGFLAGS      += -debug all
#DEBUGFLAGS      += -debug-parameters
#DEBUGFLAGS      += -fstack-security-check
#DEBUGFLAGS      += -ftrapuv

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
LIBLAPACK   = -L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl

# GKlib
GKLIB_LIB = -L$(INTEL_GKLIB_LIB_LOC) -lGKlib
GKLIB_INC = -I$(INTEL_GKLIB_INC_LOC)

# METIS
#METIS_HOME = $(METIS_ROOT)
LIB_METIS = -L$(INTEL_METIS_LIB_LOC) -lmetis $(GKLIB_LIB)
INC_METIS = -I$(INTEL_METIS_INC_LOC) $(GKLIB_INC)

# PARMETIS
#PARMETIS_HOME = $(PARMETIS_ROOT)
LIB_PARMETIS = -L$(INTEL_PARMETIS_LIB_LOC) -lparmetis
INC_PARMETIS = -I$(INTEL_PARMETIS_INC_LOC)

# STRUMPACK
STRUMPACKLIB = -L$(INTEL_STRUMPACK_LIB64_LOC) -lstrumpack $(LIB_PARMETIS) $(LIB_METIS) $(LIBLAPACK)
STRUMPACKINC = -I$(INTEL_STRUMPACK_INC_LOC) $(INC_PARMETIS) $(INC_METIS)

# MUMPS
#MUMPS_HOME = $(MUMPS_ROOT)
LIB_MUMPS  = -L$(INTEL_MUMPS_LIB_LOC) -lzmumps -ldmumps -lmumps_common $(LIB_PARMETIS) $(LIB_METIS)
INC_MUMPS  = $(INTEL_MUMPS_INC_LOC) $(INC_SCOTCH) $(INC_PARMETIS) $(INC_METIS)
ORDLIB     = -L$(INTEL_MUMPS_ROOT_LOC)/PORD/lib -lpord

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
FFTW_HOME = $(INTEL_FFTW_ROOT_LOC)
LIBFFTW   = -L$(INTEL_FFTW_LIB_LOC) -lfftw3_mpi -lfftw3_omp -lfftw3_threads -lfftw3
INC_FFTW  = -I$(INTEL_FFTW_INC_LOC)

# HDF5
HDF5INCLUDE = -I$(HDF5_INC)
HDF5LIB     = $(HDF5_LIB)/libhdf5_hl_fortran.a \
              $(HDF5_LIB)/libhdf5_fortran.a \
              $(HDF5_LIB)/libhdf5_hl_f90cstub.a \
              $(HDF5_LIB)/libhdf5_f90cstub.a \
              $(HDF5_LIB)/libhdf5_hl.a \
              $(HDF5_LIB)/libhdf5_tools.a \
              $(HDF5_LIB)/libhdf5.a $(ZLIB_NG_LIB)/libz.a -ldl
              
```
