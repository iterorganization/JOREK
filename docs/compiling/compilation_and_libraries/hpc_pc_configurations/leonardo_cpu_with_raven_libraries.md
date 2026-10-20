---
title: "Leonardo with Raven Libraries"
nav_order: 4
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---
# Compiling JOREK on Leonardo-CPU with Intel + Oneapi and Raven libaries

The libraries have been copied from Raven and placed in /leonardo_work/FUA38_MHD_0/libraries_from_raven.

If you cannot access it, please download it from https://datashare.mpcdf.mpg.de/s/5KJucm9EnHbsltA.

Please modify the LIBROOT in the Makefile.inc with the correct path.

  ## bashrc
  ```bash
# Load modules for setting-up the environment
# INTEL ONEAPI
module load pkgconf/1.9.5-gr2lx34 zlib-ng/2.1.4--oneapi--2023.2.0
module load intel-oneapi-compilers
module load intel-oneapi-mpi/2021.10.0
module load intel-oneapi-mkl/2023.2.0
module load hdf5/1.14.3--intel-oneapi-mpi--2021.10.0--oneapi--2023.2.0
module load python/3.10.8--gcc--8.5.0
 
# Setting for the new UTF-8 terminal support
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
```

## Makefile.inc 

```bash
# Configuration file for jorek
# model directory
MODEL             = model711
FC                = mpiifort -fc=ifx
CC                = mpiicx
CXX               = mpiicpx
 
FFLAGS += -mcmodel=medium -O2
FFLAGS += -xHost # for Intel
# FFLAGS += -march=skylake-avx512 -axCORE-AVX512 # for AMD
EXTRA_FLAGS += -lstdc++
COMPILER_FAMILY   = intel
# Never use STRUMPACK and PASTIC at the same time 
USE_PASTIX        = 1
USE_PASTIX_MURGE  = 0
USE_PASTIX6       = 0
USE_STRUMPACK     = 0
USE_MUMPS         = 0
USE_BLOCK         = 1
USE_HDF5          = 1
USE_FFTW          = 1
USE_MKL           = 0
USE_DIRECT_CONSTRUCTION=0
USE_BICGSTAB      = 0
USE_INTSIZE64     = 0
DEBUG             = 0
LIBS += -liomp5 -pthread -ldl -lm
 
LIBROOT = /leonardo_work/FUA38_MHD_0/libraries_from_raven
LIBDIR = ${LIBROOT}/ihol_intel_32bits
LIB64DIR = ${LIBROOT}/64bits
 
ifeq (1, $(USE_PASTIX))
 
 MKLLIB += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
 MKLINC = -I${MKLROOT}/include
 PASTIX_DIR = $(LIBDIR)/pastix_5.2.3/install
 SCOTCH_DIR = $(LIBDIR)/scotch_5.1.12
 ifeq (1, $(USE_INTSIZE64))
   PASTIX_DIR = $(LIBDIR)/pastix_5.2.3_i64/install
   SCOTCH_DIR = $(LIBDIR)/scotch_5.1.12_i64
 endif
 
 LIB_PASTIX = -L$(PASTIX_DIR) -lpastix -L$(SCOTCH_DIR)/lib -lscotch -lscotcherr
 INC_PASTIX = -I$(PASTIX_DIR)/include
 
 INC_PASTIX += $(MKLINC)
 LIB_PASTIX += $(MKLLIB)
 
endif
 
ifeq (1, $(USE_MUMPS))
# --- MUMPS
 MKLLIB += -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64
 MKLINC = -I${MKLROOT}/include
 MUMPS_DIR = $(LIBDIR)/MUMPS_5.4.1
 LIB_MUMPS = -L$(MUMPS_DIR)/lib -ldmumps -lmumps_common -lpord
 INC_MUMPS = $(MUMPS_DIR)/include
 
 LIB_MUMPS += $(MKLLIB)
 INC_MUMPS_EXTRA += $(MKLINC)
 
endif
 
ifeq (1, $(USE_STRUMPACK))
 MKLLIB= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64
 MKLINC = -I${MKLROOT}/include
 GKLIB_HOME = $(LIBDIR)/GKlib/install
 METIS_HOME = $(LIBDIR)/METIS/install
 PARMETIS_HOME = $(LIBDIR)/ParMETIS/install
 STRUMPACK_HOME=$(LIBDIR)/STRUMPACK-7.1.3/install
 
 ifeq (1, $(USE_INTSIZE64))
   METIS_HOME = $(LIB64DIR)/METIS_64/install
   PARMETIS_HOME = $(LIB64DIR)/ParMETIS_64/install
   STRUMPACK_HOME=$(LIB64DIR)/STRUMPACK_64/install
 endif
 
 STRUMPACKINC = -I$(STRUMPACK_HOME)/include
 STRUMPACKLIB = -L$(STRUMPACK_HOME)/lib64 -lstrumpack
 
 STRUMPACKINC += -I$(METIS_HOME)/include -I$(GKLIB_HOME)/include -I$(PARMETIS_HOME)/include
 STRUMPACKLIB += -L$(METIS_HOME)/lib -lmetis -L$(GKLIB_HOME)/lib -lGKlib -L$(PARMETIS_HOME)/lib -lparmetis
 
 STRUMPACKINC += $(MKLINC)
 STRUMPACKLIB += $(MKLLIB)
 
 DEFINES += -DNEWSPK
endif
 
# --- FFTW Library
# LIBFFTW           = -L$(FFTW_HOME)/lib -lfftw3
# INC_FFTW          = -I$(FFTW_HOME)/include
INC_FFTW          = -I$(MKLROOT)/include/fftw
 
# --- HDF5 Library
HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = -L$(HDF5_HOME)/lib -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz            
```
