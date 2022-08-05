# Makefile.inc file for compiling JOREK on the Marconi100 Supercomputer
# including CUDA FORTRAN kernels

# --- Select physics model
MODEL             = model307


# --- Compiler and options
FC                = mpipgifort
CC                = mpipgicc
CXX               = mpipgic++

##FLAGS += -g

COMPILER_FAMILY=pg

DEFINES += -DCUDA_KERNELS

# --- Various switches
USE_PASTIX              = 1
USE_STRUMPACK           = 0
USE_MUMPS               = 1
USE_PASTIX_MURGE        = 0
USE_COMPLEX_PRECOND     = 0
USE_HDF5                = 1
USE_FFTW                = 1
USE_DIRECT_CONSTRUCTION = 0
USE_BLOCK               = 1
USE_MPI                 = 1

LIBDIR = /m100/home/userexternal/iholod00/bin
LIBS += -pthread 

##################

ifeq (1, $(USE_PASTIX))

 PASTIX_DIR=$(LIBDIR)/pastix_5.2.3-fun/install

# This should work but leaves a trailing -L so will not compile
# LIB_PASTIX = $(shell $(PASTIX_DIR)/pastix-conf --libs)
# Use output from pastix-conf directly
LIB_PASTIX = -L/m100/home/userexternal/iholod00/bin/pastix_5.2.3-fun/src/../install -lpastix -lm -lrt -L/m100/home/userexternal/iholod00/bin/scotch_6.0.8//lib -lscotch -lscotcherrexit -lhwloc -lpthread -L/cineca/prod/opt/libraries/lapack/3.9.0/gnu--8.4.0/lib -llapack -L/cineca/prod/opt/libraries/blas/3.8.0/gnu--8.4.0/lib -lblas
 INC_PASTIX = $(shell $(PASTIX_DIR)/pastix-conf --incs)

endif

##################

ifeq (1, $(USE_STRUMPACK))

 MKLLIB = -L${SCALAPACK_LIB} -lscalapack -L${BLAS_LIB} -lblas -L${LAPACK_LIB} -llapack
 MKLINC = -I${LAPACK_INC}

 MKLLIB += -lmpi_ibm_usempi -lmpi_ibm_mpifh -lgfortran -lquadmath -lrt -ldl

 PARMETIS_HOME=$(LIBDIR)/parmetis-4.0.3--spectrum_mpi/build
 STRUMPACK_HOME=$(LIBDIR)/STRUMPACK/install

 STRUMPACKINC = -I$(STRUMPACK_HOME)/include
 STRUMPACKLIB = $(STRUMPACK_HOME)/lib/libstrumpack.a

 STRUMPACKINC += -I$(PARMETIS_HOME)/include
 STRUMPACKLIB += $(PARMETIS_HOME)/libparmetis/libparmetis.a

 STRUMPACKINC += -I$(PARMETIS_HOME)/metis/include
 STRUMPACKLIB += $(PARMETIS_HOME)/libmetis/libmetis.a

 STRUMPACKINC += $(MKLINC)
 STRUMPACKLIB += $(MKLLIB)

 DEFINES += -DNEWSPK

endif

##################

ifeq (1, $(USE_MUMPS))

 MUMPS_HOME = $(LIBDIR)/MUMPS_5.2.1
 PARMETIS_HOME=$(LIBDIR)/parmetis-4.0.3/build

 MKLLIB += -L${SCALAPACK_LIB} -lscalapack -L${BLAS_LIB} -lblas -L${LAPACK_LIB} -llapack

 LIB_MUMPS = $(MUMPS_HOME)/lib/libdmumps.a $(MUMPS_HOME)/lib/libmumps_common.a
 INC_MUMPS = $(MUMPS_HOME)/include

 INC_MUMPS_EXTRA += -I$(PARMETIS_HOME)/include
 ORDLIB += $(PARMETIS_HOME)/libparmetis/libparmetis.a

 INC_MUMPS_EXTRA += -I$(PARMETIS_HOME)/metis/include
 ORDLIB += $(PARMETIS_HOME)/libmetis/libmetis.a
 
 LIB_MUMPS += $(MKLLIB)

endif



####################


LIBFFTW           = ${FFTW_LIB}/libfftw3.a
INC_FFTW          = -I${FFTW_INC}

HDF5INCLUDE = ${HDF5_INC}
HDF5LIB     =-L${HDF5_LIB} -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -L${ZLIB_LIB} -lz
