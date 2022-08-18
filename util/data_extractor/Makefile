# --- Makefile for jorek processing with hdf5

# --- Compiler
#FCOMP = pgf90
#FCOMP = f95
#FCOMP = mpif90
FCOMP = mpiifort -fpp

# --- With DEBUG?
#FCOMP = $(FCOMP) -O0 -g -traceback -ftrapuv -debug all -debug-parameters -fstack-security-check -fpe0 -check all,noarg_temp_created -check bounds -check uninit -init=snan -gen-interfaces -warn-interfaces

# --- HDF5 libs
HDF5_INC  = $(HDF5_HOME)/include
HDF5_LIB  = $(HDF5_HOME)/lib/libhdf5hl_fortran.a \
            $(HDF5_HOME)/lib/libhdf5_hl.a \
            $(HDF5_HOME)/lib/libhdf5_fortran.a \
            $(HDF5_HOME)/lib/libhdf5.a \
            $(SZIP_LIB)/libsz.a -lz

all:
	$(FCOMP) -I$(HDF5_INC) -DUSE_HDF5 -c hdf5_io_module.f90
	$(FCOMP) -I$(HDF5_INC) -DUSE_HDF5 -c fortran_process_jorek.f90
	$(FCOMP) -I$(HDF5_INC) fortran_process_jorek.o hdf5_io_module.o -o fortran_process $(HDF5_LIB)
	rm -f *.o *.mod
	
