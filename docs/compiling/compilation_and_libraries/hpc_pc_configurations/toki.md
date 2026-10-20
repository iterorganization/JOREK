---
title: "TOK-I Cluster"
nav_order: 9
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---

# Run JOREK on the TOK-I Cluster at IPP Garching

> **This page was recently updated after the update of the cluster to SLES15. If you find anything that is incorrect with the new system, please correct it.**

- [Information about the cluster in the TOK Wiki](https://tok.ipp.mpg.de/wiki/TOK/index.php/Computing-TOKI)

```bash
ssh toki01.bc.rzg.mpg.de
...
ssh toki08.bc.rzg.mpg.de
```

## Modules required

```bash
module load git intel/2023.1.0.x impi/2021.9 mkl fftw-mpi gnuplot hdf5-mpi
```

Add to your `~/.bashrc` or `~/.alias` file.

```bash
module load impi-interactive
```

It is recommended **not** to add `impi-interactive` to your `.bashrc`, as this may stop simulations on the TOK cluster from working.

## Compiling JOREK and Libraries

- [Installation/compilation of libraries and JOREK](../compiling.md) works in the normal way
- You can simply use the libraries available via:

```bash
module load pastix
```

> **Warning:** There seems to be an issue when JOREK is run on clusters at IPP (TOK-I, TOK, DRACO, and COBRA) with more than one MPI task per harmonic. The run hangs at the first PASTIX call during the very first time-step. In such a situation, the combination of `Pastix_5.2.3` and `Scotch_5.1.12b` has been seen to fix this issue.

The following `Makefile.inc` file should allow you to compile JOREK:

```makefile
# JOREK configuration file for compilation on the TOK cluster at IPP Garching
# see also https://www.jorek.eu/wiki/doku.php?id=draco
 
MODEL             = model307
 
FC                = mpiifort
CC                = mpiicc
CXX               = mpiicpc
 
COMPILER_FAMILY   = intel
 
USE_PASTIX        = 1
USE_STRUMPACK     = 0
USE_MUMPS         = 1
USE_PASTIX_MURGE  = 0
USE_BLOCK         = 1
USE_HDF5          = 1
USE_FFTW          = 1
USE_MKL           = 0
 
 
LIBDIR = /tokp/work/ihol/bin/intel
 
ifeq (1, $(USE_MUMPS))
 MKLLIB += -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64
 MKLINC = -I${MKL_HOME}/include
 
 MUMPSDIR = $(LIBDIR)/MUMPS_5.2.1
 LIB_MUMPS = -L$(MUMPSDIR)/lib -ldmumps -lmumps_common -lpord
 INC_MUMPS = $(MUMPSDIR)/include
 
 LIB_MUMPS += $(MKLLIB)
 INC_MUMPS += $(MKLINC)
endif
 
 
ifeq (1, $(USE_PASTIX))
 MKLLIB += -L${MKL_HOME}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
 MKLINC = -I${MKL_HOME}/include
 
 PASTIX_DIR = $(LIBDIR)/pastix_5.2.3/install
 SCOTCH_DIR = $(LIBDIR)/scotch_5.1.12
 
 LIB_PASTIX = -L$(PASTIX_DIR) -lpastix -L$(SCOTCH_DIR)/lib -lscotch -lscotcherrexit
 INC_PASTIX = -I$(PASTIX_DIR)/include
 
 INC_PASTIX += $(MKLINC)
 LIB_PASTIX += $(MKLLIB)
 
endif
 
ifeq (1, $(USE_STRUMPACK))
 MKLLIB = -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64
 MKLINC = -I${MKLROOT}/include
 
 METIS_HOME=$(LIBDIR)/metis_git/install
 PARMETIS_HOME=$(LIBDIR)/ParMETIS_git/install
 GKLIB_HOME=$(LIBDIR)/GKlib_git/install
 
 STRUMPACK_HOME=$(LIBDIR)/STRUMPACK-7.1.3/install
 
 STRUMPACKINC = -I$(STRUMPACK_HOME)/include
 STRUMPACKLIB = -L$(STRUMPACK_HOME)/lib64 -lstrumpack
 
 STRUMPACKINC += -I$(METIS_HOME)/include -I$(GKLIB_HOME)/include
 STRUMPACKLIB += -L$(METIS_HOME)/lib -lmetis -L$(GKLIB_HOME)/lib -lGKlib
 
 STRUMPACKINC += -I$(PARMETIS_HOME)/include
 STRUMPACKLIB += -L$(PARMETIS_HOME)/lib -lparmetis
 
 STRUMPACKINC += $(MKLINC)
 STRUMPACKLIB += $(MKLLIB)
 
 DEFINES += -DNEWSPK
endif
 
# --- FFTW Library
#LIBFFTW           = -Wl -rpath -L$(FFTW_HOME)/lib -lfftw3
LIBFFTW           = $(FFTW_HOME)/lib/libfftw3.a
INC_FFTW          = -I$(FFTW_HOME)/include
 
# --- HDF5 Library
HDF5INCLUDE = $(HDF5_HOME)/include
#HDF5LIB     = -Wl -rpath -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
HDF5LIB     = $(HDF5_HOME)/lib/libhdf5hl_fortran.a  $(HDF5_HOME)/lib/libhdf5_hl.a  $(HDF5_HOME)/lib/libhdf5_fortran.a $(HDF5_HOME)/lib/libhdf5.a -L$(HDF5_HOME)/lib  -lz
```

## File systems

- Best is to use `/tokp/work/<user>/` for your runs
- Also possible:
  - `/tokp/scratch/<user>/` — not backed up
  - `/ptmp1/work/<user>/`
  - `/ptmp1/scratch/<user>/` — not backed up

## Batch submission

_Not available._

## Interactive Runs

- Check with `top` beforehand that the machine is not already heavily used by others
- Do not run very long simulations
- Since each node has 24 cores, a rule of thumb is that OpenMP threads × MPI tasks should not exceed 20

```bash
export OMP_NUM_THREADS=<number-of-OpenMP-threads>
nice mpirun -n <mpi-tasks> ./jorek_model199 < ./input | tee logfile
```

When running on TOK-I, if you get the error below, specify the location of the JOREK binary:

```bash
env: ‘jorek_model199’: No such file or directory
```

### Without MPI

To run axisymmetric simulations, use only one MPI task:

```bash
mpirun -n 1 ./jorek_model199 < ./input | tee logfile
```

Since MPI is always included in JOREK, this is the closest one can get to running without MPI.

### Interactive regression tests

To run regression tests interactively, set:

```bash
export MPIRUN="mpirun -n 1"
```

and then perform them in the usual way using:

```bash
./non_regression_tests/run_test.sh
```
