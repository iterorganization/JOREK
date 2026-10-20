---
title: "Raven"
nav_order: 11
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---

# Raven system at MPCDF

## Basic information

Documentation:

- **[Raven user guide](https://docs.mpcdf.mpg.de/doc/computing/raven-user-guide.html)**
- [Raven at MPCDF](https://www.mpcdf.mpg.de/services/supercomputing/raven)

**Compute nodes:**

- currently 360 nodes
- 2×36 cores per node
- memory per node: 240 GB

### Accessing the cluster

First log into the gate machine:

```bash
ssh gate.mpcdf.mpg.de
```

From the gate machine:

```bash
ssh raven[01-02]i
```

### File systems

AFS is accessible from the login nodes.

Both `/u` and `/ptmp` file systems are accessible from the login and compute nodes.

- Your data can be stored under `/u/<user>`; a quota of 2.5 TB applies  
  (check via `'/usr/lpp/mmfs/bin/mmlsquota raven_u'`)
- `/ptmp/<user>` should be used as temporary storage

**No automatic backup** is done for either file system.

Note that files on `/ptmp/` are **automatically deleted** if they have not been touched for 12 weeks.

## Compiling JOREK

- Add the following lines to your `.bashrc` file:

```bash
module purge
module load intel/2025.1 impi/2021.15
module load mkl/2025.1
module load fftw-mpi
module load hdf5-mpi
module load hwloc/2.2

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${INTEL_HOME}/compiler/latest/linux/compiler/lib/intel64/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MKLROOT}/lib/intel64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HDF5_ROOT}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FFTW_HOME}/lib
```

- In order to use pre-compiled external solver libraries, please send a request to `ihor.holod@ipp.mpg.de` indicating your Raven user name.

- Example `Makefile.inc`  
  (note that PaStiX and STRUMPACK should not be chosen simultaneously):

```makefile
MODEL             = model303
FC                = mpiifx
CC                = mpiicx
CXX               = mpiicpx
 
FFLAGS += -O3 -mcmodel=medium
EXTRA_FLAGS += -lstdc++
COMPILER_FAMILY   = intel
# Never use STRUMPACK and PASTIX at the same time 
USE_PASTIX        = 0
USE_PASTIX_MURGE  = 0
USE_PASTIX6       = 0
USE_STRUMPACK     = 1
USE_MUMPS         = 1
USE_BLOCK         = 1
USE_HDF5          = 1
USE_FFTW          = 1
USE_MKL           = 0
USE_DIRECT_CONSTRUCTION=0
USE_BICGSTAB=0
LIBS += -liomp5 -pthread -ldl -lm
 
LIBDIR = /raven/u/ihol/lib/intel

ifeq (1, $(USE_PASTIX))

 MKLLIB += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
 MKLINC = -I${MKLROOT}/include
 PASTIX_DIR = $(LIBDIR)/pastix_5.2.3/install
 SCOTCH_DIR = $(LIBDIR)/scotch_5.1.12

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
 STRUMPACK_HOME=$(LIBDIR)/STRUMPACK-7.2.0/install

 STRUMPACKINC = -I$(STRUMPACK_HOME)/include
 STRUMPACKLIB = -L$(STRUMPACK_HOME)/lib64 -lstrumpack

 STRUMPACKINC += -I$(METIS_HOME)/include -I$(GKLIB_HOME)/include -I$(PARMETIS_HOME)/include
 STRUMPACKLIB += -L$(METIS_HOME)/lib -lmetis -L$(GKLIB_HOME)/lib -lGKlib -L$(PARMETIS_HOME)/lib -lparmetis

 STRUMPACKINC += $(MKLINC)
 STRUMPACKLIB += $(MKLLIB)

 DEFINES += -DNEWSPK
endif

# --- FFTW Library
LIBFFTW           = -L$(FFTW_HOME)/lib -lfftw3
INC_FFTW          = -I$(FFTW_HOME)/include

# --- HDF5 Library
HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
```

## Submitting jobs

Below is an example SLURM job script:

```bash
#!/bin/bash -l
#SBATCH -D ./
#SBATCH -J jorek
#SBATCH -C icelake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=36
#SBATCH --time=00:10:00

source $HOME/.bashrc

export OMP_NUM_THREADS=36
export OMP_PLACES=cores

srun ./jorek_model303 < intear > logout.out
```

## Running interactively

For running interactively on Raven, go to the interactive nodes by:

```bash
ssh raven-i
```

or

```bash
ssh raven[03-06]
```

There are limits on:

- cores (`<= 8`)
- memory (`<= 256000M`)
- time (`<= 2 hours`)

(My jobs also started with more than the recommended number of cores, but the user guide warns against it.)

```bash
srun -n N_TASKS -p interactive --time=TIME_LESS_THAN_2HOURS --mem=MEMORY_LESS_THAN_32G ./EXECUTABLE
```

Example:

```bash
export OMP_NUM_THREAD=4
nice srun -n 2 -p interactive --time=00:30:00 --mem=16GB ./jorek_model600 < input_jorek | tee logfile
```

It can take a short time to allocate the resources.
