---
title: "Viper"
nav_order: 12
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---

# Viper supercomputer system at MPCDF

- 768 CPU compute nodes  
  (AMD EPYC Genoa 9554 CPUs with 128 cores and at least 512 GB RAM per node, 4.9 PFlop/s theoretical peak performance (FP64))
- 228 GPU nodes comprising 456 APUs  
  (to be deployed in the course of 2024)

- <https://docs.mpcdf.mpg.de/doc/computing/viper-user-guide.html>

## Login via gateway

```bash
ssh <user>@gate.mpcdf.mpg.de
```

```bash
ssh <user>@viper.mpcdf.mpg.de
```

## File systems

The file system `/u` (a symbolic link to `/viper/u`) is designed for permanent user data (source files, config files, etc.). The size of `/u` is 1.2 PB. Your home directory is in `/u`. The default disk quota in `/u` is 1.0 TB, and the file quota is 256K files. You can check your disk quota in `/u` with:

```bash
/usr/lpp/mmfs/bin/mmlsquota viper_u
```

The file system `/ptmp` (a symbolic link to `/viper/ptmp`) is designed for batch job I/O (12 PB, no system backups). Files in `/ptmp` that have not been accessed for more than 12 weeks will be removed automatically. The period of 12 weeks may be reduced if necessary, with prior notice.

As a current policy, no quotas are applied on `/ptmp`.

The `/r` file system (a symbolic link to `/ghi/r`) stages archive data. It is available only on the login nodes `viper.mpcdf.mpg.de` and on the interactive nodes `viper-i.mpcdf.mpg.de`.

Each user has a subdirectory `/r/<initial>/<userid>` to store data. For efficiency, files should be packed into tar files, with a size of about 1 GB to 1 TB, before archiving them in `/r`. Please avoid archiving many small files. When the file system `/r` becomes filled above a certain threshold, files will be transferred from disk to tape, beginning with the largest files that have not been used for the longest time.

Please do not use the file system `/tmp` or `$TMPDIR` for scratch data. Instead, use `/ptmp`, which is accessible from all Viper cluster nodes.

In cases where an application really depends on node-local storage, please use the directories from the environment variables `JOB_TMPDIR` and `JOB_SHMTMPDIR`, which are set individually for each Slurm job.

## Compiling JOREK

Add the following lines to your `.bashrc` file:

```bash
module purge
module load intel/2024.0 impi/2021.11 mkl/2024.0
module load hdf5-mpi/1.14.1
module load fftw-mpi/3.3.10

export LD_LIBRARY_PATH=${FFTW3_ROOT}/lib:${HDF5_ROOT}/lib:${MKLROOT}/lib:${LD_LIBRARY_PATH}
```

Example `Makefile.inc`:

```makefile
MODEL             = model303

FC                = mpiifx
CC                = mpiicx
CXX               = mpiicpx

FLAGS += -O2 -march=znver4
COMPILER_FAMILY   = intel

USE_PASTIX        = 0
USE_PASTIX_MURGE  = 0
USE_PASTIX6       = 0
USE_STRUMPACK     = 1
USE_MUMPS         = 0
USE_BLOCK         = 1
USE_HDF5          = 1
USE_FFTW          = 1
USE_MKL           = 0

LIBDIR = /u/ihol/lib/intel

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

## Running JOREK

Example batch job script:

```bash
#!/bin/bash -l
#SBATCH -D ./
#SBATCH -J JOREK
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores

export EXECUTABLE="./jorek_model303"

srun ${EXECUTABLE} <inxflow |& tee logfile.out
```
