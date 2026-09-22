---
title: "EUROfusion Gateway"
nav_order: 14
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---
## The Hardware resources of the Gateway independent cluster

9 virtual nodes including the following Gateway user servers:

- ITM MySQL database `itmcatalog`: `itmmysql1.eufus.eu`
- `gforge.eufus.eu`
- `portal.eufus.eu`
- `www.eufus.eu`
- `ticket.eufus.eu`
- 1 node `x3550 M5` with FC HBA as AFS file server of cell `eufus.eu`
- 1 node `x3550 M5` with FC HBA as login node and spare of AFS server
- 4 nodes `thinkSystem SD530` as login nodes with advanced graphics adapter NVIDIA K80
- 16 nodes `thinkSystem SD530` as compute nodes with 192 GB of RAM
- 8 nodes `thinkSystem SD530` as compute nodes with 384 GB of RAM
- 10 TB of AFS file server storage area
- 15 PB of backup storage area for AFS/GPFS
- 600 TB of MARCONI GSS as high-performance storage area

The `thinkSystem SD530` node is based on 2 × 24-core Intel Xeon 8160 (SkyLake) processors at 2.10 GHz.

## Accessing the system

- The following link describes **how to get access** to the Gateway machine: <https://wiki.eufus.eu/doku.php>
- See also the material of a mini CodeCamp:
  - [here](https://box.pionier.net.pl/d/6a7d79934c2b4a4db7e2/) (Password: `miniC0d3C4mp5`)
  - [here](https://docs.psnc.pl/pages/viewpage.action?pageId=70884347)

The machine can then be accessed with an SSH client:

```bash
ssh -CXY <username>@login.eufus.eu
```

## Compiling JOREK

The following modules need to be loaded, for example by including them in `~/.bashrc_profile`:

```bash
module purge
module load cineca
module load intel/pe-xe-2017--binary
module load intelmpi/2017--binary
module load mkl/2017--binary
module load fftw/3.3.4--intelmpi--2017--binary
module load szip/2.1--gnu--6.1.0
module load zlib/1.2.8--gnu--6.1.0
module load hdf5/1.8.17--intelmpi--2017--binary
```

The `Makefile.inc` for compiling JOREK is as follows:

```bash
# --- Select physics model
MODEL             = model303

# --- Compiler and options
FC                = mpiifort
CC                = mpiicc
CXX               = mpiicpc

FFLAGS += -O3
FFLAGS += -vecabi=compat -mcmodel=medium

# --- Various switches
USE_PASTIX              = 1
USE_PASTIX_MURGE        = 0
USE_COMPLEX_PRECOND     = 0
USE_STRUMPACK           = 0
USE_HDF5                = 1
USE_FFTW                = 1
USE_MUMPS               = 0
USE_DIRECT_CONSTRUCTION = 0
USE_BLOCK               = 1

LIBDIR = /afs/eufus.eu/user/g/g2iholod/public
LIBS += -liomp5 -pthread -ldl -lm

ifeq (1, $(USE_MUMPS))
# --- MUMPS
 MKLLIB= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64
 MKLINC = -I${MKL_HOME}/include

 MUMPSDIR = $(LIBDIR)/MUMPS_5.2.1
 LIB_MUMPS = -L$(MUMPSDIR)/lib -ldmumps -ldmumps -lmumps_common -lpord
 INC_MUMPS = $(MUMPSDIR)/include

 LIB_MUMPS += $(MKLLIB)
 INC_MUMPS += $(MKLINC)
endif

##################

ifeq (1, $(USE_PASTIX))
 MKLLIB = -L${MKL_HOME}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
 MKLINC = -I${MKL_HOME}/include

 INC_PASTIX += $(MKLINC)
 LIB_PASTIX += $(MKLLIB)

 PASTIX_DIR=$(LIBDIR)/pastix_5.2.3/install
 SCOTCH_DIR=$(LIBDIR)/scotch_5.1.12

 LIB_PASTIX += $(PASTIX_DIR)/libpastix.a $(SCOTCH_DIR)/lib/libscotch.a $(SCOTCH_DIR)/lib/libscotcherr.a
 INC_PASTIX += -I$(PASTIX_DIR)

endif

##################

ifeq (1, $(USE_STRUMPACK))
 MKLLIB= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64
 MKLINC = -I${MKL_HOME}/include

 PARMETIS_HOME=$(LIBDIR)/parmetis-4.0.3/build
 STRUMPACK_HOME=$(LIBDIR)/STRUMPACK_MKL/install

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

LIBFFTW           = $(FFTW_LIB)/libfftw3f.a
INC_FFTW          = -I$(FFTW_INC)

HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     =-L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -L$(SZIP_LIB) -lsz -lz
```

## Submitting jobs

Example batch job submission script:

```bash
#!/bin/bash

#SBATCH --job-name=jorek
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:10:00
#SBATCH --partition=gw

export OMP_NUM_THREADS=24

srun ./jorek_model303 < ./inxflow > logfile.out
```
