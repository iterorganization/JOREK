---
title: "Pitagora-CPU"
nav_order: 6
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---

# Running JOREK on Pitagora CPU partition hosted by CINECA

## Getting Access

1.  Account Creation via <https://userdb.hpc.cineca.it/> - Create New
    User
2.  Then fill out the HPC-Related information and upload your scanned
    ID/passport
3.  Contact the respective PI to add you to the project

## README

- (For Leonardo) - Slides of the webinar offered by CINECA on the 18th
  of February, 2025 ![slides_qa.pdf](/slides_qa.pdf){.align-center}
- E-Mail Address of CINECA Helpdesk: superc@cineca.it
- The present project related to the EUROFUSION TSVV-F project is called: **FUPA1_REDISMHD**
- Documentation:
  <https://docs.hpc.cineca.it/hpc/pitagora.html#pitagora-card>
- Slides from Cineca (access, modules, etc):
  <https://learn.cineca.it/course/view.php?id=2121>![pitagora_ef_new.pdf](/pitagora_ef_new.pdf){.align-center}c/pitagora.html#pitagora-card)`

## Login

    ssh -Y <user>@login.pitagora.cineca.it

which establishes a connection to one of the available login nodes. You
can also indicate explicitly the login nodes:

    login01-ext.pitagora.cineca.it
    login02-ext.pitagora.cineca.it
    login03-ext.pitagora.cineca.it
    login04-ext.pitagora.cineca.it
    login05-ext.pitagora.cineca.it
    login06-ext.pitagora.cineca.it

Login nodes with odd number (01,03,05) are similar to the CPU partition
of the compute nodes: 2 × AMD EPYC 9745
(<https://www.amd.com/en/products/processors/server/epyc/9005-series/amd-epyc-9745.html>),
128 cores total, 768 GiB DDR5 RAM.

Login nodes with even number (02,04,06) match the compute nodes with
GPU: 2 × Intel Xeon Gold 6548Y+
(<https://www.intel.com/content/www/us/en/products/sku/237564/intel-xeon-gold-6548y-processor-60m-cache-2-50-ghz/specifications.html>),
32 cores total, 512 GiB DDR5 RAM, 1 × NVIDIA H100 NVL GPU
(<https://www.nvidia.com/en-us/data-center/h100/>).

Compile your code on the correct login if you want it to be executed on
a specific compute node partition.

## Software Available

As usual `module avail/load/list/show` allows to access the available
software.

Some non-preinstalled softwares (gnuplot, paraview, visit) are available
at /pitagora_work/FUPA1_MHD/JOREK_LIBRARIES/.

To use the gnuplot, set:

    export PATH=/pitagora_work/FUPA1_MHD/JOREK_LIBRARIES/gnuplot/bin:$PATH

The visit and paraview (copied from Viper) work well on Intel CPU nodes
(login02/04/06), but fail to start on AMD CPU nodes (login01/03/05).

The visit can be enabled by:

    export PATH=/pitagora_work/FUPA1_MHD/JOREK_LIBRARIES/mpcdf/soft/RHEL_9/packages/x86_64/visit/3.4.2/bin:$PATH

The paraview can be enabled by two steps:

    export PATH=/pitagora_work/FUPA1_MHD/JOREK_LIBRARIES/mpcdf/soft/RHEL_9/packages/x86_64/paraview/5.11.2/bin:$PATH
    alias paraview='paraview --mesa' # to disable hardware rendering, as H100 is only for computing

## Compiling JOREK with INTEL oneapi

Below are some examples of the .bashrc and Makefile.inc setup.

- Add the following lines to your `.bashrc` file

<!-- -->

    module purge
    module load git
    module load intel-oneapi-mpi/
    module load intel-oneapi-compilers
    module load intel-oneapi-mkl/2024.0.0--intel-oneapi-mpi--2021.12.1
    module load hdf5/1.14.3--intel-oneapi-mpi--2021.12.1--oneapi--2024.1.0

- The **Makefile.inc** example (note that PaStiX and STRUMPACK should
  not be chosen simultaneously).
- Note that there are two sets of libraries compiled on AMD and Intel
  login nodes, respectively.
- Pastix is recompiled with intel-oneapi-compilers-classic (icc and
  ifort), Strumpack is recompiled with intel-oneapi-compilers (icx and
  ifx). Therefore, STRUMPACK is recommended.
- The MUMPS (use_MUMPS) and the 64-bit (USE_INTSIZE64) libraries have
  not been recompiled yet.
- Only the full MHD models (without particle) have been tested so far
  (until 02.06.2025) for these libraries.

``` Make
MODEL             = model711
FC                = mpiifx
CC                = mpiicx
CXX               = mpiicpx
 
FFLAGS += -mcmodel=medium -O2
FFLAGS += -g -traceback -check bounds  # error tracing

FFLAGS += -march=skylake-avx512 -axCORE-AVX512 # for AMD
# FFLAGS += -xHost                             # for Intel
# FFLAGS += -march=sapphirerapids
# FFLAGS += -double-size 64 -integer-size 32 #-DFUNNELED
# FFLAGS += -fpp -mcmodel=medium -axCORE-AVX512 -align array64byte -qopt-zmm-usage=high

EXTRA_FLAGS += -lstdc++
COMPILER_FAMILY   = intel
# Never use STRUMPACK and PASTIC at the same time 
USE_PASTIX        = 0
USE_PASTIX_MURGE  = 0
USE_PASTIX6       = 0
USE_STRUMPACK     = 1
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
 
# TODO: MUMPS for USE_MUMPS needs to be compiled
# TODO: 64bits libraries for USE_INTSIZE64 need to be compiled
# LIB64DIR =

LIBDIR = /pitagora_work/FUPA1_MHD/JOREK_LIBRARIES/libraries_recompile_amd
# LIBDIR = /pitagora_work/FUPA1_MHD/JOREK_LIBRARIES/libraries_recompile_intel
 
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
 STRUMPACK_HOME=$(LIBDIR)/STRUMPACK/install

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
INC_FFTW          = -I$(MKLROOT)/include/fftw
 
# --- HDF5 Library
HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = -L$(HDF5_HOME)/lib -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz

```

## Compiling JOREK with GNU openmp

Below are some examples of the .bashrc and Makefile.inc setup.

\* Add the following lines to your `.bashrc` file

    module purge

    # GCC + OPENMPI
    module load gcc/12.3.0
    module load openmpi/4.1.6--gcc--12.3.0
    module load openblas/0.3.26--gcc--12.3.0
    module load netlib-scalapack/2.2.0--openmpi--4.1.6--gcc--12.3.0
    module load fftw/3.3.10--openmpi--4.1.6--gcc--12.3.0
    module load boost/1.85.0--openmpi--4.1.6--gcc--12.3.0
    module load hdf5/1.14.3--openmpi--4.1.6--gcc--12.3.0
    module load metis/5.1.0--gcc--12.3.0
    module load parmetis/4.0.3--openmpi--4.1.6--gcc--12.3.0

    # Other module loads
    module load python/3.11.7

``` Make

MODEL             = model600
FC                = mpif90
CC                = mpicc
CXX               = mpicxx

FFLAGS_OMP       = -fopenmp
FFLAGS           = -w -fallow-argument-mismatch -fdefault-real-8 -fdefault-double-8 -DFUNNELED
FFLAGS          := $(FFLAGS) -O2 -msse2 -march=native

DEBUGFLAGS  = -g -p

FFLAGS_FIXEDFORM := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NOBOUNDS  := $(FFLAGS)               $(FFLAGS_OMP)
FFLAGS           := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)

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
USE_MKL           = 0
USE_DIRECT_CONSTRUCTION=0
USE_BICGSTAB      = 0
USE_INTSIZE64     = 0
DEBUG             = 0

LIBDIR_GNU = /pitagora_work/FUPA1_MHD/JOREK_libraries_gnu

ifeq (1, $(USE_MUMPS))

 # METIS
 LIB_METIS = -L$(METIS_LIB) -lmetis
 INC_METIS = -I$(METIS_INCLUDE)

 # PARMETIS
 LIB_PARMETIS = -L$(PARMETIS_LIB) -lparmetis
 INC_PARMETIS = -I$(PARMETIS_INCLUDE)

 MUMPS_DIR = $(LIBDIR_GNU)/MUMPS_v541_recompiled
 LIB_MUMPS  = -L$(MUMPS_DIR)/lib -lzmumps -ldmumps -lmumps_common $(LIB_PARMETIS) $(LIB_METIS) -L${NETLIB_SCALAPACK_LIB} -lscalapack -L${OPENBLAS_LIB} -lopenblas

 INC_MUMPS  = $(MUMPS_DIR)/include $(INC_SCOTCH) $(INC_PARMETIS) $(INC_METIS)
 ORDLIB     = -L$(MUMPS_DIR)/PORD/lib -lpord
endif

ifeq (1, $(USE_STRUMPACK))
 STRUMPACK_HOME= $(LIBDIR_GNU)/STRUMPACK_gnu
 STRUMPACKINC = -I$(STRUMPACK_HOME)/include
 STRUMPACKLIB = -L$(STRUMPACK_HOME)/lib64 -lstrumpack

 STRUMPACKINC += -I${METIS_HOME}/include -I${PARMETIS_HOME}/include
 STRUMPACKLIB += -L${METIS_HOME}/lib -lmetis -L${PARMETIS_HOME}/lib -lparmetis -L${NETLIB_SCALAPACK_LIB} -lscalapack -L${OPENBLAS_LIB} -lopenblas

endif

# --- FFTW Library
INC_FFTW  = -I${FFTW_INCLUDE}
LIBFFTW = -L${FFTW_LIB} -lfftw3_mpi -lfftw3

# --- HDF5 Library
HDF5INCLUDE = ${HDF5_HOME}/include
HDF5LIB     = -L${HDF5_HOME}/lib -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz

```

## Compiling STARWALL with INTEL oneapi

Below are the .bashrc and config.in files required to compile STARWALL
(tested by D. Bonfiglio).

- Add the following lines to your `.bashrc` file

<!-- -->

    module purge
    module load git
    module load intel-oneapi-mpi/
    module load intel-oneapi-compilers
    module load intel-oneapi-mkl/2024.0.0--intel-oneapi-mpi--2021.12.1

- Compile STARWALL (make clean; make -j 8) using this `config.in` file

``` config.in
# config.in file for the LINUX systems
MKL_HOME = $(INTEL_ONEAPI_MKL_HOME)/mkl/latest

include $(FILES_MK)
FC = mpif90 -f90=ifort -fc=ifort
FFLAGS = -autodouble -I$(OBJ_DIR) -module $(OBJ_DIR) -g

ifeq ($(DEBUG),1)
  FFLAGS += -O0 -warn all,nounused -check all,noarg_temp_created -debug all    \
  -debug-parameters -fstack-security-check -ftrapuv -traceback
else
  FFLAGS += -O3 -g
endif

FPPFLAGS =

LD = $(FC)
LDFLAGS = $(FFLAGS)

DPREF = -D
FPPFLAGS += $(DPREF)$(COMPILER) $(DPREF)$(OSTYPE)

########################################################################

MKLLIB = $(MKL_HOME)/lib/intel64
LAPACK_LIB = $(MKLLIB)/libmkl_sequential.a -Wl,--start-group                   \
  $(MKLLIB)/libmkl_intel_lp64.a $(MKLLIB)/libmkl_sequential.a                  \
  $(MKLLIB)/libmkl_core.a $(MKLLIB)/libmkl_scalapack_lp64.a                    \
  $(MKLLIB)/libmkl_blacs_intelmpi_lp64.a    -Wl,--end-group -lpthread

LIB = $(LAPACK_LIB) #-lX11 -lpmapi
########################################################################

all: $(EXECOS_MAIN)

$(EXECOS_MAIN): $(MAIN_OBJ) $(OS_MK) $(FILES_MK)
        $(LD) $(LDFLAGS) -o $@ $(MAIN_OBJ) $(addprefix -L,$(LIB_DIR)) $(LIB)
        mv $@ $(OBJ_DIR)

$(OBJ_DIR)/%.o: %.f90
        $(FC) -fpp $(FPPFLAGS) $(FFLAGS) -c -o $@ $<
```

### Example job script

- The sbatch with a job script seems still not working perfectly (until
  02.07.2025), which returns the following error information: *"sbatch:
  error: Batch job submission failed: Unexpected message received" if
  there are some loaded modules*.
- Please run ***module purge*** before submitting the job via sbatch.
- Load all the modules you need inside the jobscript itself.
- (kuan-wen) For kinetic-MHD codes that uses STRUMPACK as the solver,
  the total number of MPI ranks must be the power of 2. E.g. 2\^6=64.
  Otherwise, the solver will hang forever after the particle loop and
  projection step, because STRUMPACK uses distributed symbolic
  factorization.

``` bash
#!/bin/bash

#SBATCH -J jorek
#SBATCH -p cpu

#SBATCH --nodes=5
#SBATCH --ntasks-per-node=16                   # number of tasks per node
#SBATCH --cpus-per-task=16                     # number of cores per task max 256 cores per node

#SBATCH --mem=760000                           # memory per node

#SBATCH --time=24:00:00
#SBATCH --output=./jorek.%j.out
#SBATCH --error=./jorek.%j.err


### Request e-mail notification
#SBATCH --mail-type=FAIL   # can be BEGIN, END, TIME_LIMIT, FAIL, ALL, NONE
#SBATCH --mail-user=xxxx@ipp.mpg.de

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

## For Intel compiler
module load intel-oneapi-mpi/
module load intel-oneapi-compilers
module load intel-oneapi-mkl/2024.0.0--intel-oneapi-mpi--2021.12.1
module load hdf5/1.14.3--intel-oneapi-mpi--2021.12.1--oneapi--2024.1.0

# This line is critical (at least for PastiX, tbc for Strumpack).
# This infiniband option was changed to default ``export FI_PROVIDER=verbs`` by Cineca team,
# which causes runs to hang forever after a few time-steps (usually during the PastiX solve)
export FI_PROVIDER=MLX

# srun ./jorek_model711 < ./input | tee logfile
mpirun -np $SLURM_NTASKS ./jorek_model711 < ./input | tee logfile
```

- Alternatively, the test job can also be submitted via salloc
  interactively, such as:

``` bash
salloc --nodes=5 --ntasks=80 --cpus-per-task=16 --partition=cpu bash -c 'export OMP_NUM_THREADS=16; mpirun -np 80 ./jorek_model711 < ./input | tee logfile'
```

\<color #ed1c24\> If your simulations are hanging at the matrix
factorization when using Intel and STRUMPACK, try adding the following
line to your jobscript: \</color\>

``` bash
unset UCX_TLS
```

*Use with care, in a few cases this was seen to change simulation
results!*

## Accounting

\<color #ed1c24\> Please be mindful about the CPU consumption, we're
burning up too fast so far! You can check the usage with the saldo
commands below. \</color\>

To check the used resources on Leonardo DCGP, use

`saldo -b –dcgp`

To check more details (containing the usage of all the users) one can
use

`saldo -ra PROJECT_NAME –dcgp`
