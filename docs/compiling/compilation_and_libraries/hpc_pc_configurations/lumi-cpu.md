---
title: "LUMI-CPU"
nav_order: 5
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---
 
# Running JOREK on LUMI CPU partition (LUMI-C)

## Getting Access
- Contact the respective PI to add you to the project
- Set up [SSH keys](https://docs.lumi-supercomputer.eu/firststeps/SSH-keys/)
  - `ssh-keygen -t ed25519 -f /home/username/.ssh/id_rsa_lumi` or `ssh-keygen -t rsa -b 4096 -f /home/username/.ssh/id_rsa_lumi`
  - They suggested not to leave the passphrase empty.
  - After generating the key pair, register the public key in your [MyAccessID user profile](https://mms.myaccessid.org/fed-apps/profile/).

## README/Useful links

- [LUMI First steps](https://docs.lumi-supercomputer.eu/firststeps/)
- [LUMI-C hardware overview](https://docs.lumi-supercomputer.eu/hardware/lumic/)
- E-Mail Address of Helpdesk: [LUMI User Support](https://lumi-supercomputer.eu/user-support/need-help/general/)

## Login

```bashrc
ssh -i <path-to-private-key> <username>@lumi.csc.fi
```

Example:

```bashrc
ssh -i ~/.ssh/id_rsa_lumi username@lumi.csc.fi
```

LUMI has several login nodes, for reliability and for sharing the interactive workload. The name lumi.csc.fi points automatically to one of the login nodes - the IP address to which it resolves will belong to one of:

| Login node name | IP Address |
|----------------|------------|
| lumi-uan01.csc.fi | 193.167.209.163 |
| lumi-uan02.csc.fi | 193.167.209.164 |
| lumi-uan03.csc.fi | 193.167.209.165 |
| lumi-uan04.csc.fi | 193.167.209.166 |

## Compiling JOREK by GNU (under tests)
Only mumps solver is available now.

.bashrc configuration:

```bash
test -s ~/.alias && . ~/.alias || true

# Variables
#export DEBUG=1
export WORK="/projappl/project_465002663"
export SCRATCH="/scratch/project_465002663"
export EBU_USER_PREFIX="$WORK/EasyBuild"
export MPIRUN="srun -N"
export OMP_NUM_THREADS=128
#export OMP_STACKSIZE=24M
export OMP_DISPLAY_ENV=true
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export PMI_NO_PREINITIALIZE=y
export CRAY_OMP_CHECK_AFFINITY=FALSE
export MPICH_MAX_THREAD_SAFETY=funneled
export KMP_AFFINITY=compact
#GNU Compiler Variables
export MUMPS_ROOT_GNU="$EBU_USER_PREFIX/SW/LUMI-25.03/C/MUMPS/5.6.1-cpeGNU-25.03-OpenMP"
export MUMPS_LIB_GNU="$MUMPS_ROOT/lib"
export MUMPS_INC_GNU="$MUMPS_ROOT/include"
export STRUMPACK_ROOT_GNU="$WORK/wangchiz/STRUMPACK/install"
export STRUMPACK_LIB_GNU="$STRUMPACK_ROOT/lib64"
export STRUMPACK_INC_GNU="$STRUMPACK_ROOT/include"
#export PASTIX_ROOT_GNU="$WORK/sommariv/PASTIX_523_cpeGNU_ptscotch_1004205"
#export GKLIB_ROOT_GNU="$WORK/sommariv/GKlib_master_cpeGNU_openmp_14042025"
#export GKLIB_LIB_GNU="$GKLIB_ROOT_GNU/lib"
#export GKLIB_INC_GNU="$GKLIB_ROOT_GNU/include"
export METIS_ROOT_GNU="$EBROOTPARMETIS"
export METIS_LIB_GNU="$METIS_ROOT_GNU/lib64"
export METIS_INC_GNU="$METIS_ROOT_GNU/include"
export PARMETIS_ROOT_GNU="$EBROOTPARMETIS"
export PARMETIS_LIB_GNU="$PARMETIS_ROOT_GNU/lib"
export PARMETIS_INC_GNU="$PARMETIS_ROOT_GNU/include"
export SCOTCH_ROOT_GNU="$EBROOTSCOTCH"
export SCOTCH_LIB_GNU="$SCOTCH_ROOT_GNU/lib"
export SCOTCH_INC_GNU="$SCOTCH_ROOT_GNU/include"
#export HWLOC_ROOT_GNU="$WORK/sommariv/HWLOC_2112_cpeGNU"
#export HWLOC_LIB_GNU="$HWLOC_ROOT_GNU/lib"
#export HWLOC_INC_GNU="$HWLOC_ROOT_GNU/include"

# Module load GNU
module purge
module load LUMI/25.03
module load partition/C
module load cpeGNU
module load zlib
module load wget
module load util-linux
module load systools
module load syslibs
module load nano
module load libreadline
module load libtirpc
module load libunistring
module load gettext
module load cURL
module load Vim
module load Tcl
module load ParMETIS
module load GLib
module load Boost
module load cray-hdf5-parallel
module load zlib
module load cray-libsci
module load cray-fftw
module load cray-python
module load SCOTCH
module load EasyBuild-user
module load MUMPS/5.6.1-cpeGNU-25.03-OpenMP
module load gnuplot
```
Makefile.inc configuration:

```bashrc
# Configuration file for jorek
# model directory
MODEL = model600

# Fortran compiler 
FC     = ftn
CC     = cc
CXX    = CC

FFLAGS_OMP = -fopenmp
FFLAGS     = -O2 -fdefault-real-8 -fdefault-double-8 -mcmodel=medium
FFLAGS     += -mavx2 -fallow-argument-mismatch -DFUNNELED

# --- DEBUGGING
#   --- Debug symbols for debuggers
#DEBUGFLAGS  :=
#DEBUGFLAGS  = -g #-fsave-loopmark -ffpe-trap=fp,underflow,denormal

FFLAGS_FIXEDFORM := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NOBOUNDS  := $(FFLAGS)               $(FFLAGS_OMP)
FFLAGS           := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)

#FFLAGS_NO_OMP    = $(FFLAGS) $(DEBUGFLAGS)


# Solvers dependencies
USE_HIPS         = 0
USE_PASTIX       = 0
USE_PASTIX_MURGE = 0
USE_PASTIX6      = 0
USE_STRUMPACK    = 1
USE_MUMPS        = 1
USE_WSMP         = 0
USE_FFTW         = 1
USE_HDF5         = 1
USE_BOOST        = 0
USE_STD_BESSELK  = 1
USE_TASKLOOP     = 0

# SCOTCH
SCOTCH_HOME = $(EBROOTSCOTCH)
LIB_SCOTCH  = -L$(SCOTCH_HOME)/lib64 -lptesmumps -lptscotch -lptscotchparmetisv3 -lptscotcherr -lptscotcherrexit -lesmumps -lscotch -lscotchmetisv3 -lscotcherr -lscotcherrexit -lscotch_group
INC_SCOTCH  = -I$(SCOTCH_HOME)/include

# METIS
METIS_HOME = $(EBROOTPARMETIS)
LIB_METIS  = -L$(METIS_HOME)/lib64 -lmetis
INC_METIS  = -I$(METIS_HOME)/include

# PARMETIS
PARMETIS_HOME = $(EBROOTPARMETIS)
LIB_PARMETIS = -L$(PARMETIS_HOME)/lib64 -lparmetis
INC_PARMETIS = -I$(PARMETIS_HOME)/include

# MUMPS
MUMPS_HOME = $(MUMPS_ROOT_GNU)
LIB_MUMPS  = -L$(MUMPS_HOME)/lib -lzmumps -ldmumps -lmumps_common $(LIB_PARMETIS) $(LIB_METIS) $(LIB_SCOTCH)
INC_MUMPS  = $(MUMPS_HOME)/include $(INC_PARMETIS) $(INC_METIS) $(INC_SCOTCH)
#ORDLIB     = -L$(MUMPS_HOME)/PORD/lib -lpord $(LIB_PARMETIS) $(LIB_METIS) $(LIB_SCOTCH)


# STRUMPACK
STRUMPACK_HOME = $(STRUMPACK_ROOT_GNU)
STRUMPACKLIB   = -L$(STRUMPACK_HOME)/lib64 -lstrumpack $(LIB_PARMETIS) $(LIB_METIS) $(LIB_SCOTCH)
STRUMPACKINC   = -I$(STRUMPACK_HOME)/include $(INC_PARMETIS) $(INC_METIS) $(INC_SCOTCH)

# BOOST
BOOST_HOME = $(EBROOTBOOST)
LIB_BOOST  = -L$(BOOST_HOME)/lib64 -lboost_math_tr1
INC_BOOST  = -I$(BOOST_HOME)/boost/include 

# LIBFFTW
#FFTW_HOME = $(FFTW_ROOT)
FFTW_HOME = /opt/cray/pe/fftw/3.3.10.10/x86_milan
LIBFFTW   = -L$(FFTW_HOME)/lib -lfftw3f_mpi -lfftw3f_omp -lfftw3f -lfftw3_mpi -lfftw3_omp -lfftw3_threads -lfftw3
INC_FFTW  = -I$(FFTW_HOME)/include

# HDF5
HDF5INCLUDE = $(HDF5_ROOT)/include
HDF5LIB     = $(HDF5_ROOT)/lib/libhdf5hl_fortran_parallel.so \
              $(HDF5_ROOT)/lib/libhdf5hl_fortran.so \
              $(HDF5_ROOT)/lib/libhdf5_hl_parallel.so \
	      $(HDF5_ROOT)/lib/libhdf5_hl.so \
              $(HDF5_ROOT)/lib/libhdf5_parallel.so \
	      $(HDF5_ROOT)/lib/libhdf5_fortran.so \
              $(HDF5_ROOT)/lib/libhdf5.so $(EBROOTZLIB)/lib64/libz.a -ldl

#LIBLAPACK   = -L$(CRAY_LIBSCI_PREFIX)/lib -lsci_gnu_mpi_mp -lsci_gnu_mp
```

## File Systems
Information about data storage can be found in [LUMI storage](https://docs.lumi-supercomputer.eu/storage/). The $HOME directory can contain up to 20 GB of data, the $WORK directory can hold 50 GB. There are **no backups** of any storage systems of LUMI.

You can check the memory and file usage quotas of your projects with the following command:

```bash
lumi-workspaces
```

Copying files between different UNIX-like systems can be done with the scp command. This command, which stands for Secure Copy Protocol, allows you to transfer files between a local host and a remote host or between two remote hosts. The basic syntax of the scp command is the following (example):

```bash
scp -i ~/.ssh/id_rsa_lumi username@lumi.csc.fi:/projappl/project_number/user_folder/logfile.out .
```

## Batch System
- Basic commands:
```bash
sbatch <jobscript>                          # Submit a job script
sbatch -d afterok:<jobid> <jobscript>       # Submit a job script to start after <jobid> is finished

squeue -t running                           # List all currently running jobs
squeue -u <username>                        # List all own jobs

scancel <jobid>                             # Delete a queued or running batch job
```
- Information on partitions can be found [here](https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/partitions/).

- Example job script
```bash
#!/bin/bash
#SBATCH -J fixed
#SBATCH -A project_number
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --time=00:10:00
#SBATCH --mem=0
#SBATCH --exclusive

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load LUMI/25.03
module load partition/C
module load cpeGNU
module load zlib
module load wget
module load util-linux
module load systools
module load syslibs
module load nano
module load libreadline
module load libtirpc
module load libunistring
module load gettext
module load cURL
module load Vim
module load Tcl
module load ParMETIS
module load GLib
module load Boost
module load cray-hdf5-parallel
module load zlib
module load cray-libsci
module load cray-fftw
module load cray-python
module load SCOTCH
module load EasyBuild-user
module load MUMPS/5.6.1-cpeGNU-25.03-OpenMP

ulimit -s unlimited

srun ./jorek_model600 < input_jorek > logfile.out
```

## Accounting
Avaialble partitions can be found on [this website](https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/partitions/). 
Information on billing is [here](https://docs.lumi-supercomputer.eu/runjobs/lumi_env/billing/).

It is recommended to regularly check one's allocation and usage using the following command: 
```bash
lumi-allocations
```
