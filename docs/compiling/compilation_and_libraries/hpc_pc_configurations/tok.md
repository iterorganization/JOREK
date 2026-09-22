---
title: "TOK Cluster"
nav_order: 10
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---

# Run JOREK on the TOK Linux cluster (installed 2025)

## Basic information

- [Information about the TOK cluster](https://hpc.wiki.ipp.mpg.de/tok_batch_cluster_in_garching)

**Two login nodes are available** (only for compilation and submission of jobs):

```bash
ssh tok[11-12].ipp.mpg.de
```

**Pool of interactive nodes** (to run small tests):

For running small interactive tests, the TOK-I cluster `toki01...toki12.ipp.mpg.de` is available with slightly older hardware.

**Compute nodes:**

- 72 nodes
- 2× AMD EPYC 9754 CPUs per node, each with 128 cores and 2.25 GHz base clock speed
- Memory per node: 768 GB
- Operating system: SLES15

## Modules required

Basic set for compilation:

```bash
module load intel/2025.1 impi/2021.15 mkl/2025.1 fftw-mpi/3.3.10 hdf5-mpi/1.14.1
```

Additional module for running interactively (without SLURM):

```bash
module load impi-interactive
```

Additional modules for post-processing analysis, visualization, data, and code sharing:

```bash
module load visit
module load paraview
module load git
module load gnuplot
module load texlive
module load imagemagick
module load datashare
```

You can either add these modules to your `~/.bashrc` file (not recommended, as clashes with other environments can be hard to debug), or use the `module save/restore` functionality.

## Compiling JOREK and Libraries

> **Warning:** There seems to be an issue when JOREK is run on clusters at IPP (TOK-I, TOK, DRACO, and COBRA) with more than one MPI task per harmonic. The run hangs at the first PASTIX call during the very first time-step. In such a situation, the combination of `Pastix_5.2.3` and `Scotch_5.1.12b` has been seen to fix this issue.

- [Installation/compilation of libraries and JOREK](../compiling.md) works in the normal way.
- You can simply use the libraries already compiled by Ihor Holod (you may need to ask Matthias Hoelzl for access) via the following `Makefile.inc` file:

```makefile
# JOREK configuration file for compilation on the TOK cluster at IPP Garching
# from: https://jorek.eu/wiki/doku.php?id=garching-tok-cluster
 
MODEL             = model199
 
FC                = mpiifx
CC                = mpiicx
CXX               = mpiicpx
 
FFLAGS += -O2 -march=native
 
COMPILER_FAMILY   = intel
 
USE_PASTIX        = 0
USE_STRUMPACK     = 1
USE_MUMPS         = 0
USE_PASTIX_MURGE  = 0
USE_BLOCK         = 1
USE_HDF5          = 1
USE_FFTW          = 1
USE_MKL           = 0
 
# folder owned by Matthias Hoelzl - ask for access if needed
LIBDIR = /tokp/work/mhoelzl/bin-of-ihol/intel/
 
ifeq (1, $(USE_MUMPS))
 MKLLIB= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 \
         -Wl,-rpath,${MKLROOT}/lib/intel64
 MKLINC = -I${MKL_HOME}/include
 
 MUMPSDIR = $(LIBDIR)/MUMPS_5.2.1
 LIB_MUMPS = -L$(MUMPSDIR)/lib -ldmumps -lmumps_common -lpord
 INC_MUMPS = $(MUMPSDIR)/include
 
 LIB_MUMPS += $(MKLLIB)
 INC_MUMPS += $(MKLINC)
endif
 
 
ifeq (1, $(USE_PASTIX))
 MKLLIB = -L${MKL_HOME}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
         -Wl,-rpath,${MKLROOT}/lib/intel64
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
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 \
         -Wl,-rpath,${MKLROOT}/lib/intel64
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
LIBFFTW           = -L$(FFTW_HOME)/lib -lfftw3 -Wl,-rpath,$(FFTW_HOME)/lib
INC_FFTW          = -I$(FFTW_HOME)/include
 
# --- HDF5 Library 
HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz  -Wl,-rpath,$(HDF5_HOME)/lib
```

## File systems

- Best is to use `/tokp/work/<user>` (your home directory) for your simulations
- Your quota usage should be 2 TB in `/tokp/work/<user>` by default, and this data will automatically be backed up
- **Be careful when using scratch directories, as there are no automatic backups**

## Batch submission

- Submit a job:

```bash
sbatch jobscript
```

- Check job status:

```bash
squeue -u <username>
```

- Cancel a job:

```bash
scancel <jobid>
```

- List available queues (partitions):

```bash
sinfo
```

### Example job scripts

#### Serial job

```bash
#!/bin/bash -l
#SBATCH -J jorek_testing      #Job name
#SBATCH -o ./%x.%j.out        #stdout (%x=jobname, %j=jobid)
#SBATCH -e ./%x.%j.err        #stderr (%x=jobname, %j=jobid)
#SBATCH -D ./                 #Initial working directory
#SBATCH --partition=s.tok     #Queue/Partition
#SBATCH --qos=s.tok.standard  #Quality of Service: s.tok.short, s.tok.standard, s.tok.long, tok.debug
#SBATCH --nodes=1             #Total number of nodes  
#SBATCH --ntasks-per-node=1   #MPI tasks per node
#SBATCH --cpus-per-task=1     #CPUs per task for OpenMP
#SBATCH --mem-per-cpu 1GB     #Set memory requirement (default: 1 GB, max: 5GB)
#SBATCH --time=12:00:00       #Wall clock limit
##
#SBATCH --mail-type=end       #Send mail, e.g. for begin/end/fail/none
#SBATCH --mail-user=YOUR_USER_ID@rzg.mpg.de  #Mail address

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
# For pinning threads correctly:
export OMP_PLACES=cores 

# Run the program:
srun ./jorek_model199 < ./intear
```

#### OpenMP job

```bash
#!/bin/bash -l
#SBATCH -J jorek_testing_omp      #Job name
#SBATCH -o ./%x.%j.out        #stdout (%x=jobname, %j=jobid)
#SBATCH -e ./%x.%j.err        #stderr (%x=jobname, %j=jobid)
#SBATCH -D ./                 #Initial working directory
#SBATCH --partition=p.tok.openmp     #Queue/Partition
#SBATCH --qos=p.tok.openmp.48h       #Quality of Service: p.tok.openmp.2h, p.tok.openmp.24h, p.tok.openmp.48h, tok.debug
#SBATCH --nodes=1             #Total number of nodes  
#SBATCH --ntasks-per-node=1   #MPI tasks per node
#SBATCH --cpus-per-task=4     #CPUs per task for OpenMP
#SBATCH --mem-per-cpu 1GB     #Set memory requirement (default: 1 GB, max: 5GB)
#SBATCH --time=12:00:00       #Wall clock limit
##
#SBATCH --mail-type=end       #Send mail, e.g. for begin/end/fail/none
#SBATCH --mail-user=YOUR_USER_ID@rzg.mpg.de  #Mail address

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
# For pinning threads correctly:
export OMP_PLACES=cores 

# Run the program:
srun ./jorek_model199 < ./intear | tee logfile
```

#### OpenMP + MPI job

```bash
#!/bin/bash -l
#SBATCH -J jorek_testing           #Job name
#SBATCH -o ./%x.%j.out             #stdout (%x=jobname, %j=jobid)
#SBATCH -e ./%x.%j.err             #stderr (%x=jobname, %j=jobid)
#SBATCH -D ./                      #Initial working directory
#SBATCH --partition=p.tok          #Queue/Partition
#SBATCH --qos=p.tok.48h            #Quality of Service: p.tok.2h, p.tok.48, tok.debug
#SBATCH --nodes=2                  #Total number of nodes  
#SBATCH --ntasks-per-node=1        #MPI tasks per node
#SBATCH --cpus-per-task=2          #CPUs per task for OpenMP
#SBATCH --mem 90GB                 #Set mem./node requirement (default: 63000MB, max: 190GB)
#SBATCH --time=12:00:00            #Wall clock limit
##
#SBATCH --mail-type=end            #Send mail, e.g. for begin/end/fail/none
#SBATCH --mail-user=YOUR_USER_ID@rzg.mpg.de  #Mail address

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
# For pinning threads correctly:
export OMP_PLACES=cores 

# Run the program:
srun ./jorek_model199 < ./intear | tee logfile
```

More example job scripts for the TOK cluster can be found [here](https://hpc.wiki.ipp.mpg.de/tok_batch_cluster_in_garching).

## VNC Viewer

The setup of a Virtual Network Connection server is outlined [here](https://www.mpcdf.mpg.de/services/computing/software/vnc-at-the-mpcdf).

In addition to these commands, it has been observed that visualization software can crash on the new systems, seemingly due to errors related to required fonts. It is therefore recommended to add the following flag when setting up the VNC server:

```bash
vncserver -fp "built-ins"
```

## Might be useful

### Depending jobs

- Add the following to your `~/.bashrc` file (it might not exist, so create it if necessary):

```bash
function qsubchain() {
if [ $# -lt 2 ]; then
  echo "Usage: qsubchain jobscript1 jobscript2 [...]"
  return
fi

thejobids=""
depend=""
while [ $# -gt 0 ]; do
  thejobid=`sbatch $depend $1 | sed -e 's/Submitted batch job *//'`
  thejobids="$thejobids $thejobid"
  depend="--dependency=afterok:$thejobid"
  shift
done
echo "JobIDs: $thejobids"
unset depend thejobid thejobids
}
```

- You can then submit two or more job scripts at once, where the second job only starts when the first one has finished, and so on:

```bash
$ qsubchain jobscript1 jobscript2
JobIDs:  43693 43694
```
