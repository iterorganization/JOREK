---
title: "ITER Cluster"
nav_order: 13
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---

# ITER cluster Generalities

See:

- <https://confluence.iter.org/pages/viewpage.action?pageId=267026649> for access
- <https://confluence.iter.org/display/IMP/ITER+Computing+Cluster> for more detailed information

# Compiling

You need to prepare a specific environment by loading appropriate modules and exporting some environment variables.

You can prepare a dedicated `.sh` file, such as the one below:

```bash
#!/bin/bash

module purge

module load IMAS/3.41.0-2024.04-intel-2023b
module load VSCode/1.88.1
module load SCOTCH/7.0.4-iimpi-2023b
module load MUMPS/5.7.2-intel-2023b-parmetis
module load HDF5/1.14.3-iimpi-2023b
module load PaStiX/5.2.3-intel-2023b
module load FFTW/3.3.10-GCC-13.2.0
module load STRUMPACK

export PASTIX_HOME=$EBROOTPASTIX
export MUMPS_HOME=$EBROOTMUMPS
export METIS_HOME=$EBROOTPARMETIS
export SCOTCH_HOME=$EBROOTSCOTCH
export MKL_HOME=$EBROOTIMKL
export HDF5_HOME=$EBROOTHDF5
export LANG=C
export JOREK_HOST=iter-hpc
export compilethreads=8
export MAKEFLAGS="-j$compilethreads"
export PRERUN="export OMP_NUM_THREADS=8"
export MPIRUN="mpirun -np "
export BATCHCOMMAND="qsub"
export CXXFLAGS=-O0 # problem with stdio library on ITER http://gcc.1065356.n8.nabble.com/g-4-8-fails-with-Ox-option-td953876.html
```

You can either load it before compilation via:

```bash
source env2024.sh
```

or load it directly at each startup by updating your `.bashrc` file, for example if your `env2024.sh` is located in `~/public`:

```bash
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

if [ -f $HOME/.setup.sh ]; then
        . $HOME/.setup.sh
fi

# User specific aliases and functions

export PATH=.:$PATH
export CC=mpiicx
export FC=mpiifort
export F90=ifort
export plots=~/scripts_hub/scripts_javier/plot_scripts/

alias l.='ls -d .* --color=auto'
alias ll='ls -l --color=auto'
alias ls='ls --color=auto'
alias sqe='squeue -u $USER -o "%.8i %.2t %.10M %.6D %Z %R"'
alias q10='squeue -p gen10'

if [ -z "$PS1" ]; then
    # prompt var is not set, so this is *not* an interactive shell
    return
fi

# Decide OS and load correct environment file
source /etc/os-release
if [[ "$ID" == "rhel" ]]; then
    echo "RedHat"
    source ~/public/env2024.sh
elif [[ "$ID" == "centos" && "$VERSION_ID" == "8" ]]; then
#    echo 2
    echo "CentOS 8"
    source ~/public/env2020.sh
elif [[ "$ID" == "centos" && "$VERSION_ID" == "7" ]]; then
    echo "CentOS 7"
    source ~/public/env.sh
else
    echo "Unknown OS or version"
fi


export MPLBACKEND=TkAgg
export PATH="~/.local/bin:$PATH"

ulimit -s unlimited
```

If you already started your session, the best thing to do may be to re-run the directives in both files. From your shell in the `$HOME` directory, run:

```bash
source .bashrc
source ./public/env2024.sh
```

Then go to the JOREK code directory:

```bash
cd <jorek_directory>
```

and start the compilation from scratch:

```bash
make clean
make cleanall
make
```

# Running

It is recommended to run via the SLURM scheduler.

You have to choose the partition where to run your program, which may be `gen10` or `gen11` for production runs. For debug runs, you may use the dedicated partitions `gen10_debug` or `gen11_debug`.

A typical batch script could look like this:

```bash
#!/bin/bash -l
#SBATCH --time=01:00:00
#SBATCH --job-name=JOREK_TM
#SBATCH --exclusive
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=14
#SBATCH --nodes=2
#SBATCH --partition=gen11_debug
#SBATCH --error myjob.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SLURM_CPU_BIND="rank"  
#export I_MPI_PMI_VALUE_LENGTH_MAX=512
export I_MPI_THREAD_LEVEL_DEFAULT=multiple
srun ./jorek_model600 < input > output
```

Here, `input` is your configuration file and `output` is the text file where you want to save the runtime log information.

For proper resource allocation, see the ITER-SDCC guide.

Remember to include any required extra files in your working directory, for example equilibrium profiles or the wall response file, usually named `starwall-response.dat`.

# Compiling CARIDDI

If CARIDDI does not compile with the environment above, do the following:

```bash
module purge
module load intel/2023b
module load ScaLAPACK/2.2.0-gompi-2023b-fb
make clean
make
```
