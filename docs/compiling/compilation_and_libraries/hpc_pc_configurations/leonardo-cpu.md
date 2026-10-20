---
title: "Leonardo-CPU"
nav_order: 1
parent: "HPC & PC Configurations"
grand_parent: "Compilation and Libraries"
layout: default
render_with_liquid: false
---

# Running JOREK on Leonardo CPU partition (DCGP)

## Getting Access

- Account Creation via https://userdb.hpc.cineca.it/ - Create New User
- Then fill out the HPC-Related information and upload your scanned ID/passport
- Contact the respective PI to add you to the project

## README

- [LEONARDO User Guide](https://wiki.u-gov.it/confluence/display/SCAIUS/LEONARDO+User+Guide)
- [Setup for Eurofusion users](https://wiki.u-gov.it/confluence/display/SCAIUS/EUROfusion+users%3A+Marconi+and+Leonardo+environments)
- Slides of the webinar offered by CINECA on the 18th of February, 2025 [slides_qa.pdf](slides_qa.pdf)
- E-Mail Address of CINECA Helpdesk: superc@cineca.it

## Login

All the login nodes (4 icelake, no-gpu) have an identical environment and can be reached with SSH (Secure Shell) protocol using the "collective" hostname:

```bash
ssh -Y <user>@login.leonardo.cineca.it
```
which establishes a connection to one of the available login nodes. You can also indicate explicitly the login nodes:

```bash
login01-ext.leonardo.cineca.it
login02-ext.leonardo.cineca.it 
login05-ext.leonardo.cineca.it
login07-ext.leonardo.cineca.it
```
**Follow the instructions at [this link](https://wiki.u-gov.it/confluence/display/SCAIUS/FAQ#FAQ-Ikeepreceivingtheerrormessage%22WARNING:REMOTEHOSTIDENTIFICATIONHASCHANGED!%22evenifImodifyknown_hostfile) if you get the error message “WARNING: REMOTE HOST IDENTIFICATION HAS CHANGED!”**

## Software Available
As usual module avail/load/list/show allows to access the available software

To use the file manager midnight commander, put these two lines in your .bashrc:

```bash
export PATH=/leonardo/pub/userexternal/dbonfigl/mc/bin:$PATH
export LD_LIBRARY_PATH=/leonardo/pub/userexternal/dbonfigl/mc/shares/
mc:$LD_LIBRARY_PATH
```

## Compiling JOREK (As of 2025.02: Under test)
Below are some examples of the .bashrc and Makefile.inc setup.

[Using GNU + OPENMPI](leonardo_cpu_with_gnu_openmpi.md) (Most stable setup so far & libraries publicly accessible)

[Using Intel + OneAPI](leonardo_cpu_with_intel_oneapi.md) (Under tests & requiring access to the FUAL8_MHD project)

[Using Intel + OneAPI + Raven Libraries](leonardo_cpu_with_raven_libraries.md) (Under tests & using the copied libraries from Raven)

## File Systems
Most important are $HOME (user specific) and $WORK (project specific), see the LEONARDO User Guide for details. BE AWARE THAT $HOME IS THE ONLY FILE SYSTEM THAT GETS BACKED UP! (suspended for now)

One can check quota with:

bash
cindata

## Batch System
- Basic commands:
```bash
sbatch <jobscript>                          # Submit a job script
sbatch -d afterok:<jobid> <jobscript>       # Submit a job script to start after <jobid> is finished

squeue -t running                           # List all currently running jobs
squeue -u <username>                        # List all own jobs

scancel <jobid>                             # Delete a queued or running batch job
```

- Example job script
Example job script for 2*8 MPI tasks on 2 different nodes with 14 OpenMP threads each. Note: You need to add _0 to your project name, e.g., FUA38_MHD_0

```bash
#!/bin/bash
  
#SBATCH -J test_jobscript
#SBATCH -A your_project
#SBATCH -p dcgp_fua_prod
#SBATCH --qos=<qos_name>

#SBATCH --nodes=2
#SBATCH --cores-per-socket=56                  # 56*2 sockets
#SBATCH --ntasks-per-node=8                    # number of tasks per node
#SBATCH --cpus-per-task=14                     # number of cores per task max 32 cores per node

#SBATCH --mem=494000                           # memory per node
#SBATCH --gres=tmpfs:10g

#SBATCH --time=24:00:00
#SBATCH --output=job.out
#SBATCH --error=job.err


### Request e-mail notification
#SBATCH --mail-type=FAIL   # can be BEGIN, END, TIME_LIMIT, FAIL, ALL, NONE
#SBATCH --mail-user=your_email_address

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE='50M'  ## Optional, depending on the size of the job
export I_MPI_THREAD_LEVEL_DEFAULT=funneled #multiple

srun --cpu-bind=cores --cpus-per-task=$SLURM_CPUS_PER_TASK ./jorek_model600 < jorek_namelist > logfile
```
## Accounting
To check the used resources on Leonardo DCGP, use:

```bash
saldo -b --dcgp
```
To check more details (containing the usage of all the users) one can use:

```bash
saldo -ra PROJECT_NAME --dcgp
```
Available partitions (see also [here](https://wiki.u-gov.it/confluence/display/SCAIUS/EUROfusion+users%3A+Marconi+and+Leonardo+environments)):

- Serial partition lrd_all_serial. 2 dedicated login-type nodes. Max 4 physical core per job, max walltime 4 hours. Jobs are free (compute resources are not accounted). Good for running jorek2_postproc for example.

- Production partition dcgp_fua_prod. 258 compute nodes dynamically allocated. Max 16 nodes per job, max walltime 24 hours. For big production runs (min 17 full nodes, max 64 nodes) you can specify the QOS dcgp_qos_fuabprod.

- Debug partition dcgp_fua_dbg. 2 compute nodes dynamically allocated. Max 2 nodes, max walltime 10 min.

- Low priority queue on the production partition. Max 16 nodes, max walltime 8 hours. Jobs are free (compute resources are not accounted). This can be used by specifying the FUA38_LOWPRIO_0 account (ask CINECA to be added to this 'project') and the QOS qos_fualowprio. Good to save project resources. Typically jobs start in the weekend after their submission.
