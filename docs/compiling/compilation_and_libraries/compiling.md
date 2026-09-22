---
title: "Compile"
nav_order: 1
parent: "Compilation and Libraries"
grand_parent: "Compiling and Running"
layout: default
render_with_liquid: false
---

# Compiling JOREK

**This page explains the libraries, build options, and other components required to compile JOREK. Before compiling any dependencies yourself, check whether prepared libraries are already available as modules on your HPC system.**

**Ask colleagues who work on the relevant machine which modules and configurations they recommend, and consult the [HPC & PC Configurations](hpc_pc_configurations.md) section for system-specific instructions.**

- **Check the code out from the git repository** as explained [here](../get_jorek.md).
- **Compile JOREK** as explained at the end of this page.

## Specialized Guides

- [Compiling on macOS](/JOREK/compiling/compilation_and_libraries/hpc_pc_configurations/macos.html)
- [Compiling MUMPS (+METIS)](/JOREK/compiling/compilation_and_libraries/compiling_mumps.html)
- [Compiling with 64-bit integers](/JOREK/compiling/compilation_and_libraries/compiling_int64.html)
- [Compiling all dependencies without Intel MKL](/JOREK/compiling/compilation_and_libraries/compiling_wo_mkl.html)

## SCOTCH

*Scotch is a software package for graph and mesh/hypergraph partitioning, graph clustering, and sparse matrix ordering.*

- [Download the source code here](https://files.inria.fr/pastix/releases/) or use [the version with `libesmumps` included](http://gforge.inria.fr/frs/download.php/file/28978/scotch_5.1.12b_esmumps.tar.gz)
- Currently version 5.1.12b is recommended; recently also version 6.0.4 worked well but is not yet suitable if several MPI tasks are required per toroidal harmonic (feature to be added with 6.1.x)
- Unpack with

```bash
tar -xvzf scotch_5.1.12b.tar.gz
```

- Create the configuration file `Makefile.inc` in the `src/` folder
  - Generic examples can be found in `src/Make.inc/`
  - [Suitable example for use with JOREK (Scotch 5.1.12b)](exp_scotch.md)
- Compile the library:

```bash
make clean
make -j 8 scotch
make -j 8 ptscotch
```

## PaStiX (before version 6.x)

*PaStiX is a high performance parallel solver for very large sparse linear systems based on direct methods. It provides also an adaptive blockwise iLU(k) factorization that can be used as a parallel preconditioner.*

- [Download the source code here](https://files.inria.fr/pastix/releases/)
- **Note that there is a bug in PaStiX revisions 3519 to 4030 leading to large memory consumption!**
- Unpack with:

```bash
tar -xvjf pastix_5.2.3.tar.bz2
```

- Create the configuration file `config.in` in the `src/` folder
  - Generic examples can be found in `src/config/`
  - [Suitable example for use with JOREK (PaStiX revision 4492)](exp_pastix.md)
  - Do not forget to uncomment the memory usage statistics, located under Options in `config.in`
- Compile the library:

```bash
make clean
make -j 8
make -j 8 install
```

- It is possible to validate the PaStiX installation after it has been compiled:

```bash
# simplest example
make examples
mpirun -np <Nprocs> ./example/bin/simple -rsa matrix/small.rsa -t <nthreads> -v 4

# Calling PaStiX with separated steps (as in Jorek)
mpirun -np <Nprocs> ./example/bin/step-by-step -rsa matrix/small.rsa -t <nthreads> -v 4

# Calling PaStiX with several communicators (as in Jorek)
mpirun -np <Nprocs> ./example/bin/multi-comm -rsa matrix/small.rsa -t <nthreads> -v 4
```

## PaStiX 6

The new version of PaStiX (with new features such as block-low-rank compression).

- [Download the source code here](https://gforge.inria.fr/frs/?group_id=186) (release versions)
- [Download latest develop version here](https://gitlab.inria.fr/solverstack/pastix) (see Readme for git download instructions)
- [Download the tar.gz of the version used for the benchmarks and BLR test from 2019/07](assets/compiling/pastix6.tar.gz)
- Unless downloaded from git, unpack with:

```bash
tar -xvzf pastix-6.0.2.tar.gz
```

- Go to the PaStiX 6 base directory, create a `build` subdirectory and build with CMake. Several options can be set, most of which will be set automatically but you might need to specify your ordering library (Scotch or Metis) and the compiler used to compile it, as well as the installation prefix (where PaStiX 6 should be installed). In this case, the command line argument is as follows (replacing `...` with the appropriate paths):

```bash
mkdir build
cd build
cmake .. -DSCOTCH_DIR=... -DCMAKE_C_COMPILER=... -DCMAKE_INSTALL_PREFIX=...
```

- Just like previous versions, PaStiX 6 can now be installed from the PaStiX home directory with

```bash
make clean
make -j 8
make -j 8 install
```

## FFTW

*JOREK also uses FFTW for fast Fourier transformations. This library is usually available on machines via `module load fftw`.* See also [here](preprocessor.md) to include FFTW.

## Compile JOREK

We use `Make` to compile JOREK, with dependency generation via `util/makedepend`.
See [Makefile](makefile.md) for more information.

- **Create the JOREK configuration file `trunk/Makefile.inc` for your system**, examples are available in `trunk/Make.inc/`.
- Add the following lines to link to the libraries:

```make
# SCOTCH Library: Automatically included via the pastix-conf script calls below

# PaStiX Library:
USE_PASTIX        = 1
USE_PASTIX_MURGE  = 1  # Different way for matrix assembly, the choice is up to you
PASTIX_HOME       = <path-to-pastix>/install/
LIB_PASTIX        = `$(PASTIX_HOME)/pastix-conf --libs`
LIB_PASTIX_MURGE  = `$(PASTIX_HOME)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS   = `$(PASTIX_HOME)/pastix-conf --blas`
INC_PASTIX        = `$(PASTIX_HOME)/pastix-conf --incs`

# FFTW Library: FFTW_HOME might be set by "module load fftw", otherwise do so manually
LIBFFTW           = $(FFTW_HOME)/lib/libfftw3.a
```

- For PaStiX after revision 3945, additionally the following is required (see also [preprocessor flags](preprocessor.md)):

```make
FFLAGS := $(FFLAGS) -DFUNNELED -DWORLDWAR2
```

- In case of **PaStiX 6.x**, instead of the parameters shown above, `USE_PASTIX6=1`, `PASTIX6_HOME`, `LIB_PASTIX6`, `LIB_PASTIX6_BLAS`, and `INC_PASTIX6` need to be set correctly. In that case, `USE_PASTIX` must be 0.

- Set the [hard-coded parameters](hard-coded_parameters.md) correctly. For example:

```bash
./util/config.sh model=303 n_tor=3 n_period=8 n_plane=8   # set hard-coded parameters
./util/config.sh                                          # check hard-coded parameters
```

- Compile the main program and/or diagnostics tools:

```bash
make                       # compile main program jorek_modelxxx
make jorek2vtk
make jorek2_postproc
...
```

- Parallel compilation is possible, e.g. with 8 threads:

```bash
make -j 8
```

### Special Compile Options for Single Source Files

- It is now possible to specify special compile options (added to the normal `FFLAGS`) for individual source files. An example for the file `models/model333/element_matrix_fft.f90` (to put in `Makefile.inc`):

```make
.obj/element_matrix_fft.o: FFLAGS+=-O3
```

### Remarks

- **When switching the physics model or hard-coded parameters**, you should always call a `make clean` before recompiling. When dependencies between files differ for different models (preprocessor directives), this is not sufficient and a `make cleanall` is required. `./util/config.sh` has been adapted to do this automatically.
- **Compiler Optimization:** (Preliminary -- will be updated)
  - Tests in Garching on our local Linux Cluster and on our supercomputer Hydra showed that the best optimization options seem to be simply `-fast -DUSE_BLOCK` without anything else like `-O2 -xSSE4.2` etc. Both systems are equipped with Intel Ivy Bridge processors very similar to Helios.
