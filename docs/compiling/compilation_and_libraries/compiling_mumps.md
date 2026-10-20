---
title: "Compile MUMPS"
nav_order: 6
parent: "Compilation and Libraries"
grand_parent: "Compiling and Running"
layout: default
render_with_liquid: false
---

# Compiling MUMPS

This document is a short guide to installing MUMPS and its dependencies (parallel version).

## Dependencies

- BLAS
- MPI
- BLACS
- ScaLAPACK

You will probably also want to install SCOTCH or METIS for graph partitioning, especially if using Block-Low-Rank compression in MUMPS.

## BLAS

Usually BLAS can be installed through your package manager or will already be available as a module.

## MPI

MPI comes in different versions (`mpich`, `mpich2`, `openmpi`) and is usually already available on your system, or can be installed through your package manager.

## BLACS

> The BLACS (Basic Linear Algebra Communication Subprograms) project is an ongoing investigation whose purpose is to create a linear algebra oriented message passing interface that may be implemented efficiently and uniformly across a large range of distributed memory platforms.

See [netlib](http://www.netlib.org/blacs/) for more information.

If BLACS is not available on your system already (e.g. through MKL):

- Download the MPI version [here](http://www.netlib.org/blacs/mpiblacs.tgz)
- Run:
  ```bash
  tar zxvf mpiblacs.tgz
  ```
- Run:
  ```bash
  cp BMAKES/Bmake.MPI-LINUX Bmake.inc
  ```
- Edit `Bmake.inc` to set `BTOPdir` correctly
- See [this page](https://www.open-mpi.org/faq/?category=mpi-apps#blacs) to build with OpenMPI instead of MPICH
- Build the MPI version:
  ```bash
  make mpi
  ```
- Libraries will be placed in `BLACS/LIB`

## ScaLAPACK

> ScaLAPACK is a library of high-performance linear algebra routines for parallel distributed memory machines. ScaLAPACK solves dense and banded linear systems, least squares problems, eigenvalue problems, and singular value problems.

If ScaLAPACK is not available on your system already (e.g. through MKL):

- Download [scalapack-2.0.2.tgz](http://www.netlib.org/scalapack/scalapack-2.0.2.tgz) or a newer version
- Run:
  ```bash
  tar zxvf scalapack-2.0.2.tgz
  ```
- Ignore the BLACS directory included in the ScaLAPACK install
- Run:
  ```bash
  cp SLmake.inc.example SLmake.inc
  ```
- See [this page](https://www.open-mpi.org/faq/?category=mpi-apps#scalapack) for instructions on editing `SLmake.inc`
- Compile ScaLAPACK:
  ```bash
  make
  ```
- The output file will be `libscalapack.a` in the current directory

## METIS

> METIS is a set of serial programs for partitioning graphs, partitioning finite element meshes, and producing fill reducing orderings for sparse matrices.

- Download [METIS 5.1.0](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz) (serial version) or [ParMETIS 4.0.3](http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz) (parallel version)
- See the [METIS website](http://glaros.dtc.umn.edu/gkhome/views/metis) for installation instructions, documentation, and possibly newer versions
- Run:
  ```bash
  tar xvzf parmetis-4.0.3.tar.gz
  ```
- Configure:
  ```bash
  make config prefix=PARMETIS_ROOT
  ```
  replacing `PARMETIS_ROOT` accordingly
- Compile:
  ```bash
  make
  ```
- Install:
  ```bash
  make install
  ```
- Note that if you downloaded ParMETIS but also need METIS, it is placed in the `metis` subdirectory. Go to this subdirectory and simply repeat the steps above (with the root directory changed accordingly).

## Compiling MUMPS itself

> MUMPS is a parallel sparse direct solver.

It can optionally be used in the JOREK solver and is used in the particle routines for impurity transport.

- [Download the source code here](http://mumps.enseeiht.fr/MUMPS_5.2.0.tar.gz)
- Unpack with:
  ```bash
  tar zxvf MUMPS_5.2.0.tar.gz
  ```
- Copy the correct `Makefile.inc` from the `Make.inc` directory for your architecture. These Makefiles also contain instructions on adapting them to your needs. The `INSTALL` file contains further installation instructions.
- Edit `Makefile.inc` to contain the paths to BLACS and ScaLAPACK, for example:

```makefile
BLACS_DIR = /home/daan/source/BLACS/LIB
SCALAP  = -llapack -L$(BLACS_DIR) /home/daan/source/scalapack-2.0.2/libscalapack.a \
          $(BLACS_DIR)/blacs_MPI-LINUX-0.a \
          $(BLACS_DIR)/blacsF77init_MPI-LINUX-0.a \
          $(BLACS_DIR)/blacsCinit_MPI-LINUX-0.a
#SCALAP  = -lscalapack -lblacs  -lblacsF77init -lblacsCinit
```

- To use SCOTCH ordering, set `LSCOTCHDIR` to scotch with `esmumps`
  - For SCOTCH compilation, see the [main compilation guide](/JOREK/compiling/compilation_and_libraries/compiling.html)
  - See [this page](http://mumps.enseeiht.fr/index.php?page=faq#2) for information on SCOTCH ordering
- To use METIS ordering, set `LMETISDIR` and `IMETIS` to the appropriate paths and uncomment `LMETIS`
  - For METIS compilation, see above
- Compile the library in all precisions:
  ```bash
  make all
  ```
- In the `lib/` folder you will find `cmumps`, `dmumps`, `mumps_common`, and `pord` libraries
