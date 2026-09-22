---
title: "Compile without Intel MKL"
nav_order: 7
parent: "Compilation and Libraries"
grand_parent: "Compiling and Running"
layout: default
render_with_liquid: false
---

# Compiling JOREK and its dependencies without access to Intel-MKL libraries

It is not entirely straightforward to compile JOREK with GNU compilers when dependencies like ScaLAPACK are not provided by the support team. Here is one possible approach.

We assume:

- you have access to a relatively recent GNU compiler, e.g. 7.5
- you have access to CMake
- you may want to compile older versions of SCOTCH and PASTIX, for which the latest version of OpenMPI (4.0) does not work
- you are working in the directory:

```text
/home/user/lib/
```

## HDF5

First, you need to compile the HDF5 libraries.

You will need to download the code from:

- <https://www.hdfgroup.org/downloads/hdf5/source-code/>

Do not try to use the code from the HDF Git repository; this does not work properly for the Fortran libraries.

Once you have unpacked the HDF5 tarball, do the following:

```bash
cd ./hdf5/
mkdir install
./configure --prefix=/home/user/lib/hdf5/install/ --enable-fortran
make -j8
make -j8 install
```

## OpenMPI

Next, download the OpenMPI sources from:

- <https://www.open-mpi.org/software/ompi/v3.0/>

Unpack the tarball, then run:

```bash
cd ./openmpi
mkdir install
./configure --prefix=/home/user/lib/openmpi/install/
make -j8
make -j8 install
export LD_LIBRARY_PATH=/home/user/lib/openmpi/install/lib:$LD_LIBRARY_PATH
export PATH=/home/user/lib/openmpi/install/bin:$PATH
```

## OpenBLAS

Next, build the OpenBLAS library:

```bash
git clone https://github.com/xianyi/OpenBLAS.git
export CC=mpicc
export FC=mpifort
export CXX=mpicxx
cd OpenBLAS
mkdir build install
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=/home/user/lib/OpenBLAS/install
make -j8
make -j8 install
```

## LAPACK

Next, build the LAPACK library:

```bash
git clone https://github.com/Reference-LAPACK/lapack.git
mkdir build install
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=/home/user/lib/lapack/install
make -j8
make -j8 install
```

## ScaLAPACK

Next, build the ScaLAPACK library:

```bash
git clone https://github.com/Reference-ScaLAPACK/scalapack.git
mkdir build install
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=/home/user/lib/scalapack/install
make -j8
make -j8 install
```

## PASTIX-5

SCOTCH compilation is straightforward once you have set `PATH` and `LD_LIBRARY_PATH` to include the OpenMPI compilers.

For the PASTIX solver, you will need to have the following inside `config.in`:

```makefile
BLAS_HOME=/home/user/lib/scalapack/install/lib
BLASLIB = -L$(BLAS_HOME) -lscalapack
```

Note: this may not even be required, as long as linking is done properly in JOREK.

## JOREK

For JOREK, you will need to use these lines instead of the usual MKL settings in your `Makefile.inc`:

```makefile
SCA_HOME  = /home/user/lib/scalapack/install
BLA_HOME  = /home/user/lib/OpenBLAS/install
LAP_HOME  = /home/user/lib/lapack_GNU/install
LIBLAPACK = -L$(SCA_HOME)/lib -lscalapack -L$(LAP_HOME)/lib64 -llapack -L$(BLA_HOME)/lib64 -lopenblas
SCALAP    = -L$(SCA_HOME)/lib -lscalapack -L$(LAP_HOME)/lib64 -llapack -L$(BLA_HOME)/lib64 -lopenblas
```
