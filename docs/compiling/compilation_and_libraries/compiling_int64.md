---
title: "Compile with 64-bit Integers"
nav_order: 5
parent: "Compilation and Libraries"
grand_parent: "Compiling and Running"
layout: default
render_with_liquid: false
---

# Compiling JOREK and solvers using long integers

Compiling the code with long integers allows larger matrices to be used.

Typically, a sparse solver would use a distributed matrix, and so by increasing the number of MPI processes, you reduce the matrix indexing size on each process.

However, there are two issues with this:

1. By increasing MPI processes, you increase overall memory usage.
2. Some solvers, like PASTIX, do not have a distributed matrix option.

STRUMPACK has a distributed matrix option, which is used in the Block-Jacobi harmonic matrices, but not in `solve_strumpack_all.f90` at the moment.

In any case, being able to use only 1 MPI process per harmonic, and all other cores as OMP, allows reducing the memory requirements. On memory-intensive architectures like Intel Optane, SGI UV-2000, or HPE SuperdomeFlex, this is crucial.

Here are some outlines on how to compile the solvers and JOREK with long integers.

## METIS & PARMETIS

METIS and PARMETIS simply require changing the integer size in the include files:

- `include/metis.h`
- `include/parmetis.h`

The CMake command for METIS on Marconi is (using Intel compilers):

```bash
cmake ../ -DGKLIB_PATH="../GKlib" -DCMAKE_INSTALL_PREFIX="../install" -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_Fortran_COMPILER=mpiifort
```

## SCOTCH

SCOTCH requires adding the compile flag

```makefile
CFLAGS = other_flags -DINTSIZE64
```

in `Makefile.inc`.

## PASTIX

In Pastix-5, you need to change the section that specifies the integer size, to choose:

```makefile
VERSIONINT = _int64
CCTYPES = -DINTSSIZE64
```

## STRUMPACK

In STRUMPACK, you will need to change the line

```cpp
template<typename scalar_t,typename integer_t>
```

to

```cpp
template<typename scalar_t,typename integer_t=int64_t>
```

in the files:

- `StrumpackSparseSolver.hpp`
- `StrumpackSparseSolverBase.hpp`
- `StrumpackSparseSolverMPIDist.hpp`

For STRUMPACK, on Marconi, the CMake commands before compiling are:

```bash
export ScaLAPACKLIBS="-L$MKLROOT/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl"
```

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_C_COMPILER=mpiicc -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_CXX_FLAGS=-std=c++14 -DSTRUMPACK_USE_OPENMP=ON -DTPL_ENABLE_BPACK=OFF -DTPL_ENABLE_ZFP=OFF -DTPL_ENABLE_SLATE=OFF -DTPL_METIS_INCLUDE_DIRS=/marconi_work/FUA34_ELM-UK/spamela/STRUMPACK_tries/try_intel_64/metis-5.1.0/install/include -DTPL_METIS_LIBRARIES="-L/marconi_work/FUA34_ELM-UK/spamela/STRUMPACK_tries/try_intel_64/metis-5.1.0/install/lib -lmetis" -DTPL_SCALAPACK_LIBRARIES="$ScaLAPACKLIBS" ../
```

Of course, you will need to update the locations of METIS/PARMETIS accordingly.

## JOREK

In JOREK, simply use the directive

```makefile
USE_INTSIZE64 = 1
```

in `Makefile.inc`.
