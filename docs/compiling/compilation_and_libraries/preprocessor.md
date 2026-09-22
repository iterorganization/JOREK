---
title: "Preprocessor Flags"
nav_order: 4
parent: "Compilation and Libraries"
grand_parent: "Compiling and Running"
layout: default
render_with_liquid: false
---

# Preprocessor Directives

JOREK uses the **fpp preprocessor** ([Information from Intel](https://software.intel.com/de-de/node/510269)). `fpp` processes the source files before they are actually compiled. This allows for conditional compilation. **Preprocessor directives** are put into the code (starting with `#`) and allow testing for **preprocessor symbols** (like `USE_BLOCK`, `JOREK_MODEL`).

## Two examples for preprocessor directives

- Encapsulate model-specific parts:

```fortran
#if JOREK_MODEL == 500
use mgi_module
#endif
```

- Implement different solutions for the same problem, e.g. for the matrix analysis by the solver:

```fortran
#ifdef USE_BLOCK
  pastix_iparm(IPARM_DOF_NBR) = block_size
#else
  pastix_iparm(IPARM_DOF_NBR) = 1
#endif
```

## Setting preprocessor symbols

The preprocessor symbols like `USE_BLOCK` are normally defined by command line options like

```text
-DUSE_BLOCK
```

for the Intel Fortran compiler. In JOREK, we have defined switches for this purpose, such that it is sufficient to set

```make
USE_BLOCK = 1
```

in `Makefile.inc` for the most important switches (see table below for an overview).

## List of important preprocessor symbols

The following table summarizes the most important preprocessor symbols and how they are set.

| Preprocessor Symbol | Set in Makefile.inc | Explanation |
|----------|--------|--------|
| USE_MUMPS | USE_MUMPS = 1 | Use the MUMPS Solver |
| USE_PASTIX | USE_PASTIX = 1 | Use the PaStiX Solver (**recommended**) |
| USE_MURGE | USE_PASTIX_MURGE = 1 | Use the MURGE interface of PaStiX |
| USE_HIPS | USE_HIPS = 1 | Use the HIPS Solver |
| USE_WSMP | USE_WSMP = 1 | Use the WSMP Solver |
| USE_FFTW | USE_FFTW = 1 | Use the FFTW library for fast Fourier transform (**recommended**) |
| USE_BLOCK | USE_BLOCK = 1 | Make use of the block structure of the matrix for analyzing it (**recommended**) |
| USE_HDF5 | USE_HDF5 = 1 | Use the HDF5 data format |
| MEMTRACE | FFLAGS := $(FFLAGS) -DMEMTRACE | Output memory usage information in `trace*` files |
| COMPARE_ELEMENT_MATRIX | FFLAGS := $(FFLAGS) -DCOMPARE_ELEMENT_MATRIX | Compare `element_matrix` and `element_matrix_fft` directly for debugging |
| JOREK_MODEL | set automatically | Allows preprocessor directives like `#if JOREK_MODEL == 500` |
| FULL_MHD | set automatically | For model 710 |
| FUNNELED | FFLAGS := $(FFLAGS) -DFUNNELED | Flag for PaStiX to not call MPI from OpenMP parallel code (**required on many machines**) |
| WORLDWAR2 | FFLAGS := $(FFLAGS) -DWORLDWAR2 | Flag for PaXtiX (**required for `pastix_release_3945` and newer**, see also compiling) |

For the libraries, you need to set additional variables in `Makefile.inc`, for the FFTW library for instance `LIBFFTW` and `INC_FFTW`. Check the Makefile to find out which variables you need to specify for the libraries.
