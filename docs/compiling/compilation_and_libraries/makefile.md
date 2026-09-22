---
title: "Build System"
nav_order: 2
parent: "Compilation and Libraries"
grand_parent: "Compiling and Running"
layout: default
render_with_liquid: false
---

# The build system

We use `make` to compile the code.
The build system consists of three parts:

- `Makefile`
- `Makefile.inc`
- `defaults.mk`

and two supporting scripts

- `util/makedepend`
- `util/obj_deps`

If you wish to understand more about the build system it is recommended to read the source also.

## The main Makefile

This is the file read by your version of `make`.

It contains the default rules (`clean` etc.) and a list of directories to look for source files in.
For each of the found source files (matching `*.f90`) we generate a dependency file (`.d`) with `util/obj_deps` containing build and link-time dependencies for this source file.
A line will also be added if this source file provides a `PROGRAM`.

`make` will, on first reading the file, automatically try to build all of the `.d` files.
After this, the dependency graph is constructed and your programs can be compiled.
Programs will thus be autodiscovered based on the presence of the `PROGRAM` keyword in the file.

## defaults.mk

This file contains default settings, flags and rules for building Fortran sources.
Some flags are defined, based on which compiler you use, and the build templates are converted into rules.

## Makefile.inc

In this file you can set which libraries to use (`HDF5`, `HIPS`, `MUMPS`, `WSMP`, `PASTIX`, `MURGE`, `FFTW`, `LAPACK`, `BLAS`).
You can also set specific compile flags, select which model to use and provide link options.
Many examples are provided in `Make.inc/`

## util/makedepend

This is a simple Perl script, inspired by `sfmakedepend`, to read through a Fortran source file and look for

- `use` statements
- `#include`
- `include`
- `module`
- `program`
- `subroutine calls`
- `external functions`

After looking for all external procedures to call and modules to load, the program will try to localize the files in the passed list of directories, using `find`.
For this to work, the filename must match the procedure name.
Some known translations are performed (such as `pplot -> ppplib`). The current list of transformations is not exhaustive and should probably be amended in the future.

## util/obj_deps

This is an even simpler Perl script, which reads through a tree of `.d` dependency files to create a list of dependencies required for linking.
This is necessary because we need to keep different sets of build and link dependencies, and `make` cannot do this for us.
You will probably not need to alter this file.
