---
title: "rst_bin2hdf5"
nav_order: 19
render_with_liquid: false
parent: "Running JOREK"
---

## rst_bin2hdf5
- See also [hdf5-tools](hdf5-tools.md)

This tool is used to convert jorek00XXX.rst files to hdf5 format.
The first argument is the file to convert. If no argument is present, it will use jorek_restart.rst.
Example usage:

  rst_bin2hdf5 jorek00001.rst

Other options:

  ```
  rst_bin2hdf5 -v --verbose
  rst_bin2hdf5 -h --help
  ```

##### Compiling
To compile this, add `USE_HDF5 = 1` to your `Makefile.inc`, and add the following lines showing the linker where the libraries are. (`$H5DIR` is an existing environment variable)

  ```
  #HDF5
  HDF5LIB = -L$(H5DIR)/lib/ -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lhdf5_hl
  HDF5INCLUDE = $(H5DIR)/include
  ```

## rst_hdf52bin
This tool performs the opposite process as the one outlined above.
