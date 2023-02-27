# Using ParaView Catalyst in-situ visualization with JOREK

## Why use in-situ visualization?

Whilst the processing power of supercomputers has grown rapidly over recent 
years, the corresponding growth in I/O bandwidth and data storage has not kept
pace. This has led to an increased need to run diagnostics and perform analysis 
pipelines for HPC simulations at runtime. [ParaView 
Catalyst](https://www.paraview.org/hpc-insitu/) is one such in-situ framework
that allows performing visualization and data analysis during a simulation.

Although JOREK has relatively modest resource requirements and disk capacity or
I/O bandwidth limitations may not be a concern for most users, there are still
reasons to consider trying out Catalyst. For example, once the pipelines are set
up, the workflow can be more efficient as one does not need to post-process the
data after running a simulation. One can also monitor simulations better.

In this page we describe how to use the in-situ visualization capabilities of
ParaView Catalyst with a JOREK model.

## Prerequisites

> **Note** This adaptor uses the new Catalyst2 framework which is documented 
> [here](https://catalyst-in-situ.readthedocs.io/en/latest/introduction.html) 
> and is available in ParaView v5.9 or greater. This is different to the 
> previous Catalyst framework which was the only one available in older versions 
> of ParaView (sometimes referred to as Catalyst1 or Legacy Catalyst).

In order to build JOREK with Catalyst enabled, you need to compile and link with
a development version of a _Catalyst implementation_ (i.e. a library that
implements the _Catalyst API_). If there is a version of ParaView on the target
system, it might have been built with this (Note that the pre-built versions of
ParaView available on the ParaView website do not come with the necessary
development header files). You can check by navigating to the ParaView
installation directory and checking if there is a `catalyst.hpp` header file
under the following path:
```
/path/to/paraview/include/catalyst-2.0/catalyst.hpp
```
If not, you will need to build a suitable version of a Catalyst implementation.
It is easiest to build the stub Catalyst implementation (rather than building
ParaView from scratch which is significantly bigger and more complicated).
Documentation on how to that can be found 
[here](https://catalyst-in-situ.readthedocs.io/en/latest/build_and_install.html).
Note that, there is also a 
[Spack package available](https://packages.spack.io/package.html?name=libcatalyst)
if you prefer.

Note that it is not necessary to build JOREK with the same Catalyst 
implementation you intend to use. We will switch to using the ParaView Catalyst
implementation at runtime as the stub implementation does nothing.

You will need ParaView 5.9 or greater to run JOREK with ParaView's Catalyst
implementation. This can be downloaded from the [ParaView
website](https://www.paraview.org/download/).

## Building JOREK with Catalyst

> **Note** Currently only model199 has Catalyst support enabled.

First make sure you are on a branch with Catalyst support. This will have the
file `jorek/diagnostics/mod_catalyst_adaptor.f90`. Add the following lines to 
your JOREK `Makefile.inc` file:
```makefile
USE_CATALYST = 1
CATALYST_HOME = /path/to/catalyst/implementation
CATALYSTLIB = -L$(CATALYST_HOME)/lib64 -lcatalyst
CATALYSTINCLUDE = $(CATALYST_HOME)/include/catalyst-2.0
```

Note that there are some C++ files that will need to be compiled so make sure
you have set an appropriate C++ compiler in the CXX variable in your 
`Makefile.inc` file.

Build JOREK as normal with a command such as
```
make -j 8
```

## Running JOREK with Catalyst

You will need to provide a Catalyst Python script to describe a Catalyst
pipeline. Some example scripts are provided in the `jorek/diagnostics/catalyst`
directory. See below for some documentation on how to generate a script for a 
new pipeline. To pass a script to the code, set the namelist parameter 
`catalyst_scripts` to a string that is the path to the script. Multiple scripts
can be provided by separating paths with a `:` (make sure your paths do not 
include `:`). For example
```fortran
 catalyst_scripts = "/path/to/catalyst_script1.py:/path/to/catalyst_script2.py"
```

In order to load the ParaView Catalyst implementation at runtime and not use the
stub implementation (which does nothing), you will need to set the 
following environment variables
```bash
export CATALYST_IMPLEMENTATION_NAME=paraview
export CATALYST_IMPLEMENTATION_PATHS=/path/to/paraview/lib/catalyst/
```
Note that the installation of ParaView has the necessary Catalyst support if you 
can see the file `libcatalyst-paraview.so` under `paraview/lib/catalyst/` 
(this is what Catalyst will look for if you set the environment variables as
above). If this doesn't exist, you may need to ask for ParaView to be rebuilt 
with Catalyst2 (`libcatalyst` in the Spack package) support.

Assuming everything is working, you should see the following output when you run 
the simulation (Catalyst is initialized towards the end of the initialization 
phase of the code)
```
...
**********************************
* Initializing Catalyst pipeline *
**********************************
Catalyst info:
catalyst: 
  version: "2.0"
  abi_version: "2"
  implementation: "paraview"
  use_mpi: 1
...
```

If there are issues loading ParaView's Catalyst implementation, you can set the
`CATALYST_DEBUG` environment variable to a non-empty value to log the search and
loading procedure.

You can control the verbosity of ParaView's output by setting the enviroment
variable `PARAVIEW_LOG_CATALYST_VERBOSITY` appropriately. For example
```bash
export PARAVIEW_LOG_CATALYST_VERBOSITY=OFF # print nothing
export PARAVIEW_LOG_CATALYST_VERBOSITY=ERROR # only print errors
export PARAVIEW_LOG_CATALYST_VERBOSITY=WARNING # print errors and warning
export PARAVIEW_LOG_CATALYST_VERBOSITY=INFO # print errors, warning and info
export PARAVIEW_LOG_CATALYST_VERBOSITY=TRACE # max verbosity, print everything
```

## Generating a ParaView Catalyst script

The easiest way to generate a suitable ParaView Catalyst Python script using the
GUI is to load a JOREK VTK file generated from data from a previous simulation.
Once you have loaded the file into ParaView,
**make sure you rename the source to 'grid'** (by right clicking on it and
selecting 'Rename'). This is necessary so that Catalyst can replace the data
source in the script with the simulation data

Next apply filters and adjust the camera view as desired. Note that not all
variables may be available during the simulation.

In order to generate a valid Catalyst script, it is then necessary to add an 
'Extractor' which will, for example, save an image or VTK file. This can be done
using the menu options, Extractors → Image → PNG.

Finally, to generate the Catalyst script, use the menu options File → 
Save Catalyst State. By default, the output directory for Extracts will be
'dataset'.

## Implementation details

Currently, only model199 has Catalyst support added. However, the only part that
is model specific are the changes to the `initialise_parameters` subroutine so
it is straightforward to add support to another model.

The adaptor currently creates and passes a 2D grid that is similar to that
created by the `jorek2vtk` program. In particular, JOREK elements are subdivided
into `catalyst_nsub * catalyst_nsub = 5 * 5` elements and the components of the
variables with respect to the toroidal harmonics are interpolated onto this new
grid and then projected onto the planar slice (`i_plane = 1` i.e. `phi = 0`).
Currently, only the variables in `mod_model_settings` are interpolated onto this
new grid so there are `n_var` variables passed.