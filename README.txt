*******************************************************
*                    JOREK2                           *
*                                                     *
* copyright : Guido Huysmans (Association Euratom/CEA *
*******************************************************

For a documentation about the code, see jorek.eu/wiki


=======================================================
 Building with CMake (recommended)
=======================================================

Each supported machine is a CMake preset in CMakePresets.json. To list them:

  cmake --list-presets

Pick the preset for your machine (e.g. viper-intel, pitagora-intel; add
-debug for a debug build). Each preset's "description" field states the exact
`module load` line to run first.

  # 1. load the modules listed in the preset's description
  # 2. configure once
  cmake --preset <preset>

  # 3. build (build tree is build/<preset>/)
  cmake --build build/<preset> -j 16                        # default: main solver
  cmake --build build/<preset> --target jorek2_diagno -j 16 # a specific program
  cmake --build build/<preset> --target help                # list all targets

Every .f90 file with a `program` statement becomes a target named after the
file; the built executables land in build/<preset>/. Set the model and options
per preset in CMakePresets.json (JOREK_MODEL, JOREK_USE_*), not in Makefile.inc.


=======================================================
 Building with Makefile.inc (DEPRECATED)
=======================================================

NOTE: the Makefile.inc / `make` build below is deprecated. Prefer the CMake
build above. The instructions are kept for reference only.

to build jorek2:

- edit the Makefile.inc file to match your environment
  (some examples are available in the directory configs)
  
  choose a model : model199 (reduced MHD, no v_parallel)
                   model300 (reduced MHD, with V_parallel)

- edit the model dependent mod_parameters.f90 file to define
  n_tor    : the number of toroidal harmnonics (sin+cos) 
  n_plane  : the number of toroidal planes (for FFT)
  n_period : the periodicity of the torus
  
- make all (in top directory)

to execute jorek2:

mpirun -np n_mpi jorek2_model199 < namelist/model199/intear

this is a simple testcase (example output is available in test_output)

the number of mpi's should be a multiple of the number of harmonics ((n_tor-1)/2+1)


Topic: Timing

We use the "r3_info" package from R. Richter, SGI, Oct 2006 (trunk/timing) to compute elapsed time in JOREK.
Ask G. Latu (guillaume.latu@cea.fr) for more details.
To use it, we need
 - Fortran preprocessing options: -cpp -DUSE_R3_INFO_MPI -DUSE_R3_INFO (see trunk/configs/config.hpcff)
 - define a MPICC compiler and preprocessing options: -DUSE_R3_INFO_MPI -DUSE_R3_INFO (see trunk/configs/config.hpcff)
 - need #include "r3_info.h" statement (warning: C syntax with #) in Fortran file (see trunk/jorek2_main.f90 for example)
