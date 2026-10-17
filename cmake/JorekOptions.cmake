# Cache options replacing the switches that used to live in Makefile.inc.
# Every JOREK_USE_* here maps onto exactly one -DUSE_* preprocessor symbol,
# with the same meaning it had in defaults.mk.

set(JOREK_MODEL "model600" CACHE STRING
  "Physics model to build (directory name under models/, e.g. model600)")

# --- Physics / numerics switches -------------------------------------------
option(JOREK_USE_DOMM "Use Dommaschk potentials (no FE correction of n.B on boundary)" ON)
option(JOREK_USE_EXT_FIELD "Use external magnetic field from a gvec2jorek file (requires JOREK_USE_DOMM=OFF)" OFF)
option(JOREK_USE_BLOCK "Enable blocked matrix assembly" ON)
option(JOREK_USE_DIRECT_CONSTRUCTION "Construct the matrix directly" OFF)
option(JOREK_USE_BICGSTAB "Use BiCGSTAB instead of GMRES as the iterative solver" OFF)
option(JOREK_USE_COMPLEX_PRECOND "Use a complex-valued preconditioner" OFF)
option(JOREK_USE_INTSIZE64 "Link against the 64-bit-integer builds of the solver libraries" OFF)
option(JOREK_USE_TASKLOOP "Use OpenMP taskloop constructs" OFF)
option(JOREK_USE_R3_INFO "Enable the r3_info timing instrumentation" OFF)

# Element localisation strategy: the r-tree is the default, the other two are
# mutually exclusive alternatives selected by preprocessor symbol.
option(JOREK_USE_QUADTREE "Use a quadtree instead of an r-tree for element localisation" OFF)
option(JOREK_USE_NO_TREE "Use brute-force element localisation (no spatial tree)" OFF)

# --- Sparse solvers ---------------------------------------------------------
option(JOREK_USE_STRUMPACK "Build with the STRUMPACK sparse direct solver" ON)
option(JOREK_USE_MUMPS "Build with the MUMPS sparse direct solver" OFF)
option(JOREK_USE_PASTIX "Build with the PaStiX 5 sparse direct solver" OFF)
option(JOREK_USE_PASTIX_MURGE "Use the PaStiX MURGE interface" OFF)
option(JOREK_USE_PASTIX6 "Build with the PaStiX 6 sparse direct solver" OFF)
option(JOREK_USE_WSMP "Build with the WSMP sparse direct solver" OFF)
option(JOREK_USE_HIPS "Build with the HIPS sparse direct solver" OFF)

# --- I/O and optional integrations -----------------------------------------
option(JOREK_USE_HDF5 "Build with HDF5 input/output (strongly recommended)" ON)
option(JOREK_USE_FFTW "Build with FFTW (or the MKL FFTW3 interface)" ON)
option(JOREK_USE_MKL "Use Intel MKL for BLAS/LAPACK/ScaLAPACK/FFTW" OFF)
option(JOREK_USE_IMAS "Build the IMAS/IDS coupling (communication/IMAS)" OFF)
option(JOREK_USE_CATALYST "Build the ParaView Catalyst in-situ adaptor" OFF)

# --- Bessel function backend (mutually exclusive) ---------------------------
option(JOREK_USE_STD_BESSELK "Use the C++ standard library for modified Bessel functions" OFF)
option(JOREK_USE_BOOST "Use Boost for modified Bessel functions" OFF)

# --- What to build ----------------------------------------------------------
option(JOREK_BUILD_TESTS "Compile the FRUIT unit-test modules into the core library" OFF)
option(JOREK_BUILD_EXAMPLES "Compile the particle examples and benchmarks" OFF)

# --- External library locations ---------------------------------------------
# Presets fill these in; they are plain paths so a site can point them anywhere.
set(JOREK_LIB_DIR "" CACHE PATH "Root of the pre-built JOREK library collection")
set(STRUMPACK_ROOT "" CACHE PATH "STRUMPACK install prefix")
set(METIS_ROOT "" CACHE PATH "METIS install prefix")
set(PARMETIS_ROOT "" CACHE PATH "ParMETIS install prefix")
set(GKLIB_ROOT "" CACHE PATH "GKlib install prefix")
set(MUMPS_ROOT "" CACHE PATH "MUMPS install prefix")
set(PASTIX_ROOT "" CACHE PATH "PaStiX install prefix")
set(SCOTCH_ROOT "" CACHE PATH "Scotch install prefix")
set(JOREK_DIERCKX_LIBRARIES "" CACHE STRING "DIERCKX libraries (only needed by eqdsk2jorek)")

# Escape hatch for sources that do not compile in a given configuration.
# CMake compiles every source in the tree into the core library, whereas the
# Make build only compiled what the requested program transitively referenced --
# so CMake also compiles files that are unreachable, and therefore unmaintained.
set(JOREK_EXCLUDE_SOURCES "" CACHE STRING
  "Regexes (matched against absolute source paths) to exclude from the build")
set(JOREK_BLAS_LAPACK_LIBRARIES "" CACHE STRING "Override the detected BLAS/LAPACK libraries")
set(JOREK_SCALAPACK_LIBRARIES "" CACHE STRING "Override the detected ScaLAPACK/BLACS libraries")

# --- Consistency checks -----------------------------------------------------
if(JOREK_USE_QUADTREE AND JOREK_USE_NO_TREE)
  message(FATAL_ERROR "JOREK_USE_QUADTREE and JOREK_USE_NO_TREE are mutually exclusive")
endif()

if(JOREK_USE_STD_BESSELK AND JOREK_USE_BOOST)
  message(FATAL_ERROR "JOREK_USE_STD_BESSELK and JOREK_USE_BOOST are mutually exclusive")
endif()

# defaults.mk only honoured USE_EXT_FIELD when USE_DOMM was off; make that
# an error rather than silently dropping the flag.
if(JOREK_USE_EXT_FIELD AND JOREK_USE_DOMM)
  message(FATAL_ERROR
    "JOREK_USE_EXT_FIELD requires JOREK_USE_DOMM=OFF (the external field replaces "
    "the Dommaschk potentials)")
endif()

if(NOT JOREK_USE_HDF5)
  message(WARNING "JOREK_USE_HDF5=OFF: restart and diagnostic I/O will be limited")
endif()

if(NOT (JOREK_USE_STRUMPACK OR JOREK_USE_MUMPS OR JOREK_USE_PASTIX
        OR JOREK_USE_PASTIX6 OR JOREK_USE_WSMP OR JOREK_USE_HIPS))
  message(FATAL_ERROR "No sparse direct solver selected; enable at least one of "
    "JOREK_USE_STRUMPACK / JOREK_USE_MUMPS / JOREK_USE_PASTIX / JOREK_USE_PASTIX6")
endif()
