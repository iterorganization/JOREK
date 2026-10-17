# Per-compiler flags, replacing the COMPILER_FAMILY blocks of defaults.mk.
#
# Two things are non-negotiable for JOREK and are set here rather than left to
# the site config:
#   * Fortran sources must be run through the C preprocessor (-cpp / -fpp).
#     The sources are named .f90, so no compiler does this automatically.
#   * Default reals must be promoted to 8 bytes (-fdefault-real-8 / -r8).
#     The code assumes it; a build without it is silently wrong.

if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL CMAKE_C_COMPILER_ID)
  message(FATAL_ERROR
    "Fortran (${CMAKE_Fortran_COMPILER_ID}) and C (${CMAKE_C_COMPILER_ID}) compilers "
    "must come from the same vendor")
endif()

set(_common "")   # applied to every Fortran source
set(_debug "")    # applied to Fortran sources in Debug builds only

# Relaxes the implicit-typing diagnostics for the vendored fixed-form .f sources,
# which declare nothing explicitly; the Make build had the same carve-out in
# F77FLAGS. Applied per source file by JorekSources.
set(JOREK_Fortran_FIXED_FORM_OPTIONS "")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

  list(APPEND _common
    -cpp -ffree-line-length-none
    -fdefault-real-8 -fdefault-double-8
    # JOREK passes arrays where scalars are declared in several MPI calls;
    # gfortran >= 10 rejects that by default.
    -fallow-argument-mismatch
    -Wall -Wextra -Wno-unused-variable -Wno-tabs
    -Wcharacter-truncation -Wsurprising)

  list(APPEND _debug
    -fcheck=all
    -ffpe-trap=invalid,zero,overflow
    -finit-real=snan -finit-integer=12345678
    -fimplicit-none)

  # gfortran makes -fimplicit-none an error, so without this the fixed-form
  # sources do not compile at all in a Debug build.
  set(JOREK_Fortran_FIXED_FORM_OPTIONS -fno-implicit-none)

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")

  # -warn and -check take their keyword list as a separate argument -- ifx warns
  # (#10158) and drops the option for a `-warn=...` spelling -- so keep the two
  # tokens together in a `SHELL:` group, which is also what stops CMake from
  # de-duplicating them apart. -init is the opposite case and wants `-init=`.
  list(APPEND _common
    -fpp -r8 -align
    "SHELL:-warn all,nointerfaces,nounused,noexternal"
    # #5462 "global name too long" fires all over the tree: ifx truncates mangled
    # symbols, keeping the distinguishing tail, and the only fix would be
    # renaming modules. No -warn keyword covers it short of -warn nogeneral.
    -diag-disable=5462)

  list(APPEND _debug
    -traceback
    "SHELL:-check all,noarg_temp_created"
    -ftrapuv -fpe0 -init=snan
    -implicitnone)

  # Kills #6717 from both -warn all and -implicitnone. Two plain list elements,
  # not a `SHELL:` group: source-file COMPILE_OPTIONS pass `SHELL:` through to
  # the compiler verbatim, and are not de-duplicated anyway.
  set(JOREK_Fortran_FIXED_FORM_OPTIONS -warn nodeclarations)

  # The semianalytical models build very large automatic arrays in the element
  # matrix routines and overflow the stack without this.
  if(JOREK_SEMIANALYTICAL)
    list(APPEND _common -heap-arrays)
  endif()

else()
  message(WARNING
    "Untested Fortran compiler '${CMAKE_Fortran_COMPILER_ID}'. You must supply "
    "preprocessing and 8-byte-default-real flags yourself via CMAKE_Fortran_FLAGS.")
endif()

add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${_common}>")
add_compile_options(
  "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:${_debug}>")

# A plain-text record of the flags for version.h. Generator expressions cannot
# be evaluated at configure time, so keep the flags as a string here rather than
# trying to recover them from the COMPILE_OPTIONS directory property.
set(_record ${_common})
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  list(APPEND _record ${_debug})
endif()
string(REPLACE ";" " " JOREK_Fortran_FLAGS_RECORD "${_record}")
# `SHELL:` is a CMake marker, not something the compiler ever sees.
string(REPLACE "SHELL:" "" JOREK_Fortran_FLAGS_RECORD "${JOREK_Fortran_FLAGS_RECORD}")

find_package(OpenMP REQUIRED COMPONENTS Fortran C CXX)
