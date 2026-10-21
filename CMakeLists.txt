cmake_minimum_required(VERSION 3.20)

# The version header embeds the compile flags, so languages must be enabled
# before we can generate it. Fortran is the primary language; C and C++ are
# needed for tools/*.c, the r-tree, the sparse-solver shims and Catalyst.
project(JOREK
  LANGUAGES Fortran C CXX
  DESCRIPTION "Non-linear extended MHD code for toroidal X-point geometries")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR
    "In-source builds are not supported. Configure into a separate directory, "
    "e.g. cmake --preset pitagora-intel")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(JorekOptions)         # cache options, replaces Makefile.inc switches
include(JorekModelConfig)     # derive -D flags from the model settings files
include(JorekCompilerFlags)   # per-compiler flag sets, replaces defaults.mk
include(JorekDependencies)    # MPI, HDF5, FFTW, BLAS/LAPACK, sparse solvers
include(JorekVersion)         # git-derived version.h
include(JorekSources)         # source discovery, program/library split

# --- Usage requirements shared by every JOREK target ------------------------
add_library(jorek_common INTERFACE)
target_include_directories(jorek_common INTERFACE
  "${CMAKE_BINARY_DIR}/generated"        # version.h
  "${JOREK_GENERATED_INCLUDE_DIR}"       # rhs/amat/aux_automatic.h
  "${CMAKE_CURRENT_SOURCE_DIR}/tools"    # r3_info.h
  "${CMAKE_CURRENT_SOURCE_DIR}/models"
  ${JOREK_INCLUDE_DIRS})
target_compile_definitions(jorek_common INTERFACE ${JOREK_DEFINES})
target_link_libraries(jorek_common INTERFACE ${JOREK_LIBRARIES})

# --- The code generator (semianalytical models only) ------------------------
#
# Built before, and independently of, the core library -- see JorekSources for
# why it cannot simply link against it. Its objects are compiled a second time
# into jorek_core; that is a dozen files, and it keeps the graph acyclic.
if(JOREK_SEMIANALYTICAL)
  add_library(jorek_codegen STATIC ${JOREK_CODEGEN_SOURCES})
  add_dependencies(jorek_codegen jorek_version_header)
  target_link_libraries(jorek_codegen PRIVATE jorek_common)
  # Its own module directory: these sources are also compiled into jorek_core,
  # and the two must not fight over the same .mod files.
  set_target_properties(jorek_codegen PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules_codegen")

  # algexpr2fort evaluates the model's symbolic equations and writes the
  # RHS/AMAT/AUX Fortran that mod_elt_matrix{,_fft}.f90 #include. It hard-codes
  # its output path as "models/model<NNN>/<name>.h" relative to the working
  # directory, so run it from the build tree to keep the source tree clean.
  add_executable(algexpr2fort "${CMAKE_CURRENT_SOURCE_DIR}/models/algexpr2fort.f90")
  target_link_libraries(algexpr2fort PRIVATE jorek_codegen jorek_common)
  set_target_properties(algexpr2fort PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules_codegen"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    LINKER_LANGUAGE Fortran)

  set(JOREK_GENERATED_HEADERS
    "${JOREK_GENERATED_INCLUDE_DIR}/rhs_automatic.h"
    "${JOREK_GENERATED_INCLUDE_DIR}/amat_automatic.h"
    "${JOREK_GENERATED_INCLUDE_DIR}/aux_automatic.h")

  add_custom_command(
    OUTPUT ${JOREK_GENERATED_HEADERS}
    COMMAND "$<TARGET_FILE:algexpr2fort>"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/generated"
    DEPENDS algexpr2fort
    COMMENT "Generating element-matrix code for ${JOREK_MODEL}"
    VERBATIM)
  add_custom_target(jorek_generated_code DEPENDS ${JOREK_GENERATED_HEADERS})
endif()

# --- The core library -------------------------------------------------------
#
# The Make build discovered each executable's objects by scanning for `call`
# statements and linking whatever file was named after the callee. Here every
# non-program source instead goes into one static archive, and the linker
# extracts what each executable actually references. Keeping it a single archive
# matters: JOREK is full of bare external subroutines, so any split into two
# libraries would leave them referencing each other across the link line.
add_library(jorek_core STATIC ${JOREK_LIB_SOURCES})
add_dependencies(jorek_core jorek_version_header)
target_link_libraries(jorek_core PUBLIC jorek_common)
target_include_directories(jorek_core PUBLIC "${CMAKE_Fortran_MODULE_DIRECTORY}")

if(JOREK_SEMIANALYTICAL)
  add_dependencies(jorek_core jorek_generated_code)
endif()

# --- Executables ------------------------------------------------------------
#
# Only the main solver is built by default, mirroring `make`. Every other
# program is a named target; `make most` builds the regression-tested subset.
foreach(_prog_src IN LISTS JOREK_PROGRAM_SOURCES)
  get_filename_component(_prog "${_prog_src}" NAME_WE)

  if(_prog STREQUAL "jorek2_main")
    add_executable(jorek2_main "${_prog_src}")
    set_target_properties(jorek2_main PROPERTIES
      OUTPUT_NAME "jorek_model${JOREK_MODEL_NUMBER}")
  else()
    add_executable(${_prog} EXCLUDE_FROM_ALL "${_prog_src}")
  endif()

  target_link_libraries(${_prog} PRIVATE jorek_core)
  set_target_properties(${_prog} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    # The core library mixes Fortran, C and C++ objects, and CMake would
    # otherwise pick the C++ driver to link -- which does not supply the
    # Fortran runtime's main(). Every JOREK program is a Fortran program.
    LINKER_LANGUAGE Fortran)
endforeach()

if(TARGET eqdsk2jorek AND JOREK_DIERCKX_LIBRARIES)
  target_link_libraries(eqdsk2jorek PRIVATE ${JOREK_DIERCKX_LIBRARIES})
endif()

# The set of diagnostics that must keep compiling for regression testing to
# pass. Kept in sync with the `most:` target of the legacy Makefile.
set(JOREK_MOST_TARGETS
  jorek2_connection2
  jorek2_connection_stan
  jorek2_connection_flux_aligned
  jorek2_diagno
  jorek2_diagno_spi
  jorek2_fieldlines_vtk
  jorek2_four
  CARIDDI_wall_curr
  jorek2_poincare
  jorek2_powers
  jorek2_target2vtk
  jorek2vtk_3d
  jorek2vtk
  jorek2vtk_GaussVortTerms
  test_gvec2jorek_import
  jorek_to_helena
  new_diag_demo
  jorek2_postproc
  rst_bin2hdf5
  rst_hdf52bin
  jorek2_main)
if(JOREK_USE_IMAS)
  list(APPEND JOREK_MOST_TARGETS jorek2_IDS)
endif()

set(_most_existing "")
foreach(_t IN LISTS JOREK_MOST_TARGETS)
  if(TARGET ${_t})
    list(APPEND _most_existing ${_t})
  else()
    message(WARNING "'most' target references unknown program '${_t}'")
  endif()
endforeach()
add_custom_target(most DEPENDS ${_most_existing})

include(JorekSummary)
