# Source discovery.
#
# The Make build listed source directories in DIRS and then used
# util/makedepend + util/obj_deps to work out, per executable, which objects to
# link -- by scanning for `call foo(...)` and pulling in a file named foo.f90.
# CMake needs none of that: Fortran module dependencies are tracked natively,
# and the link closure falls out of putting every non-program source into a
# static archive and letting the linker extract what it references.

# Directories searched for sources, non-recursively (as in the Makefile's DIRS).
set(JOREK_SOURCE_DIRS
  .
  benchmarks
  communication
  core
  datatypes
  diagnostics
  diagnostics/new_diag
  diagnostics/postproc
  elements
  grids
  grids/grid_utils
  matrix
  models
  models/${JOREK_MODEL}
  numerics
  particles
  particles/diagnostics
  particles/postprocessors
  particles/postprocessors/camera
  particles/postprocessors/filters
  particles/postprocessors/lens
  particles/postprocessors/lights
  particles/postprocessors/spectra
  particles/postprocessors/utils
  particles/initialisers
  particles/projection_functions
  particles/pushers
  # Not optional despite the name: the synthetic-light postprocessors in
  # particles/postprocessors/lights `use mod_particle_common_test_tools`, which
  # lives here. The Make build always had this directory in DIRS, so the
  # coupling went unnoticed.
  particles/tests
  plots
  refinement
  solvers
  tools
  tools/fruit
  vacuum)

# The IMAS coupling only compiles against the ITER Access Layer.
if(JOREK_USE_IMAS)
  list(APPEND JOREK_SOURCE_DIRS communication/IMAS)
endif()

# The FRUIT test suites are not part of a normal build; the Make build only
# compiled them on demand too. The suites in particles/tests `use` modules from
# particles/benchmarks, so enabling one requires the other.
if(JOREK_BUILD_TESTS)
  list(APPEND JOREK_SOURCE_DIRS
    communication/tests
    core/tests
    diagnostics/tests
    elements/tests
    grids/tests
    particles/postprocessors/tests
    tools/tests
    particles/benchmarks/projection
    particles/benchmarks/pusher
    particles/benchmarks/pusher_cartesian)
endif()

if(JOREK_BUILD_EXAMPLES)
  list(APPEND JOREK_SOURCE_DIRS
    particles/examples
    particles/postprocessors/examples)
endif()

list(REMOVE_DUPLICATES JOREK_SOURCE_DIRS)

set(_all_sources "")
foreach(_dir IN LISTS JOREK_SOURCE_DIRS)
  file(GLOB _found CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/${_dir}/*.f90"
    "${CMAKE_SOURCE_DIR}/${_dir}/*.f"
    "${CMAKE_SOURCE_DIR}/${_dir}/*.c"
    "${CMAKE_SOURCE_DIR}/${_dir}/*.cpp")
  list(APPEND _all_sources ${_found})
endforeach()

# The Catalyst adaptor is guarded by USE_CATALYST but its .cpp does not compile
# at all without the Catalyst headers, so drop it rather than rely on #ifdef.
if(NOT JOREK_USE_CATALYST)
  list(FILTER _all_sources EXCLUDE REGEX "/(mod_)?catalyst_adaptor\\.(f90|cpp)$")
endif()

# communication/export_grid.f90 and communication/import_equil.f90 write and
# read type_element with unformatted derived-type I/O. For the stellarator
# models type_element gains an allocatable `chi` component, which makes that
# I/O invalid Fortran, so neither file compiles. Nothing in the tree calls
# either of them -- they are dead code, and the Make build never compiled them
# because it only built what the target program transitively referenced.
if(JOREK_STELLARATOR)
  list(APPEND JOREK_EXCLUDE_SOURCES
    "/communication/export_grid\\.f90$"
    "/communication/import_equil\\.f90$")
endif()

foreach(_pattern IN LISTS JOREK_EXCLUDE_SOURCES)
  list(FILTER _all_sources EXCLUDE REGEX "${_pattern}")
endforeach()

# algexpr2fort is a build-time code generator, not part of the solver; it is
# given its own target below.
list(FILTER _all_sources EXCLUDE REGEX "/algexpr2fort\\.f90$")

# --- Test-directory helpers required by production code ---------------------
#
# particles/tests holds two different kinds of file: the FRUIT test suites, and
# shared helper modules that non-test code genuinely depends on (the synthetic
# light postprocessors `use mod_particle_common_test_tools`). The Make build
# kept every tests/ directory in DIRS and compiled lazily, so it never had to
# tell them apart.
#
# Selecting by filename would be wrong: mod_projection_helpers_test_tools has an
# unconditional `include 'dmumps_struc.h'` and only compiles against MUMPS, yet
# nothing outside the benchmarks needs it. So pull in a tests/ file only when
# something outside tests/ actually uses a module it defines, transitively.
if(NOT JOREK_BUILD_TESTS)
  set(_pool ${_all_sources})
  set(_core ${_all_sources})
  list(FILTER _pool INCLUDE REGEX "/tests/")
  list(FILTER _core EXCLUDE REGEX "/tests/")

  # Which tests/ file provides which module?
  foreach(_src IN LISTS _pool)
    file(STRINGS "${_src}" _lines
      REGEX "^[ \t]*[Mm][Oo][Dd][Uu][Ll][Ee][ \t]+[A-Za-z]")
    foreach(_line IN LISTS _lines)
      if(_line MATCHES "^[ \t]*[Mm][Oo][Dd][Uu][Ll][Ee][ \t]+([A-Za-z][A-Za-z0-9_]*)")
        string(TOLOWER "${CMAKE_MATCH_1}" _m)
        if(NOT _m STREQUAL "procedure")
          set(_provider_${_m} "${_src}")
        endif()
      endif()
    endforeach()
  endforeach()

  # Seed the worklist with every module used outside tests/, then close over the
  # helpers that are pulled in.
  set(_wanted "")
  foreach(_src IN LISTS _core)
    if(NOT _src MATCHES "\\.f90$")
      continue()
    endif()
    file(STRINGS "${_src}" _lines REGEX "^[ \t]*[Uu][Ss][Ee][ \t]+[A-Za-z]")
    foreach(_line IN LISTS _lines)
      if(_line MATCHES "^[ \t]*[Uu][Ss][Ee][ \t]+([A-Za-z][A-Za-z0-9_]*)")
        string(TOLOWER "${CMAKE_MATCH_1}" _m)
        list(APPEND _wanted "${_m}")
      endif()
    endforeach()
  endforeach()
  list(REMOVE_DUPLICATES _wanted)

  set(_kept_helpers "")
  set(_i 0)
  list(LENGTH _wanted _n_wanted)
  while(_i LESS _n_wanted)
    list(GET _wanted ${_i} _m)
    math(EXPR _i "${_i} + 1")
    if(NOT DEFINED _provider_${_m})
      continue()
    endif()
    set(_src "${_provider_${_m}}")
    if(_src IN_LIST _kept_helpers)
      continue()
    endif()
    list(APPEND _kept_helpers "${_src}")

    file(STRINGS "${_src}" _lines REGEX "^[ \t]*[Uu][Ss][Ee][ \t]+[A-Za-z]")
    foreach(_line IN LISTS _lines)
      if(_line MATCHES "^[ \t]*[Uu][Ss][Ee][ \t]+([A-Za-z][A-Za-z0-9_]*)")
        string(TOLOWER "${CMAKE_MATCH_1}" _dep)
        if(NOT _dep IN_LIST _wanted)
          list(APPEND _wanted "${_dep}")
          math(EXPR _n_wanted "${_n_wanted} + 1")
        endif()
      endif()
    endforeach()
  endwhile()

  set(_all_sources ${_core} ${_kept_helpers})

  list(LENGTH _kept_helpers _n_helpers)
  if(_n_helpers)
    message(STATUS
      "JOREK: pulling in ${_n_helpers} helper module(s) from tests/ that non-test code depends on")
  endif()
endif()

# --- Split programs from library sources ------------------------------------
# A Fortran file with a `program` statement becomes its own executable;
# everything else goes into the core library.
set(JOREK_PROGRAM_SOURCES "")
set(JOREK_LIB_SOURCES "")
foreach(_src IN LISTS _all_sources)
  set(_is_program FALSE)
  if(_src MATCHES "\\.f90$")
    file(STRINGS "${_src}" _program_stmt
      REGEX "^[ \t]*[Pp][Rr][Oo][Gg][Rr][Aa][Mm][ \t]+[A-Za-z]"
      LIMIT_COUNT 1)
    if(_program_stmt)
      set(_is_program TRUE)
    endif()
  endif()

  if(_is_program)
    list(APPEND JOREK_PROGRAM_SOURCES "${_src}")
  else()
    list(APPEND JOREK_LIB_SOURCES "${_src}")
  endif()
endforeach()

if(NOT JOREK_PROGRAM_SOURCES MATCHES "jorek2_main\\.f90")
  message(FATAL_ERROR "jorek2_main.f90 was not found among the discovered programs")
endif()

# Module and object files are keyed by the source's basename, so two sources
# with the same name in different directories collide. The Makefile had a
# `duplicates` target to report this; resolve it at configure time instead.
#
# The Make build let the first directory in DIRS win, silently shadowing the
# others. JOREK_SOURCE_DIRS is ordered to give the same winner, so this keeps
# the two builds producing the same binary -- but say so loudly, because a
# shadowed file is almost always a mistake in the source tree.
function(_jorek_drop_shadowed list_var)
  set(_kept "")
  set(_seen "")
  foreach(_src IN LISTS ${list_var})
    get_filename_component(_name "${_src}" NAME_WE)
    if(_name IN_LIST _seen)
      list(APPEND JOREK_SHADOWED_SOURCES "${_src}")
    else()
      list(APPEND _seen "${_name}")
      list(APPEND _kept "${_src}")
    endif()
  endforeach()
  set(${list_var} "${_kept}" PARENT_SCOPE)
  set(JOREK_SHADOWED_SOURCES "${JOREK_SHADOWED_SOURCES}" PARENT_SCOPE)
endfunction()

set(JOREK_SHADOWED_SOURCES "")
set(_combined ${JOREK_LIB_SOURCES} ${JOREK_PROGRAM_SOURCES})
_jorek_drop_shadowed(_combined)
foreach(_src IN LISTS JOREK_SHADOWED_SOURCES)
  list(REMOVE_ITEM JOREK_LIB_SOURCES "${_src}")
  list(REMOVE_ITEM JOREK_PROGRAM_SOURCES "${_src}")
  get_filename_component(_name "${_src}" NAME_WE)
  file(RELATIVE_PATH _rel "${CMAKE_SOURCE_DIR}" "${_src}")
  message(WARNING
    "Duplicate source basename '${_name}': ignoring ${_rel}, which is shadowed by "
    "another file of the same name earlier in JOREK_SOURCE_DIRS. This matches what "
    "the Make build linked, but the duplicate should be removed from the repository.")
endforeach()

# --- Sources for the semianalytical code generator ---------------------------
#
# Models 180/183 build their element matrices from symbolic expressions that the
# algexpr2fort program turns into Fortran at build time. algexpr2fort is itself
# built from JOREK sources, so it cannot link the core library: that library
# contains mod_elt_matrix, which #includes the very headers algexpr2fort has not
# produced yet.
#
# Give it instead its own small library holding exactly what it needs -- the
# transitive closure of the modules it `use`s, which is about a dozen files and
# reaches neither mod_elt_matrix nor the generated headers. That keeps the
# dependency graph acyclic and leaves every executable linking a single archive.
# The names a Fortran source depends on: the modules it `use`s, plus the bare
# subroutines it `call`s. JOREK resolves the latter by filename -- a `call foo`
# is satisfied by whichever file is named foo.f90 -- which is what
# util/makedepend did and what the linker relies on.
function(_jorek_source_deps out_var file)
  set(_names "")
  if(NOT file MATCHES "\\.(f90|f)$")
    set(${out_var} "" PARENT_SCOPE)   # C/C++ pull in nothing by these rules
    return()
  endif()

  file(STRINGS "${file}" _lines)
  foreach(_line IN LISTS _lines)
    string(REGEX REPLACE "!.*$" "" _line "${_line}")   # drop comments

    if(_line MATCHES "^[ \t]*[Uu][Ss][Ee][ \t]+([A-Za-z][A-Za-z0-9_]*)")
      string(TOLOWER "${CMAKE_MATCH_1}" _m)
      list(APPEND _names "${_m}")
    endif()

    string(REGEX MATCHALL "[Cc][Aa][Ll][Ll][ \t]+[A-Za-z][A-Za-z0-9_]*" _hits "${_line}")
    foreach(_hit IN LISTS _hits)
      string(REGEX REPLACE "^[Cc][Aa][Ll][Ll][ \t]+" "" _hit "${_hit}")
      string(TOLOWER "${_hit}" _hit)
      list(APPEND _names "${_hit}")
    endforeach()
  endforeach()

  if(_names)
    list(REMOVE_DUPLICATES _names)
  endif()
  set(${out_var} "${_names}" PARENT_SCOPE)
endfunction()

function(_jorek_forward_closure out_var seed sources)
  # Index the pool twice: module name -> providing file, and basename -> file.
  foreach(_src IN LISTS sources)
    get_filename_component(_base "${_src}" NAME_WE)
    string(TOLOWER "${_base}" _base)
    if(NOT DEFINED _basefile_${_base})
      set(_basefile_${_base} "${_src}")
    endif()

    if(NOT _src MATCHES "\\.f90$")
      continue()
    endif()
    file(STRINGS "${_src}" _lines
      REGEX "^[ \t]*[Mm][Oo][Dd][Uu][Ll][Ee][ \t]+[A-Za-z]")
    foreach(_line IN LISTS _lines)
      if(_line MATCHES "^[ \t]*[Mm][Oo][Dd][Uu][Ll][Ee][ \t]+([A-Za-z][A-Za-z0-9_]*)")
        string(TOLOWER "${CMAKE_MATCH_1}" _m)
        # `module procedure` in an interface block is not a module definition.
        if(NOT _m STREQUAL "procedure")
          set(_provider_${_m} "${_src}")
        endif()
      endif()
    endforeach()
  endforeach()

  _jorek_source_deps(_pending "${seed}")
  set(_closure "")
  set(_i 0)
  list(LENGTH _pending _n)
  while(_i LESS _n)
    list(GET _pending ${_i} _name)
    math(EXPR _i "${_i} + 1")

    # A name may resolve as a module, as a filename, or both.
    set(_hits "")
    if(DEFINED _provider_${_name})
      list(APPEND _hits "${_provider_${_name}}")
    endif()
    if(DEFINED _basefile_${_name})
      list(APPEND _hits "${_basefile_${_name}}")
    endif()

    foreach(_src IN LISTS _hits)
      if(_src IN_LIST _closure)
        continue()
      endif()
      list(APPEND _closure "${_src}")

      _jorek_source_deps(_deps "${_src}")
      foreach(_dep IN LISTS _deps)
        if(NOT _dep IN_LIST _pending)
          list(APPEND _pending "${_dep}")
          math(EXPR _n "${_n} + 1")
        endif()
      endforeach()
    endforeach()
  endwhile()

  set(${out_var} "${_closure}" PARENT_SCOPE)
endfunction()

set(JOREK_GENERATED_INCLUDE_DIR "${CMAKE_BINARY_DIR}/generated/models/${JOREK_MODEL}")
file(MAKE_DIRECTORY "${JOREK_GENERATED_INCLUDE_DIR}")

set(JOREK_CODEGEN_SOURCES "")
if(JOREK_SEMIANALYTICAL)
  _jorek_forward_closure(JOREK_CODEGEN_SOURCES
    "${CMAKE_SOURCE_DIR}/models/algexpr2fort.f90" "${JOREK_LIB_SOURCES}")

  # If this ever trips, the equations module has grown a dependency on the code
  # it is supposed to generate, and the generator can no longer be built.
  foreach(_src IN LISTS JOREK_CODEGEN_SOURCES)
    if(_src MATCHES "/mod_elt_matrix(_fft)?\\.f90$")
      message(FATAL_ERROR
        "algexpr2fort transitively depends on ${_src}, which #includes the headers "
        "algexpr2fort generates. This dependency cycle must be broken in the sources.")
    endif()
  endforeach()

  list(LENGTH JOREK_CODEGEN_SOURCES _n_codegen)
  message(STATUS "JOREK: algexpr2fort code generator needs ${_n_codegen} source(s)")
endif()

# The .f sources take the relaxed diagnostics from JorekCompilerFlags. Source
# properties are scoped to the directory rather than to a target, so this one
# pass covers every target -- they are all defined in the top-level scope.
if(JOREK_Fortran_FIXED_FORM_OPTIONS)
  set(_fixed_form ${JOREK_LIB_SOURCES} ${JOREK_PROGRAM_SOURCES} ${JOREK_CODEGEN_SOURCES})
  list(FILTER _fixed_form INCLUDE REGEX "\\.f$")
  if(_fixed_form)
    list(REMOVE_DUPLICATES _fixed_form)
    set_source_files_properties(${_fixed_form} PROPERTIES
      COMPILE_OPTIONS "${JOREK_Fortran_FIXED_FORM_OPTIONS}")
  endif()
endif()

list(LENGTH JOREK_LIB_SOURCES _n_lib)
list(LENGTH JOREK_PROGRAM_SOURCES _n_prog)
message(STATUS "JOREK sources: ${_n_lib} library, ${_n_prog} program(s)")
