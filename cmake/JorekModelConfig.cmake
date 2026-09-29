# Derives preprocessor symbols from the physics model, reproducing what
# `util/config.sh -p <key>` did for defaults.mk. The model settings files are
# the single source of truth: we read them rather than duplicating their values.

set(JOREK_MODEL_DIR "${CMAKE_SOURCE_DIR}/models/${JOREK_MODEL}")
set(JOREK_MODEL_SETTINGS "${JOREK_MODEL_DIR}/mod_model_settings.f90")
set(JOREK_GLOBAL_SETTINGS "${CMAKE_SOURCE_DIR}/models/mod_settings.f90")

if(NOT IS_DIRECTORY "${JOREK_MODEL_DIR}")
  file(GLOB _available RELATIVE "${CMAKE_SOURCE_DIR}/models" "${CMAKE_SOURCE_DIR}/models/model*")
  message(FATAL_ERROR "Unknown model '${JOREK_MODEL}'. Available: ${_available}")
endif()
foreach(_f "${JOREK_MODEL_SETTINGS}" "${JOREK_GLOBAL_SETTINGS}")
  if(NOT EXISTS "${_f}")
    message(FATAL_ERROR "Missing model settings file: ${_f}")
  endif()
endforeach()

# Re-run CMake when the settings files change: they decide the compile flags.
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${JOREK_MODEL_SETTINGS}" "${JOREK_GLOBAL_SETTINGS}")

if(NOT JOREK_MODEL MATCHES "^model([0-9]+)$")
  message(FATAL_ERROR "JOREK_MODEL must look like 'model<number>', got '${JOREK_MODEL}'")
endif()
set(JOREK_MODEL_NUMBER "${CMAKE_MATCH_1}")

# Read a `<type>, parameter :: <key> = <value>` declaration out of a settings
# file.
#
# Fortran is case-insensitive and the model settings files are not consistent:
# model600 writes `with_vpar` while model303 writes `with_Vpar`. util/config.sh
# greps with -i, so match case-insensitively here too -- a case-sensitive match
# silently drops the flag and builds the wrong physics.
function(_jorek_read_setting out_var file key)
  string(TOLOWER "${key}" _key)
  file(STRINGS "${file}" _lines)
  set(${out_var} "" PARENT_SCOPE)
  foreach(_line IN LISTS _lines)
    string(TOLOWER "${_line}" _line)
    # Strip trailing comments so a value mentioned in the documentation column
    # cannot be picked up instead of the real one.
    string(REGEX REPLACE "!.*$" "" _line "${_line}")
    if(_line MATCHES "::[ \t]*${_key}[ \t]*=[ \t]*([^ \t]+)")
      set(${out_var} "${CMAKE_MATCH_1}" PARENT_SCOPE)
      return()
    endif()
  endforeach()
endfunction()

# Map a `logical, parameter :: <key> = .true.` setting onto a -D<define>.
function(_jorek_logical_define key define)
  _jorek_read_setting(_value "${JOREK_MODEL_SETTINGS}" "${key}")
  string(TOLOWER "${_value}" _value)
  if(_value STREQUAL ".true.")
    list(APPEND JOREK_DEFINES "${define}")
    set(JOREK_DEFINES "${JOREK_DEFINES}" PARENT_SCOPE)
  elseif(NOT _value STREQUAL ".false." AND NOT _value STREQUAL "")
    message(FATAL_ERROR "Could not parse '${key}' in ${JOREK_MODEL_SETTINGS} (got '${_value}')")
  endif()
endfunction()

set(JOREK_DEFINES "")

list(APPEND JOREK_DEFINES "JOREK_MODEL=${JOREK_MODEL_NUMBER}" "USE_MPI")

_jorek_logical_define(with_Vpar       WITH_Vpar)
_jorek_logical_define(with_TiTe       WITH_TiTe)
_jorek_logical_define(with_neutrals   WITH_Neutrals)
_jorek_logical_define(with_impurities WITH_Impurities)
_jorek_logical_define(with_refluid    WITH_REFluid)

# A non-zero n_mod_ext means the model is the base of a model family.
_jorek_read_setting(JOREK_N_MOD_EXT "${JOREK_MODEL_SETTINGS}" "n_mod_ext")
if(JOREK_N_MOD_EXT AND NOT JOREK_N_MOD_EXT STREQUAL "0")
  list(APPEND JOREK_DEFINES "MODEL_FAMILY")
endif()

# Polynomial orders above 3 need more Gauss points than the compiled-in default.
_jorek_read_setting(JOREK_N_ORDER "${JOREK_GLOBAL_SETTINGS}" "n_order")
if(NOT JOREK_N_ORDER)
  message(FATAL_ERROR "Could not read n_order from ${JOREK_GLOBAL_SETTINGS}")
endif()
if(NOT JOREK_N_ORDER STREQUAL "3")
  list(APPEND JOREK_DEFINES "GAUSS_ORDER=8")
endif()

# Full-MHD models take a different equation set.
set(_fullmhd_models 710 711 712 750)
if(JOREK_MODEL_NUMBER IN_LIST _fullmhd_models)
  list(APPEND JOREK_DEFINES "fullmhd")
endif()

# The stellarator models build their element matrices from symbolic expressions
# that algexpr2fort turns into Fortran at build time. STELLARATOR_MODEL also
# changes type_element (it gains an allocatable `chi`), so the two flags are
# tracked separately even though the same models set both today.
set(_semianalytical_models 180 183)
if(JOREK_MODEL_NUMBER IN_LIST _semianalytical_models)
  set(JOREK_SEMIANALYTICAL TRUE)
  set(JOREK_STELLARATOR TRUE)
  list(APPEND JOREK_DEFINES "SEMIANALYTICAL" "STELLARATOR_MODEL")
else()
  set(JOREK_SEMIANALYTICAL FALSE)
  set(JOREK_STELLARATOR FALSE)
endif()

# --- Feature switches (JorekOptions) ---------------------------------------
if(JOREK_USE_DOMM)
  list(APPEND JOREK_DEFINES "USE_DOMM")
endif()
if(JOREK_USE_EXT_FIELD)
  list(APPEND JOREK_DEFINES "USE_EXT_FIELD")
endif()
if(JOREK_USE_BLOCK)
  list(APPEND JOREK_DEFINES "USE_BLOCK")
endif()
if(JOREK_USE_DIRECT_CONSTRUCTION)
  list(APPEND JOREK_DEFINES "DIRECT_CONSTRUCTION")
endif()
if(JOREK_USE_COMPLEX_PRECOND)
  list(APPEND JOREK_DEFINES "USE_COMPLEX_PRECOND")
endif()
if(JOREK_USE_INTSIZE64)
  list(APPEND JOREK_DEFINES "INTSIZE64")
endif()
if(JOREK_USE_TASKLOOP)
  list(APPEND JOREK_DEFINES "USE_TASKLOOP")
endif()
if(JOREK_USE_QUADTREE)
  list(APPEND JOREK_DEFINES "USE_QUADTREE")
endif()
if(JOREK_USE_NO_TREE)
  list(APPEND JOREK_DEFINES "USE_NO_TREE")
endif()
if(JOREK_USE_R3_INFO)
  list(APPEND JOREK_DEFINES "USE_R3_INFO" "USE_R3_INFO_MPI")
endif()

# GMRES is the default; BiCGSTAB replaces it.
if(JOREK_USE_BICGSTAB)
  list(APPEND JOREK_DEFINES "USE_BICGSTAB")
else()
  list(APPEND JOREK_DEFINES "USE_GMRES")
endif()

if(JOREK_USE_STD_BESSELK)
  list(APPEND JOREK_DEFINES "USE_STD_BESSELK")
endif()
if(JOREK_USE_BOOST)
  list(APPEND JOREK_DEFINES "USE_BOOST")
endif()
