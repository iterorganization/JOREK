# version.h records the git revision, the compile flags and who built the code;
# constants.f90, the restart I/O and jorek2help all #include it and print it
# into the log so a run can be traced back to a build.
#
# It must be regenerated on every build, not just at configure time, because
# the git description changes as you commit. GenerateVersionHeader.cmake only
# rewrites the file when its contents actually change, so this does not force a
# rebuild of the world on every `make`.

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/generated")
set(JOREK_VERSION_HEADER "${CMAKE_BINARY_DIR}/generated/version.h")

# The header embeds the exact flags used: the site flags for this build type
# plus the JOREK ones recorded by JorekCompilerFlags.
string(TOUPPER "${CMAKE_BUILD_TYPE}" _cfg)
set(_flags "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_cfg}} ${JOREK_Fortran_FLAGS_RECORD}")

string(REPLACE ";" " " _defines_str "${JOREK_DEFINES}")
string(REPLACE ";" " " _includes_str "${JOREK_INCLUDE_DIRS}")
string(REPLACE ";" " " _libs_str "${JOREK_LIBRARIES}")

add_custom_target(jorek_version_header
  BYPRODUCTS "${JOREK_VERSION_HEADER}"
  COMMAND ${CMAKE_COMMAND}
    -DOUTPUT_FILE=${JOREK_VERSION_HEADER}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
    -DBINARY_DIR=${CMAKE_BINARY_DIR}
    -DCOMPILE_COMMAND=${CMAKE_Fortran_COMPILER}
    "-DCOMPILE_FLAGS=${_flags}"
    "-DCOMPILE_DEFINES=${_defines_str}"
    "-DCOMPILE_INCLUDES=${_includes_str}"
    "-DCOMPILE_LIBS=${_libs_str}"
    -P "${CMAKE_SOURCE_DIR}/cmake/GenerateVersionHeader.cmake"
  COMMENT "Checking JOREK version header"
  VERBATIM)

# Run it once now so the header exists before the first compile.
execute_process(COMMAND ${CMAKE_COMMAND}
  -DOUTPUT_FILE=${JOREK_VERSION_HEADER}
  -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
  -DBINARY_DIR=${CMAKE_BINARY_DIR}
  -DCOMPILE_COMMAND=${CMAKE_Fortran_COMPILER}
  "-DCOMPILE_FLAGS=${_flags}"
  "-DCOMPILE_DEFINES=${_defines_str}"
  "-DCOMPILE_INCLUDES=${_includes_str}"
  "-DCOMPILE_LIBS=${_libs_str}"
  -P "${CMAKE_SOURCE_DIR}/cmake/GenerateVersionHeader.cmake")
