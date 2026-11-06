# Writes version.h. Run as a script (cmake -P) both at configure time and on
# every build, so the recorded git revision follows the working tree.
#
# The values are Fortran character literals, so any single quote inside them
# must be doubled or the #include will not compile.

function(_git out_var)
  execute_process(
    COMMAND git ${ARGN}
    WORKING_DIRECTORY "${SOURCE_DIR}"
    OUTPUT_VARIABLE _out
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _rc)
  if(NOT _rc EQUAL 0)
    set(_out "unknown")
  endif()
  set(${out_var} "${_out}" PARENT_SCOPE)
endfunction()

function(_quote out_var value)
  string(REPLACE "'" "''" value "${value}")
  string(REPLACE "\n" " " value "${value}")
  set(${out_var} "${value}" PARENT_SCOPE)
endfunction()

_git(RCS_VERSION describe --always --dirty --abbrev)
_git(LATEST_TAG_REV rev-list --tags --max-count=1)
if(LATEST_TAG_REV STREQUAL "unknown" OR LATEST_TAG_REV STREQUAL "")
  set(JOREK_VERSION "unknown")
else()
  _git(JOREK_VERSION describe --tags "${LATEST_TAG_REV}")
endif()
_git(RCS_LABEL log -1 "--format=%s (%D)")
_git(RCS_TIME log -1 "--format=%ad")

cmake_host_system_information(RESULT COMPILE_MACHINE QUERY HOSTNAME)
if(DEFINED ENV{USER})
  set(COMPILE_USER "$ENV{USER}")
else()
  set(COMPILE_USER "unknown")
endif()
# Recording the loaded environment modules makes a build reproducible on an HPC
# system, where the module set determines which libraries were linked.
if(DEFINED ENV{LOADEDMODULES})
  set(COMPILE_MODULES "$ENV{LOADEDMODULES}")
else()
  set(COMPILE_MODULES "")
endif()

foreach(_v RCS_VERSION JOREK_VERSION RCS_LABEL RCS_TIME COMPILE_COMMAND
           COMPILE_FLAGS COMPILE_INCLUDES COMPILE_DEFINES COMPILE_LIBS
           COMPILE_USER COMPILE_MACHINE COMPILE_MODULES)
  _quote(${_v} "${${_v}}")
endforeach()

set(_content
"#define RCS_VERSION '${RCS_VERSION}'
#define JOREK_VERSION '${JOREK_VERSION}'
#define RCS_LABEL '${RCS_LABEL}'
#define RCS_TIME '${RCS_TIME}'
#define compile_command '${COMPILE_COMMAND}'
#define compile_flags '${COMPILE_FLAGS}'
#define compile_includes '${COMPILE_INCLUDES}'
#define compile_defines '${COMPILE_DEFINES}'
#define compile_libs '${COMPILE_LIBS}'
#define compile_dir '${BINARY_DIR}'
#define compile_user '${COMPILE_USER}'
#define compile_machine '${COMPILE_MACHINE}'
#define compile_modules '${COMPILE_MODULES}'
")

# Only touch the file when it changes, so rebuilding does not invalidate every
# object that includes it.
set(_existing "")
if(EXISTS "${OUTPUT_FILE}")
  file(READ "${OUTPUT_FILE}" _existing)
endif()
if(NOT _existing STREQUAL _content)
  file(WRITE "${OUTPUT_FILE}" "${_content}")
endif()
