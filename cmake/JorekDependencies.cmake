# External libraries. Everything discovered here is appended to
# JOREK_LIBRARIES and JOREK_DEFINES, which the core library exposes to all
# executables.

set(JOREK_LIBRARIES "")

# --- MPI --------------------------------------------------------------------
find_package(MPI REQUIRED COMPONENTS Fortran C)
list(APPEND JOREK_LIBRARIES MPI::MPI_Fortran MPI::MPI_C)
list(APPEND JOREK_LIBRARIES OpenMP::OpenMP_Fortran)

# JOREK compiles differently against MPI-2 and MPI-3. Fortran defines
# MPI_VERSION as an integer parameter, which cannot be used by the
# preprocessor, so the value has to be pushed in from outside. The Make build
# compiled and ran a probe program; FindMPI already reports it.
if(MPI_Fortran_VERSION MATCHES "^([0-9]+)")
  list(APPEND JOREK_DEFINES "MPI_VERSION=${CMAKE_MATCH_1}")
else()
  message(WARNING "Could not determine the MPI version; assuming MPI_VERSION=0")
  list(APPEND JOREK_DEFINES "MPI_VERSION=0")
endif()

# --- Dense linear algebra ---------------------------------------------------
# MKL provides BLAS, LAPACK, ScaLAPACK, BLACS and an FFTW3 interface in one
# package; otherwise the three are found separately.
if(JOREK_USE_MKL)
  if(NOT DEFINED ENV{MKLROOT})
    message(FATAL_ERROR "JOREK_USE_MKL=ON but MKLROOT is not set; load the MKL module")
  endif()
  set(JOREK_MKL_ROOT "$ENV{MKLROOT}")

  if(JOREK_USE_INTSIZE64)
    set(_mkl_int ilp64)
  else()
    set(_mkl_int lp64)
  endif()

  # Link the MKL interface layer as a group: the libraries are mutually
  # recursive and a plain sequential link line does not resolve.
  set(_mkl_libs
    -L${JOREK_MKL_ROOT}/lib/intel64
    -Wl,--start-group
    -lmkl_scalapack_${_mkl_int}
    -lmkl_intel_${_mkl_int}
    -lmkl_intel_thread
    -lmkl_core
    -lmkl_blacs_intelmpi_${_mkl_int}
    -Wl,--end-group
    -liomp5 -lpthread -lm -ldl)
  list(APPEND JOREK_LIBRARIES ${_mkl_libs})
  set(JOREK_BLAS_LAPACK_SUMMARY "Intel MKL (${JOREK_MKL_ROOT})")
else()
  if(JOREK_BLAS_LAPACK_LIBRARIES)
    list(APPEND JOREK_LIBRARIES ${JOREK_BLAS_LAPACK_LIBRARIES})
    set(JOREK_BLAS_LAPACK_SUMMARY "${JOREK_BLAS_LAPACK_LIBRARIES} (user-specified)")
  else()
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
    list(APPEND JOREK_LIBRARIES LAPACK::LAPACK BLAS::BLAS)
    set(JOREK_BLAS_LAPACK_SUMMARY "${LAPACK_LIBRARIES}")
  endif()

  if(JOREK_SCALAPACK_LIBRARIES)
    list(APPEND JOREK_LIBRARIES ${JOREK_SCALAPACK_LIBRARIES})
  endif()
endif()

# --- HDF5 -------------------------------------------------------------------
if(JOREK_USE_HDF5)
  # The site modules export HDF5_HOME; FindHDF5 looks at HDF5_ROOT.
  if(NOT HDF5_ROOT AND DEFINED ENV{HDF5_HOME})
    set(HDF5_ROOT "$ENV{HDF5_HOME}")
  endif()
  set(HDF5_PREFER_PARALLEL TRUE)
  find_package(HDF5 REQUIRED COMPONENTS Fortran Fortran_HL)
  list(APPEND JOREK_LIBRARIES ${HDF5_Fortran_HL_LIBRARIES} ${HDF5_Fortran_LIBRARIES})
  list(APPEND JOREK_DEFINES "USE_HDF5")
  set(JOREK_HDF5_SUMMARY "${HDF5_VERSION} (${HDF5_Fortran_INCLUDE_DIRS})")

  if(NOT HDF5_IS_PARALLEL)
    message(WARNING "The HDF5 found is serial; JOREK expects a parallel (MPI) HDF5")
  endif()
endif()

# --- FFTW -------------------------------------------------------------------
if(JOREK_USE_FFTW)
  list(APPEND JOREK_DEFINES "USE_FFTW")
  if(JOREK_USE_MKL)
    # MKL ships the FFTW3 interface, including the fftw3.f03 module file that
    # jorek2_main.f90 includes. No extra library to link.
    set(JOREK_FFTW_INCLUDE_DIR "${JOREK_MKL_ROOT}/include/fftw")
    if(NOT EXISTS "${JOREK_FFTW_INCLUDE_DIR}/fftw3.f03")
      message(FATAL_ERROR
        "MKL FFTW3 interface not found at ${JOREK_FFTW_INCLUDE_DIR}")
    endif()
    set(JOREK_FFTW_SUMMARY "MKL FFTW3 interface")
  else()
    find_path(JOREK_FFTW_INCLUDE_DIR fftw3.f03
      HINTS $ENV{FFTW_INCLUDE} $ENV{FFTW_HOME}/include
      DOC "Directory containing fftw3.f03")
    find_library(JOREK_FFTW_LIBRARY fftw3
      HINTS $ENV{FFTW_LIB} $ENV{FFTW_HOME}/lib
      DOC "The FFTW3 library")
    if(NOT JOREK_FFTW_INCLUDE_DIR OR NOT JOREK_FFTW_LIBRARY)
      message(FATAL_ERROR "JOREK_USE_FFTW=ON but FFTW was not found; load the fftw module")
    endif()
    list(APPEND JOREK_LIBRARIES ${JOREK_FFTW_LIBRARY})
    set(JOREK_FFTW_SUMMARY "${JOREK_FFTW_LIBRARY}")
  endif()
endif()

# --- Sparse direct solvers --------------------------------------------------
# Small helper: locate a library under a required install prefix and fail with
# an actionable message rather than a link error 20 minutes later.
function(_jorek_require_library out_var name)
  cmake_parse_arguments(ARG "" "" "HINTS" ${ARGN})
  find_library(${out_var} ${name}
    HINTS ${ARG_HINTS}
    PATH_SUFFIXES lib lib64
    NO_DEFAULT_PATH)
  find_library(${out_var} ${name} HINTS ${ARG_HINTS} PATH_SUFFIXES lib lib64)
  if(NOT ${out_var})
    message(FATAL_ERROR "Could not find lib${name}; searched: ${ARG_HINTS}")
  endif()
  set(${out_var} "${${out_var}}" PARENT_SCOPE)
endfunction()

if(JOREK_USE_STRUMPACK)
  if(NOT STRUMPACK_ROOT)
    message(FATAL_ERROR "JOREK_USE_STRUMPACK=ON but STRUMPACK_ROOT is not set")
  endif()
  _jorek_require_library(JOREK_STRUMPACK_LIBRARY strumpack HINTS "${STRUMPACK_ROOT}")
  find_path(JOREK_STRUMPACK_INCLUDE_DIR StrumpackSparseSolver.h
    HINTS "${STRUMPACK_ROOT}" PATH_SUFFIXES include)

  list(APPEND JOREK_LIBRARIES ${JOREK_STRUMPACK_LIBRARY})
  list(APPEND JOREK_DEFINES "USE_STRUMPACK" "NEWSPK")
  set(JOREK_STRUMPACK_SUMMARY "${JOREK_STRUMPACK_LIBRARY}")

  # STRUMPACK's sparse reordering pulls in the METIS family.
  foreach(_m METIS PARMETIS GKLIB)
    if(${_m}_ROOT)
      string(TOLOWER "${_m}" _lib)
      if(_m STREQUAL "GKLIB")
        set(_lib GKlib)
      endif()
      _jorek_require_library(JOREK_${_m}_LIBRARY ${_lib} HINTS "${${_m}_ROOT}")
      list(APPEND JOREK_LIBRARIES ${JOREK_${_m}_LIBRARY})
      find_path(JOREK_${_m}_INCLUDE_DIR NAMES metis.h parmetis.h GKlib.h
        HINTS "${${_m}_ROOT}" PATH_SUFFIXES include)
      if(JOREK_${_m}_INCLUDE_DIR)
        list(APPEND JOREK_SOLVER_INCLUDE_DIRS "${JOREK_${_m}_INCLUDE_DIR}")
      endif()
    endif()
  endforeach()
endif()

if(JOREK_USE_MUMPS)
  if(NOT MUMPS_ROOT)
    message(FATAL_ERROR "JOREK_USE_MUMPS=ON but MUMPS_ROOT is not set")
  endif()
  _jorek_require_library(JOREK_DMUMPS_LIBRARY dmumps HINTS "${MUMPS_ROOT}")
  _jorek_require_library(JOREK_MUMPS_COMMON_LIBRARY mumps_common HINTS "${MUMPS_ROOT}")
  find_library(JOREK_PORD_LIBRARY pord
    HINTS "${MUMPS_ROOT}" "${MUMPS_ROOT}/PORD" PATH_SUFFIXES lib lib64)
  find_path(JOREK_MUMPS_INCLUDE_DIR dmumps_struc.h
    HINTS "${MUMPS_ROOT}" PATH_SUFFIXES include)

  list(APPEND JOREK_LIBRARIES
    ${JOREK_DMUMPS_LIBRARY} ${JOREK_MUMPS_COMMON_LIBRARY} ${JOREK_PORD_LIBRARY})
  list(APPEND JOREK_SOLVER_INCLUDE_DIRS "${JOREK_MUMPS_INCLUDE_DIR}")
  list(APPEND JOREK_DEFINES "USE_MUMPS")
  set(JOREK_MUMPS_SUMMARY "${JOREK_DMUMPS_LIBRARY}")
endif()

if(JOREK_USE_PASTIX OR JOREK_USE_PASTIX6)
  if(NOT PASTIX_ROOT)
    message(FATAL_ERROR "PaStiX requested but PASTIX_ROOT is not set")
  endif()
  _jorek_require_library(JOREK_PASTIX_LIBRARY pastix HINTS "${PASTIX_ROOT}")
  find_path(JOREK_PASTIX_INCLUDE_DIR NAMES pastix.h pastix_fortran.h
    HINTS "${PASTIX_ROOT}" PATH_SUFFIXES include)
  list(APPEND JOREK_LIBRARIES ${JOREK_PASTIX_LIBRARY})
  list(APPEND JOREK_SOLVER_INCLUDE_DIRS "${JOREK_PASTIX_INCLUDE_DIR}")

  if(JOREK_USE_PASTIX)
    list(APPEND JOREK_DEFINES "USE_PASTIX" "MEMORY_USAGE")
  endif()
  if(JOREK_USE_PASTIX6)
    list(APPEND JOREK_DEFINES "USE_PASTIX6")
  endif()
  if(JOREK_USE_PASTIX_MURGE)
    list(APPEND JOREK_DEFINES "USE_MURGE")
  endif()

  if(SCOTCH_ROOT)
    _jorek_require_library(JOREK_SCOTCH_LIBRARY scotch HINTS "${SCOTCH_ROOT}")
    _jorek_require_library(JOREK_SCOTCHERR_LIBRARY scotcherr HINTS "${SCOTCH_ROOT}")
    list(APPEND JOREK_LIBRARIES ${JOREK_SCOTCH_LIBRARY} ${JOREK_SCOTCHERR_LIBRARY})
  endif()
else()
  # mod_pastix.f90 keeps its interface declarations even when PaStiX is absent.
  # Renaming the entry point avoids unresolved symbols at link time; this is
  # the same trick defaults.mk used.
  list(APPEND JOREK_DEFINES "pastix_fortran=fake_pastix_fortran")
endif()

if(JOREK_USE_WSMP)
  _jorek_require_library(JOREK_WSMP_LIBRARY pwsmp64 HINTS "$ENV{WSMP_HOME}")
  list(APPEND JOREK_LIBRARIES ${JOREK_WSMP_LIBRARY})
  list(APPEND JOREK_DEFINES "USE_WSMP")
endif()

if(JOREK_USE_HIPS)
  _jorek_require_library(JOREK_HIPS_LIBRARY hips HINTS "$ENV{HIPS_HOME}")
  list(APPEND JOREK_LIBRARIES ${JOREK_HIPS_LIBRARY})
  list(APPEND JOREK_DEFINES "USE_HIPS")
endif()

# --- Optional integrations --------------------------------------------------
if(JOREK_USE_BOOST)
  find_package(Boost REQUIRED)
  list(APPEND JOREK_LIBRARIES Boost::boost)
endif()

if(JOREK_USE_IMAS)
  find_path(JOREK_IMAS_INCLUDE_DIR ids_schemas.mod HINTS $ENV{IMAS_PREFIX}/include)
  find_library(JOREK_IMAS_LIBRARY imas HINTS $ENV{IMAS_PREFIX}/lib)
  if(NOT JOREK_IMAS_INCLUDE_DIR OR NOT JOREK_IMAS_LIBRARY)
    message(FATAL_ERROR "JOREK_USE_IMAS=ON but the IMAS Access Layer was not found")
  endif()
  list(APPEND JOREK_LIBRARIES ${JOREK_IMAS_LIBRARY})
  list(APPEND JOREK_SOLVER_INCLUDE_DIRS "${JOREK_IMAS_INCLUDE_DIR}")
  list(APPEND JOREK_DEFINES "USE_IMAS")
endif()

if(JOREK_USE_CATALYST)
  find_package(catalyst REQUIRED)
  list(APPEND JOREK_LIBRARIES catalyst::catalyst)
  list(APPEND JOREK_DEFINES "USE_CATALYST")
endif()

# JOREK links C++ translation units (the r-tree, the STRUMPACK/PaStiX shims)
# into a Fortran executable, so the C++ runtime must be pulled in explicitly.
list(APPEND JOREK_LIBRARIES stdc++)

# --- Collected include directories -----------------------------------------
set(JOREK_INCLUDE_DIRS
  ${JOREK_SOLVER_INCLUDE_DIRS}
  ${JOREK_STRUMPACK_INCLUDE_DIR}
  ${JOREK_FFTW_INCLUDE_DIR}
  ${HDF5_Fortran_INCLUDE_DIRS})
list(REMOVE_ITEM JOREK_INCLUDE_DIRS "")
list(REMOVE_DUPLICATES JOREK_INCLUDE_DIRS)
