cmake_minimum_required(VERSION 3.18)

project(impylib LANGUAGES C Fortran)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
include(CMakePrintHelpers)
include(F2Py)
include(FindPythonLibsNew)

if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  message(
    WARNING
      "You should call this through setup.py so that all variables are set; e.g. python setup.py develop"
  )
endif()

# cannot use find_package(Python), we need PythonLibsNew from pybind11
find_package(PythonLibsNew REQUIRED)

# in Numpy-1.22+, this becomes easier: import numpy.f2py; print(numpy.f2py.get_include())
if (NOT F2PY_INCLUDE_DIR)
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import numpy.f2py; from pathlib import Path; print(Path(numpy.f2py.__file__).parent)"
    OUTPUT_VARIABLE _f2py_directory
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  set(F2PY_INCLUDE_DIR "${_f2py_directory}/src" CACHE STRING "F2PY source directory location" FORCE)
endif()

if (NOT NUMPY_INCLUDE_DIRS)
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE _numpy_directory
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  set(NUMPY_INCLUDE_DIRS "${_numpy_directory}" CACHE STRING "NumPy source directory location" FORCE)
endif()

# set(CMAKE_Fortran_FORMAT FREE)

# Print out the discovered paths
cmake_print_variables(PYTHON_EXECUTABLE)
cmake_print_variables(PYTHON_INCLUDE_DIRS)
cmake_print_variables(F2PY_INCLUDE_DIR)
cmake_print_variables(NUMPY_INCLUDE_DIRS)

include_directories(
  ${PYTHON_INCLUDE_DIRS}
  ${NUMPY_INCLUDE_DIRS}
  ${F2PY_INCLUDE_DIR}
)
add_compile_definitions(NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

if (UNIX)
  add_compile_options(
    -O3 -Wno-uninitialized -fPIC
    $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy>
    $<$<COMPILE_LANGUAGE:Fortran>:-fno-second-underscore>
  )
  if (APPLE)
    add_compile_options(
      $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>)
  else()
    add_compile_options(
      $<$<COMPILE_LANGUAGE:Fortran>:-Wno-argument-mismatch>)
  endif()
else() # intel
  # TODO
  # set(FLAGS -fast -fpe0)
  # set(FLAGSF90 ${FLAGS} -ffree-form -Wobsolescent -fno-second-underscore)
endif()

set(impy_logging impy_openlogfile impy_closelogfile)
set(logging_source ${CMAKE_SOURCE_DIR}/src/logging.f)
set(ran_source ${CMAKE_SOURCE_DIR}/src/rangen.f)

### eposlhc
file(GLOB eposlhc_sources ${CMAKE_SOURCE_DIR}/src/EPOS-LHC/sources/*.f)
list(FILTER eposlhc_sources EXCLUDE REGEX epos_example\.f)

f2py_add_module(eposlhc
  FUNCTIONS
  aaset ainit aepos afinal hepmcstore
  getcharge idtrafo initializeepos initeposevt
  setstable setunstable xsection ${impy_logging}
  INTERFACE_SOURCES
  ${logging_source}
  ${CMAKE_SOURCE_DIR}/src/EPOS-LHC/sources/epos-bas-lhc.f
  ${CMAKE_SOURCE_DIR}/src/EPOS-LHC/sources/epos-ids-lhc.f
  ${CMAKE_SOURCE_DIR}/src/EPOS-LHC/sources/epos_interface.f
  SOURCES
  ${logging_source}
  ${eposlhc_sources}
)

### sib21
f2py_add_module(sib21
  FUNCTIONS
  sibyll sibyll_ini sib_sigma_hp sib_sigma_hair sib_sigma_hnuc
  int_nuc decsib decpar sibini sib_list isib_pid2pdg
  isib_pdg2pid pdg_ini ${impy_logging} sibhep1
  SOURCES
  ${logging_source}
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_21.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sib21aux.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
)

add_library(sib21rnd OBJECT ${ran_source})
target_compile_definitions(sib21rnd PRIVATE SIBYLL_SP)
target_compile_options(sib21rnd PRIVATE -cpp)
target_link_libraries(sib21 PRIVATE sib21rnd)
