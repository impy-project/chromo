cmake_minimum_required(VERSION 3.18)

project(impylib LANGUAGES C Fortran)

if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  message(
    WARNING
      "You should call this through setup.py so that all variables are set; e.g. python setup.py develop"
  )
endif()

find_package(Python 3.7 REQUIRED
  COMPONENTS Interpreter Development.Module NumPy)

# in Numpy-1.22+, this becomes easier: import numpy.f2py; print(numpy.f2py.get_include())
if (NOT F2PY_INCLUDE_DIR)
  execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import numpy.f2py; from pathlib import Path; print(Path(numpy.f2py.__file__).parent)"
    OUTPUT_VARIABLE _f2py_directory
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  set(F2PY_INCLUDE_DIR "${_f2py_directory}/src" CACHE STRING "F2PY source directory location" FORCE)
endif()

# Print out the discovered paths
include(CMakePrintHelpers)
cmake_print_variables(Python_EXECUTABLE)
cmake_print_variables(Python_INCLUDE_DIRS)
cmake_print_variables(F2PY_INCLUDE_DIR)
cmake_print_variables(Python_NumPy_INCLUDE_DIRS)

include_directories(${F2PY_INCLUDE_DIR})
add_compile_definitions(NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

if (UNIX)
  add_compile_options(-O3 -Wno-uninitialized -fPIC)
  add_compile_options(
    $<$<COMPILE_LANGUAGE:Fortran>:-fno-second-underscore>)
  if (APPLE)
    add_compile_options(
      $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>)
  endif()
else() # intel
  # TODO
  # set(FLAGS -fast -fpe0)
  # set(FLAGSF90 ${FLAGS} -ffree-form -Wobsolescent -fno-second-underscore)
endif()

set(IMPY_LOGGING impy_openlogfile impy_closelogfile)
set(logging_source ${CMAKE_SOURCE_DIR}/src/logging.f)

file(GLOB eposlhc_sources ${CMAKE_SOURCE_DIR}/src/EPOS-LHC/sources/*.f)
list(FILTER eposlhc_sources EXCLUDE REGEX epos_example\.f)
set(eposlhc_interface_sources
  ${CMAKE_SOURCE_DIR}/src/EPOS-LHC/sources/epos-bas-lhc.f
  ${CMAKE_SOURCE_DIR}/src/EPOS-LHC/sources/epos-ids-lhc.f
  ${CMAKE_SOURCE_DIR}/src/EPOS-LHC/sources/epos_interface.f
  ${logging_source}
)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/eposlhcmodule.c ${CMAKE_CURRENT_BINARY_DIR}/eposlhc-f2pywrappers.f
  COMMAND ${Python_EXECUTABLE} -m numpy.f2py
                  -m eposlhc
                  -h ${CMAKE_CURRENT_BINARY_DIR}/eposlhc.pyf
                  --overwrite-signature only: aaset ainit aepos afinal hepmcstore
                  getcharge idtrafo initializeepos initeposevt
                  setstable setunstable xsection ${IMPY_LOGGING} :
                  ${eposlhc_interface_sources}
                  > ${CMAKE_CURRENT_BINARY_DIR}/eposlhc.0.log 2>&1
  COMMAND ${Python_EXECUTABLE} -m numpy.f2py
                  ${CMAKE_CURRENT_BINARY_DIR}/eposlhc.pyf
                  > ${CMAKE_CURRENT_BINARY_DIR}/eposlhc.1.log 2>&1
  DEPENDS ${eposlhc_sources} ${logging_source} 
)

Python_add_library(eposlhc MODULE WITH_SOABI
  ${F2PY_INCLUDE_DIR}/fortranobject.c
  ${CMAKE_CURRENT_BINARY_DIR}/eposlhcmodule.c
  ${CMAKE_CURRENT_BINARY_DIR}/eposlhc-f2pywrappers.f
  ${eposlhc_sources}
  ${logging_source}
)
target_link_libraries(eposlhc PRIVATE Python::NumPy)
