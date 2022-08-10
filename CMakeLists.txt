cmake_minimum_required(VERSION 3.18)

project(impylib LANGUAGES C Fortran)

include(CMakePrintHelpers)

if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  message(
    WARNING
      "You should call this through setup.py so that all variables are set; e.g. python setup.py develop"
  )
endif()

# TODO this needs to be replaced, everything has to be derived from Python_EXECUTABLE
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
cmake_print_variables(Python_EXECUTABLE)
cmake_print_variables(Python_INCLUDE_DIRS)
cmake_print_variables(F2PY_INCLUDE_DIR)
cmake_print_variables(Python_NumPy_INCLUDE_DIRS)

include_directories(${F2PY_INCLUDE_DIR})
add_compile_definitions(NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

if (UNIX)
  add_compile_options(
    -O3 -Wno-uninitialized -fPIC
    # $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy>
    $<$<COMPILE_LANGUAGE:Fortran>:-fno-second-underscore>
  )

  if (APPLE)
    add_compile_options(
      $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>)
  endif()
else() # intel
  # TODO
  # set(FLAGS -fast -fpe0)
  # set(FLAGSF90 ${FLAGS} -ffree-form -Wobsolescent -fno-second-underscore)
endif()

set(impy_logging impy_openlogfile impy_closelogfile)
set(logging_source ${CMAKE_SOURCE_DIR}/src/logging.f)
set(ran_source ${CMAKE_SOURCE_DIR}/src/rangen.f)

# eposlhc
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
                  setstable setunstable xsection ${impy_logging} :
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

# sib21
set(sib21_sources
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_21.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sib21aux.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
  ${logging_source}
)
set(sib21_interface_sources ${sib21_sources})

add_custom_command(
  OUTPUT
  ${CMAKE_CURRENT_BINARY_DIR}/sib21module.c
  ${CMAKE_CURRENT_BINARY_DIR}/sib21-f2pywrappers.f
  COMMAND ${Python_EXECUTABLE} -m numpy.f2py
                  -m sib21
                  -h ${CMAKE_CURRENT_BINARY_DIR}/sib21.pyf
                  --overwrite-signature only: sibyll sibyll_ini
                  sib_sigma_hp sib_sigma_hair sib_sigma_hnuc
                  int_nuc decsib decpar sibini sib_list isib_pid2pdg
                  isib_pdg2pid pdg_ini ${impy_logging} sibhep1 :
                  ${sib21_interface_sources}
                  > ${CMAKE_CURRENT_BINARY_DIR}/sib21.0.log 2>&1
  COMMAND ${Python_EXECUTABLE} -m numpy.f2py
                  ${CMAKE_CURRENT_BINARY_DIR}/sib21.pyf
                  > ${CMAKE_CURRENT_BINARY_DIR}/sib21.1.log 2>&1
  DEPENDS ${sib21_sources} ${logging_source} 
)

add_library(sib21rnd OBJECT ${ran_source})
target_compile_definitions(sib21rnd PRIVATE SIBYLL_SP)
target_compile_options(sib21rnd PRIVATE -cpp)

Python_add_library(sib21 MODULE WITH_SOABI
  ${F2PY_INCLUDE_DIR}/fortranobject.c
  ${CMAKE_CURRENT_BINARY_DIR}/sib21module.c
  ${CMAKE_CURRENT_BINARY_DIR}/sib21-f2pywrappers.f
  ${sib21_sources}
  ${logging_source}
)
target_link_libraries(sib21 PRIVATE Python::NumPy sib21rnd)
