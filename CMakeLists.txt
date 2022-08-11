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

add_library(rangen OBJECT ${CMAKE_SOURCE_DIR}/src/rangen.f)
target_compile_options(rangen PRIVATE -cpp)

# add_library(logging OBJECT ${logging_source})

### eposlhc
file(GLOB eposlhc_sources ${CMAKE_SOURCE_DIR}/src/EPOS-LHC/sources/*.f)
list(FILTER eposlhc_sources EXCLUDE REGEX epos_example\.f)

f2py_add_module(_eposlhc
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
set(SIBYLL_COMMON_FUNCTIONS
  sibyll sibyll_ini sib_sigma_hp sib_sigma_hair sib_sigma_hnuc
  int_nuc decsib decpar sibini sib_list isib_pid2pdg
  isib_pdg2pid pdg_ini ${impy_logging})

f2py_add_module(_sib21
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep1
  SOURCES
  ${logging_source}
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_21.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sib21aux.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
)

add_library(sib21rnd OBJECT ${CMAKE_SOURCE_DIR}/src/rangen.f)
target_compile_definitions(sib21rnd PRIVATE SIBYLL_SP)
target_compile_options(sib21rnd PRIVATE -cpp)
target_link_libraries(_sib21 PRIVATE sib21rnd)

### sib23d
f2py_add_module(_sib23d
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${logging_source}
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll2.3d.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
)
target_link_libraries(_sib23d PRIVATE rangen)

### qgs01
f2py_add_module(_qgs01
  FUNCTIONS
  cqgsini sectnu xxaini psconf xxreg psaini chepevt
  ${impy_logging}
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/qgsjet/qgsjet01d.f
  ${CMAKE_SOURCE_DIR}/src/qgsjet/impy_qgs1.f
  ${logging_source}
)
target_link_libraries(_qgs01 PRIVATE rangen)

### qgsII03
set(QGSJETII_COMMON_FUNCTIONS
  cqgsini qgsect qgini qgconf qgreg chepevt ${impy_logging})

f2py_add_module(_qgsII03
  FUNCTIONS
  ${QGSJETII_COMMON_FUNCTIONS}
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/qgsjet/qgsjet-II-03.f
  ${CMAKE_SOURCE_DIR}/src/qgsjet/impy_qgsII.f
  ${logging_source}
)
target_link_libraries(_qgsII03 PRIVATE rangen)

### qgsII04
f2py_add_module(_qgsII04
  FUNCTIONS
  ${QGSJETII_COMMON_FUNCTIONS}
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/qgsjet/qgsjet-II-04.f
  ${CMAKE_SOURCE_DIR}/src/qgsjet/impy_qgsII.f
  ${logging_source}
)
target_link_libraries(_qgsII04 PRIVATE rangen)

### urqmd34
# FIXME does not compile
# urqmd34module.c:734:13: error: redefinition of 'doc_f2py_rout_urqmd34_init'
# static char doc_f2py_rout_urqmd34_init[] = "\
# urqmd34module.c:606:13: note: previous definition is here
# static char doc_f2py_rout_urqmd34_init[] = "\
# urqmd34module.c:738:18: error: redefinition of 'f2py_rout_urqmd34_init'
# static PyObject *f2py_rout_urqmd34_init(const PyObject *capi_self,
# urqmd34module.c:610:18: note: previous definition is here
# static PyObject *f2py_rout_urqmd34_init(const PyObject *capi_self,
file(GLOB urqmd34_sources ${CMAKE_SOURCE_DIR}/src/urqmd-3.4/sources/*.f)
list(APPEND urqmd34_sources
  ${CMAKE_SOURCE_DIR}/src/urqmd-3.4/sources/CFmax.f90
  ${CMAKE_SOURCE_DIR}/src/urqmd-3.4/sources/quadri.f90
  ${CMAKE_SOURCE_DIR}/src/urqmd-3.4/sources/cornelius.f90
  ${logging_source}
  ${CMAKE_SOURCE_DIR}/src/urqmd-3.4/impy_urqmd.f
)
f2py_add_module(_urqmd34
  FUNCTIONS
  urqmd init uinit set0 params uounit strini loginit
  loadwtab norm_init output cascinit nucrad urqini partname
  chepevt ptsigtot init_rmmard ${impy_logging}
  SOURCES
  ${urqmd34_sources}
)
target_link_libraries(_urqmd34 PRIVATE rangen)

### pythia6
f2py_add_module(_pythia6
  FUNCTIONS
  pyinit pyexec pytune pylist pyevnt pyevnw pystat pyedit
  pyhepc pychge pycomp pyk ${impy_logging}
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/pythia6/pythia-6.4.28.f
  ${logging_source}
)
