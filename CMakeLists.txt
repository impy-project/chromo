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

### common for all models
set(impy_logging impy_openlogfile impy_closelogfile)
set(logging_source ${CMAKE_SOURCE_DIR}/src/logging.f)
set(rangen_source ${CMAKE_SOURCE_DIR}/src/rangen.f)

add_library(rangen OBJECT ${rangen_source})
target_compile_options(rangen PRIVATE -cpp)

add_library(logging OBJECT ${logging_source})

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

### sib23c00
f2py_add_module(_sib23c00
FUNCTIONS
${SIBYLL_COMMON_FUNCTIONS} sibhep3
SOURCES
${logging_source}
${CMAKE_SOURCE_DIR}/src/sibyll/sibyll2.3c00.f
${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
)
target_link_libraries(_sib23c00 PRIVATE rangen)

### sib23c01
f2py_add_module(_sib23c01
FUNCTIONS
${SIBYLL_COMMON_FUNCTIONS} sibhep3
SOURCES
${logging_source}
${CMAKE_SOURCE_DIR}/src/sibyll/sibyll2.3c01.f
${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
)
target_link_libraries(_sib23c01 PRIVATE rangen)

### sib23c02
f2py_add_module(_sib23c02
FUNCTIONS
${SIBYLL_COMMON_FUNCTIONS} sibhep3
SOURCES
${logging_source}
${CMAKE_SOURCE_DIR}/src/sibyll/sibyll2.3c02.f
${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
)
target_link_libraries(_sib23c02 PRIVATE rangen)

### sib23c03
f2py_add_module(_sib23c03
FUNCTIONS
${SIBYLL_COMMON_FUNCTIONS} sibhep3
SOURCES
${logging_source}
${CMAKE_SOURCE_DIR}/src/sibyll/sibyll2.3c03.f
${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
)
target_link_libraries(_sib23c03 PRIVATE rangen)

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
set(urqmd34_sources
  1fluid bessel delpart getmass hepcmp iso numrec pythia6409
  siglookup upmerge addpart blockres coload detbal getspin
  hepnam ityp2pdg output string urqmd angdis boxprg dwidth
  init jdecay2 paulibl saveinfo tabinit whichres anndec
  cascinit dectim error hepchg input make22 proppot scatter
  uhmerge urqinit
)
list(TRANSFORM urqmd34_sources APPEND .f)
list(APPEND urqmd34_sources
  CFmax.f90 quadri.f90 cornelius.f90
)
list(TRANSFORM urqmd34_sources PREPEND ${CMAKE_SOURCE_DIR}/src/urqmd-3.4/sources/)
list(APPEND urqmd34_sources
  ${CMAKE_SOURCE_DIR}/src/urqmd-3.4/impy_urqmd.f
  ${logging_source}
)
f2py_add_module(_urqmd34
  FUNCTIONS
  urqmd init uinit set0 params uounit strini loginit
  loadwtab norm_init output cascinit nucrad urqini partname
  chepevt ptsigtot init_rmmard ${impy_logging}
  INTERFACE_SOURCES
  ${urqmd34_sources}
  ${rangen_source}
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
