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

# We need to overwrite the default cmake flags for Fortran to remove
# -DNDEBUG, since NDEBUG is a variable used in SIBYLL and the Ninja
# generator automatically runs the preprocessor on all Fortran files
set(CMAKE_Fortran_FLAGS_RELEASE -O3)

# CMake sets -fPIC etc. automatically, add only unusual options
if (UNIX)
  add_compile_options(
    -Wno-uninitialized
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
set(rangen_source_needs_preprocessing ${CMAKE_SOURCE_DIR}/src/rangen.fpp)
set(rangen_source rangen.f)

set_property(SOURCE rangen_source_needs_preprocessing
  PROPERTY
  Fortran_PREPROCESS True)

# WARNING this is not platform-independent and needs to be adjusted
# when new compilers should be supported
add_custom_command(
  OUTPUT ${rangen_source}
  COMMAND ${CMAKE_Fortran_COMPILER}
    -E ${rangen_source_needs_preprocessing}
    -o ${rangen_source}
)

### eposlhc
file(GLOB eposlhc_sources
  ${CMAKE_SOURCE_DIR}/src/EPOS-LHC/sources/*.f
  ${CMAKE_SOURCE_DIR}/src/EPOS-LHC/eposlhc.pyf  
)
list(FILTER eposlhc_sources EXCLUDE REGEX epos_example\.f)

f2py_add_module(_eposlhc
  FUNCTIONS
  aaset ainit aepos afinal hepmcstore
  getcharge idtrafo initializeepos initeposevt
  setstable setunstable xsection ${impy_logging}
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/EPOS-LHC/eposlhc.pyf
  ${eposlhc_sources}
  ${logging_source}
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
  ${CMAKE_SOURCE_DIR}/src/sibyll/sib21.pyf
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_21.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sib21aux.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
  ${logging_source}
)
add_library(sib21rnd OBJECT ${rangen_source_needs_preprocessing})
target_compile_definitions(sib21rnd PRIVATE SIBYLL_SP)
target_link_libraries(_sib21 PRIVATE sib21rnd)

### sib23
f2py_add_module(_sib23
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/sibyll/sib23.pyf
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll2.3.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_source}
)

### sib23c00
f2py_add_module(_sib23c00
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/sibyll/sib23c00.pyf
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll2.3c00.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_source}
)

### sib23c01
f2py_add_module(_sib23c01
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/sibyll/sib23c01.pyf
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll2.3c01.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_source}
)

### sib23c02
f2py_add_module(_sib23c02
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/sibyll/sib23c02.pyf
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll2.3c02.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_source}
)

### sib23c03
f2py_add_module(_sib23c03
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/sibyll/sib23c03.pyf
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll2.3c03.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_source}
)

### sib23d
f2py_add_module(_sib23d
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/sibyll/sib23d.pyf
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll2.3d.f
  ${CMAKE_SOURCE_DIR}/src/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_source}
)

### qgs01
f2py_add_module(_qgs01
  FUNCTIONS
  cqgsini sectnu xxaini psconf xxreg psaini chepevt
  ${impy_logging}
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/qgsjet/qgs01.pyf
  ${CMAKE_SOURCE_DIR}/src/qgsjet/qgsjet01d.f
  ${CMAKE_SOURCE_DIR}/src/qgsjet/impy_qgs1.f
  ${logging_source}
  ${rangen_source}
)

### qgsII03
set(QGSJETII_COMMON_FUNCTIONS
  cqgsini qgsect qgini qgconf qgreg chepevt ${impy_logging})

f2py_add_module(_qgsII03
  FUNCTIONS
  ${QGSJETII_COMMON_FUNCTIONS}
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/qgsjet/qgsII-03.pyf
  ${CMAKE_SOURCE_DIR}/src/qgsjet/qgsjet-II-03.f
  ${CMAKE_SOURCE_DIR}/src/qgsjet/impy_qgsII.f
  ${logging_source}
  ${rangen_source}
)

### qgsII04
f2py_add_module(_qgsII04
  FUNCTIONS
  ${QGSJETII_COMMON_FUNCTIONS}
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/qgsjet/qgsII-04.pyf
  ${CMAKE_SOURCE_DIR}/src/qgsjet/qgsjet-II-04.f
  ${CMAKE_SOURCE_DIR}/src/qgsjet/impy_qgsII.f
  ${logging_source}
  ${rangen_source}
)

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
  ${rangen_source}
)
# urqmd34 contains files which need preprocessing
set_source_files_properties(
  ${urqmd34_sources}
  PROPERTIES
  Fortran_PREPROCESS True
)
f2py_add_module(_urqmd34
  FUNCTIONS
  urqmd init uinit set0 params uounit strini loginit
  loadwtab norm_init output cascinit nucrad urqini partname
  chepevt ptsigtot init_rmmard ${impy_logging}
  SOURCES
  # Using pre-generated source results in this error
  #  Undefined symbols for architecture x86_64:
  #  "_f2pyinitterms_", referenced from:
  #  _f2py_init_terms in _urqmd34module.c.o
  # ${CMAKE_SOURCE_DIR}/src/urqmd-3.4/urqmd34.pyf
  ${urqmd34_sources}
)

### pythia6
f2py_add_module(_pythia6
  FUNCTIONS
  pyinit pyexec pytune pylist pyevnt pyevnw pystat pyedit
  pyhepc pychge pycomp pyk ${impy_logging}
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/pythia6/pythia6.pyf
  ${CMAKE_SOURCE_DIR}/src/pythia6/pythia-6.4.28.f
  ${logging_source}
)

### sophia
f2py_add_module(_sophia
  FUNCTIONS
  eventgen print_event crossection initial
  icon_pdg_sib init_rmmard toevt ${impy_logging}
  SOURCES
  ${CMAKE_SOURCE_DIR}/src/sophia/sophia.pyf
  ${CMAKE_SOURCE_DIR}/src/sophia/SOPHIA20.f
  ${CMAKE_SOURCE_DIR}/src/sophia/eventgen.f
  ${CMAKE_SOURCE_DIR}/src/sophia/sampling.f
  ${CMAKE_SOURCE_DIR}/src/sophia/inpoutput.f
  ${CMAKE_SOURCE_DIR}/src/sophia/jetset74dp.f
  ${CMAKE_SOURCE_DIR}/src/sophia/impy_sophia.f
  ${logging_source}
  ${rangen_source}
)
