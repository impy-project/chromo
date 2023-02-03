cmake_minimum_required(VERSION 3.17)

project(chromo LANGUAGES Fortran)

# Find C/CXX compilers that match the Fortran compiler
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Compiler version is something like XY.Z.N, we want only XY
  string(REGEX MATCH "[0-9]+" COMPILER_MAJOR_VERSION ${CMAKE_Fortran_COMPILER_VERSION})
  find_program(CMAKE_C_COMPILER "gcc-${COMPILER_MAJOR_VERSION}")
  find_program(CMAKE_CXX_COMPILER "g++-${COMPILER_MAJOR_VERSION}")
else()
  message(FATAL_ERROR "Unknown Fortran compiler, cannot search for C/CXX compiler")
endif()

# Delayed toolchain initialization
enable_language(C CXX)

set(BUILD_dev_dpmjetIII193 CACHE
  PATH "source directory for dev_dpmjetIII193")
option(BUILD_sib23c00 OFF)
option(BUILD_sib23c02 OFF)
option(BUILD_sib23c03 OFF)

add_subdirectory(${CMAKE_SOURCE_DIR}/src/cpp/pybind11)

set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src/cpp/pybind11/tools)
include(CMakePrintHelpers)
include(F2Py)
find_package(pybind11 CONFIG)

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

macro(use_alternate_linker_if_available linker)
  find_program(LINKER_EXECUTABLE ld.${linker})
  if(LINKER_EXECUTABLE)
    message(STATUS "Linker detected: ${linker} (${LINKER_EXECUTABLE})")
    if(CLANG AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 12.0.0)
      message(STATUS "Clang links with ${LINKER_EXECUTABLE}")
      add_link_options(-ld-path=${linker})
    elseif(GCC AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER_EQUAL 12.0.0)
      message(STATUS "GCC links with ${LINKER_EXECUTABLE}")
      add_link_options(-ld-path=${linker})
    endif()
  endif()
endmacro()

use_alternate_linker_if_available(mold)

# Print out the discovered paths
cmake_print_variables(CMAKE_BUILD_TYPE)
cmake_print_variables(PYTHON_EXECUTABLE)
cmake_print_variables(PYTHON_INCLUDE_DIRS)
cmake_print_variables(PYTHON_LIBRARIES)
cmake_print_variables(F2PY_INCLUDE_DIR)
cmake_print_variables(NUMPY_INCLUDE_DIRS)
cmake_print_variables(BUILD_sib23c00)
cmake_print_variables(BUILD_sib23c02)
cmake_print_variables(BUILD_sib23c03)
cmake_print_variables(BUILD_dev_dpmjetIII193)

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
else() # Windows MinGW!
	add_compile_options(
    -Wno-uninitialized
    -w
    $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy>
    $<$<COMPILE_LANGUAGE:Fortran>:-fno-second-underscore>
	  $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>
    $<$<COMPILE_LANGUAGE:CXX>:-fpermissive>
  )
  # TODO
  # set(FLAGS -fast -fpe0)
  # set(FLAGSF90 ${FLAGS} -ffree-form -Wobsolescent -fno-second-underscore)
endif()

### common for all models
set(chromo_definitions CHROMO)
set(fortran_dir ${CMAKE_SOURCE_DIR}/src/fortran)
set(cpp_dir ${CMAKE_SOURCE_DIR}/src/cpp)
set(f2py_dir ${CMAKE_SOURCE_DIR}/src/f2py)
# Copy fortranobject.c to get rid of deep folder nesting
configure_file(${F2PY_INCLUDE_DIR}/fortranobject.c fortranobject.c COPYONLY)
set(f2py_source fortranobject.c)

set(chromo_functions chromo_openlogfile chromo_closelogfile npyrng)
set(logging_source ${fortran_dir}/logging.f)
set(rangen_source ${fortran_dir}/rangen.fpp ${fortran_dir}/rangen.c ${fortran_dir}/normal.c)


### eposlhc
file(GLOB eposlhc_sources ${fortran_dir}/epos/sources/*.f)
list(FILTER eposlhc_sources EXCLUDE REGEX epos_example\.f)
list(FILTER eposlhc_sources EXCLUDE REGEX epos-random\.f)
list(APPEND eposlhc_sources ${fortran_dir}/epos/epos-random-dummy.f)

f2py_add_module(_eposlhc
  FUNCTIONS
  aaset ainit aepos afinal hepmcstore
  getcharge idtrafo initepos initeposevt
  xsection rmmard ${chromo_functions}
  SOURCES
  ${eposlhc_sources}
  ${logging_source}
  ${rangen_source}
  INTERFACE_SOURCES
  ${fortran_dir}/epos/sources/epos-bas-lhc.f
  ${fortran_dir}/epos/sources/epos-ids-lhc.f
  ${fortran_dir}/epos/sources/epos_interface.f
  ${fortran_dir}/epos/sources/epos-random.f
  ${rangen_source}
  ${logging_source}
  INCLUDE_DIRS
  ${fortran_dir}/epos/sources ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)


### sib21
set(SIBYLL_COMMON_FUNCTIONS
  sibyll sibyll_ini sib_sigma_hp sib_sigma_hair sib_sigma_hnuc
  int_nuc decsib decpar sibini sibhep sib_list isib_pid2pdg
  isib_pdg2pid pdg_ini ${chromo_functions})

f2py_add_module(_sib21
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} spgasdev
  SOURCES
  ${fortran_dir}/sibyll/sibyll_21.f
  ${fortran_dir}/sibyll/sib21aux.f
  ${fortran_dir}/sibyll/sibyll_init.fpp
  ${logging_source}
  ${rangen_source}
  INCLUDE_DIRS
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
  SIBYLL_21
)

### sib23
f2py_add_module(_sib23
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} gasdev
  SOURCES
  ${fortran_dir}/sibyll/sibyll2.3.f
  ${fortran_dir}/sibyll/sibyll_init.fpp
  ${logging_source}
  ${rangen_source}
  INCLUDE_DIRS
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)


### sib23c01
f2py_add_module(_sib23c01
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} gasdev
  SOURCES
  ${fortran_dir}/sibyll/sibyll2.3c01.f
  ${fortran_dir}/sibyll/sibyll_init.fpp
  ${logging_source}
  ${rangen_source}
  INCLUDE_DIRS
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)


### sib23d
f2py_add_module(_sib23d
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} gasdev
  SOURCES
  ${fortran_dir}/sibyll/sibyll2.3d.f
  ${fortran_dir}/sibyll/sibyll_init.fpp
  ${logging_source}
  ${rangen_source}
  INCLUDE_DIRS
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)

### qgs01
f2py_add_module(_qgs01
  FUNCTIONS
  cqgsini sectnu xxaini psconf xxreg psaini chepevt
  ${chromo_functions}
  SOURCES
  ${fortran_dir}/qgsjet/qgsjet01d.f
  ${fortran_dir}/qgsjet/chromo_qgs1.f
  ${logging_source}
  ${rangen_source}
  INCLUDE_DIRS
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)

### qgsII03
set(QGSJETII_COMMON_FUNCTIONS
  cqgsini qgsect qgini qgconf qgreg chepevt ${chromo_functions})

f2py_add_module(_qgsII03
  FUNCTIONS
  ${QGSJETII_COMMON_FUNCTIONS}
  SOURCES
  ${fortran_dir}/qgsjet/qgsjet-II-03.f
  ${fortran_dir}/qgsjet/chromo_qgsII.f
  ${logging_source}
  ${rangen_source}
  INCLUDE_DIRS
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)

### qgsII04
f2py_add_module(_qgsII04
  FUNCTIONS
  ${QGSJETII_COMMON_FUNCTIONS}
  SOURCES
  ${fortran_dir}/qgsjet/qgsjet-II-04.f
  ${fortran_dir}/qgsjet/chromo_qgsII.f
  ${logging_source}
  ${rangen_source}
  INCLUDE_DIRS
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
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
list(TRANSFORM urqmd34_sources PREPEND ${fortran_dir}/urqmd-3.4/sources/)
list(APPEND urqmd34_sources
  ${fortran_dir}/urqmd-3.4/chromo_urqmd.f
  ${logging_source}
  ${rangen_source}
)
f2py_add_module(_urqmd34
  FUNCTIONS
  urqmd init uinit set0 params uounit strini loginit
  loadwtab norm_init output cascinit nucrad urqini partname
  chepevt ptsigtot ${chromo_functions}
  SOURCES
  ${urqmd34_sources}
  INCLUDE_DIRS
  ${fortran_dir}/urqmd-3.4/sources ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)

### pythia6
f2py_add_module(_pythia6
  FUNCTIONS
  pyinit pyexec pytune pylist pyevnt pyevnw pystat pyedit
  pyhepc pychge pycomp pyk ${chromo_functions}
  SOURCES
  ${fortran_dir}/pythia6/pythia-6.4.28.f
  ${logging_source}
  ${rangen_source}
  INCLUDE_DIRS
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)

### pythia8
# About the license and how it affects chromo.
# Pythia-8 is licensed under GPL-v2 or higher, which means GPL-v3
# applies. An excerpt from GPL-v3:

# > A compilation of a covered work with other separate and independent
# > works, which are not by their nature extensions of the covered work,
# > and which are not combined with it such as to form a larger program,
# > in or on a volume of a storage or distribution medium, is called an
# > "aggregate" if the compilation and its resulting copyright are not
# > used to limit the access or legal rights of the compilation's users
# > beyond what the individual works permit.  Inclusion of a covered work
# > in an aggregate does not cause this License to apply to the other
# > parts of the aggregate.
#
# Our legal position is that chromo is an "aggregate". It is not adding
# functionality to Pythia-8 nor is it modifying the original source. The
# Pythia-8 library is redistributed by us under its original license, which is
# allowed. chromo is merely a way to steer Pythia-8 and to provide its output in a
# particular format. The Pythia-8 authors are free to take some of our code and
# make it GPL code, they can do that anyway, but we are still free to
# distribute the code that we write under a less restrictive license.
#
# To strengthen this position, we compile the original Pythia-8 code into
# a shared library to which we only link. This makes clear where the boundary
# between chromo and Pythia-8 is.
file(GLOB pythia8_sources ${cpp_dir}/pythia83/src/*.cc)

add_library(libpythia8 SHARED ${pythia8_sources})
set_target_properties(libpythia8 PROPERTIES OUTPUT_NAME pythia8)
target_include_directories(libpythia8 PUBLIC ${cpp_dir}/pythia83/include)
# we don't use this, we must set xmldir at runtime
target_compile_definitions(libpythia8 PRIVATE XMLDIR="")

pybind11_add_module(_pythia8 ${cpp_dir}/_pythia8.cpp)
target_link_libraries(_pythia8 PRIVATE libpythia8)
install(TARGETS libpythia8 DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})


# This part is needed when we figure out the problem with
# building Pythia8 on Windows
# if (WIN32)
# find_package(dlfcn-win32 REQUIRED)
# set(CMAKE_DL_LIBS dlfcn-win32::dl)
# target_link_libraries(libpythia8 PRIVATE ${CMAKE_DL_LIBS})
# target_link_libraries(_pythia8 PRIVATE ${CMAKE_DL_LIBS})
# target_link_libraries(_pythia8 PRIVATE "-static")
# endif()

### sophia
f2py_add_module(_sophia
  FUNCTIONS
  eventgen print_event crossection initial
  icon_pdg_sib toevt ${chromo_functions}
  SOURCES
  ${fortran_dir}/sophia/SOPHIA20.f
  ${fortran_dir}/sophia/eventgen.f
  ${fortran_dir}/sophia/sampling.f
  ${fortran_dir}/sophia/inpoutput.f
  ${fortran_dir}/sophia/jetset74dp.f
  ${fortran_dir}/sophia/chromo_sophia.f
  ${logging_source}
  ${rangen_source}
  INCLUDE_DIRS
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)

### dpmjet306
set(dpmjet306_sources
  ${fortran_dir}/dpmjet3.0-6/sources/dpmjet3.0-6.f
  ${fortran_dir}/dpmjet3.0-6/sources/phojet1.12-35c4.f
  ${fortran_dir}/dpmjet3.0-6/sources/pythia6115dpm3v1.f
  ${fortran_dir}/dpmjet3.0-6/sources/user3.0-6.f
)

f2py_add_module(_dpmjet306
  FUNCTIONS
  pho_event dt_init dt_kkinc idt_icihad dt_xsglau
  pycomp dt_initjs dt_inucas dt_evtout idt_ipdgha
  pho_init pho_setpar pho_pname pho_pmass pho_setmdl
  pho_setpdf pycomp pho_xsect pho_borncs pho_harmci pho_fitout
  pho_mcini pho_ptcut pytune pho_rregpar pho_sregpar pho_prevnt
  ipho_pdg2id ipho_id2pdg pho_harint pho_harxto pho_harxpt
  pho_setpcomb dt_phoxs dt_xshn dt_flahad
  ${chromo_functions}
  SOURCES
  ${dpmjet306_sources}
  ${rangen_source}
  ${logging_source}
  INCLUDE_DIRS
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)


### Phojet112
f2py_add_module(_phojet112
  FUNCTIONS
  pho_event pycomp pho_init pho_setpar pho_pname pho_pmass
  pho_setmdl pho_setpdf pycomp pho_xsect pho_borncs pho_harmci
  pho_fitout pho_mcini pho_ptcut pytune pho_rregpar pho_sregpar
  pho_prevnt ipho_pdg2id ipho_id2pdg pho_harint pho_harxto pho_harxpt
  pho_setpcomb ${chromo_functions}
  SOURCES
  ${dpmjet306_sources}
  ${logging_source}
  ${rangen_source}
  INCLUDE_DIRS
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)


### dpmjetIII191
file(GLOB dpmjetIII191_sources
  ${fortran_dir}/dpmjetIII-19.1/src/phojet/*.f
  ${fortran_dir}/dpmjetIII-19.1/src/pythia/*.f
  ${fortran_dir}/dpmjetIII-19.1/src/dpmjet/*.f
)

list(FILTER dpmjetIII191_sources EXCLUDE REGEX DT_RNDM\.f)
list(FILTER dpmjetIII191_sources EXCLUDE REGEX DT_RNDMST\.f)
list(FILTER dpmjetIII191_sources EXCLUDE REGEX DT_RNDMTE\.f)
list(FILTER dpmjetIII191_sources EXCLUDE REGEX PYR\.f)
# cmake does not like filenames with parantheses, but some Fortran files
# include those. As a workaround, we generate modified source files and
# renamed copies of the headers. Not elegant, but works. 🤷‍♂️


### dpmjetIII191
file(GLOB dpmjetIII191_sources
  ${fortran_dir}/dpmjetIII-19.1/src/phojet/*.f
  ${fortran_dir}/dpmjetIII-19.1/src/pythia/*.f
  ${fortran_dir}/dpmjetIII-19.1/src/dpmjet/*.f
)

list(FILTER dpmjetIII191_sources EXCLUDE REGEX DT_RNDM\.f)
list(FILTER dpmjetIII191_sources EXCLUDE REGEX DT_RNDMST\.f)
list(FILTER dpmjetIII191_sources EXCLUDE REGEX DT_RNDMTE\.f)
list(FILTER dpmjetIII191_sources EXCLUDE REGEX PYR\.f)
# cmake does not like filenames with parantheses, but some Fortran files
# include those. As a workaround, we generate modified source files and
# renamed copies of the headers. Not elegant, but works. 🤷‍♂️


# Workaround to workaround :) Windows don't like long lists
# User temporary file to keep long lists
set(temp_dpmjet_file temp.dpmjetIII191)
string(REPLACE ";"  " " temp_dpmjet_str "${dpmjetIII191_sources}")
file(WRITE ${temp_dpmjet_file} "${temp_dpmjet_str}")

execute_process(
  COMMAND
  ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/dpmjetIII191_workaround.py
  ${CMAKE_CURRENT_BINARY_DIR}
  ${fortran_dir}/dpmjetIII-19.1/include
  ${CMAKE_SOURCE_DIR}/${temp_dpmjet_file}
)
file(READ ${temp_dpmjet_file} dpmjetIII191_modded_sources)
file(REMOVE ${temp_dpmjet_file})


f2py_add_module(_dpmjetIII191
  FUNCTIONS
  pho_event dt_init dt_kkinc
  idt_icihad dt_xsglau pycomp dt_initjs dt_inucas idt_ipdgha dt_evtout
  pho_init pho_setpar poevt1 poevt2 pho_pname pho_pmass pho_setmdl
  pho_setpdf pycomp pho_xsect pho_borncs pho_harmci pho_fitout pho_mcini pho_ptcut
  pytune pho_rregpar pho_sregpar pho_prevnt ipho_pdg2id ipho_id2pdg pho_harint
  pho_harxto pho_harxpt pho_setpcomb dt_phoxs dt_xshn dt_flahad
  dt_title dt_ficonf pho_ghhias ${chromo_functions}
  SOURCES
  ${dpmjetIII191_modded_sources}
  ${fortran_dir}/dpmjetIII-19.1/common/dummies.f
  ${logging_source}
  ${rangen_source}
  INCLUDE_DIRS
  ${fortran_dir}/dpmjetIII-19.1/include/phojet
  ${fortran_dir}/dpmjetIII-19.1/include/dpmjet
  ${fortran_dir}/dpmjetIII-19.1/include/pythia
  ${CMAKE_CURRENT_BINARY_DIR}/include
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)


### phojet191
f2py_add_module(_phojet191
  FUNCTIONS
  pho_event pycomp pho_init pho_setpar poevt1
  poevt2 pho_pname pho_pmass pho_setmdl pho_setpdf
  pho_xsect pho_borncs pho_harmci pho_fitout pho_mcini
  pho_ptcut pytune pho_rregpar pho_sregpar pho_prevnt
  ipho_pdg2id ipho_id2pdg pho_harint pho_harxto pho_harxpt
  pho_setpcomb dt_ficonf pho_ghhias ${chromo_functions}
  SOURCES
  ${dpmjetIII191_modded_sources}
  ${fortran_dir}/dpmjetIII-19.1/common/dummies.f
  ${logging_source}
  ${rangen_source}
  INCLUDE_DIRS
  ${fortran_dir}/dpmjetIII-19.1/include/phojet
  ${fortran_dir}/dpmjetIII-19.1/include/dpmjet
  ${fortran_dir}/dpmjetIII-19.1/include/pythia
  ${CMAKE_CURRENT_BINARY_DIR}/include
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)



### dpmjetIII193
file(GLOB dpmjetIII193_sources
  ${fortran_dir}/dpmjetIII-19.3/src/phojet/*.f
  ${fortran_dir}/dpmjetIII-19.3/src/pythia/*.f
  ${fortran_dir}/dpmjetIII-19.3/src/dpmjet/*.f
  ${fortran_dir}/dpmjetIII-19.3/common/*.f
)

list(FILTER dpmjetIII193_sources EXCLUDE REGEX DT_RNDM\.f)
list(FILTER dpmjetIII193_sources EXCLUDE REGEX DT_RNDMST\.f)
list(FILTER dpmjetIII193_sources EXCLUDE REGEX DT_RNDMTE\.f)
list(FILTER dpmjetIII193_sources EXCLUDE REGEX PYR\.f)


f2py_add_module(_dpmjetIII193
  FUNCTIONS
  pho_event dt_init dt_kkinc idt_icihad dt_xsglau pycomp dt_initjs
  dt_inucas idt_ipdgha dt_evtout pho_init pho_setpar poevt1 poevt2
  pho_pname pho_pmass pho_setmdl pho_setpdf pycomp pho_xsect pho_borncs
  pho_harmci pho_fitout pho_mcini pho_ptcut pytune pho_rregpar pho_sregpar
  pho_prevnt ipho_pdg2id ipho_id2pdg pho_harint pho_harxto pho_harxpt
  pho_setpcomb dt_phoxs dt_xshn dt_flahad dt_title pho_ghhias
  ${chromo_functions}
  SOURCES
  ${dpmjetIII193_sources}
  ${fortran_dir}/dpmjetIII-19.3/common/dummies.f
  ${rangen_source}
  ${logging_source}
  INCLUDE_DIRS
  ${fortran_dir}/dpmjetIII-19.3/include/phojet
  ${fortran_dir}/dpmjetIII-19.3/include/dpmjet
  ${fortran_dir}/dpmjetIII-19.3/include/pythia
  ${fortran_dir}/dpmjetIII-19.3/include/flinclude
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)


f2py_add_module(_phojet193
  FUNCTIONS
  pho_event pycomp pho_init pho_setpar poevt1 poevt2
  pho_pname pho_pmass pho_setmdl pho_setpdf pycomp pho_xsect pho_borncs
  pho_harmci pho_fitout pho_mcini pho_ptcut pytune pho_rregpar pho_sregpar
  pho_prevnt ipho_pdg2id ipho_id2pdg pho_harint pho_harxto pho_harxpt
  pho_setpcomb pho_ghhias
  ${chromo_functions}
  SOURCES
  ${dpmjetIII193_sources}
  ${fortran_dir}/dpmjetIII-19.3/common/dummies.f
  ${rangen_source}
  ${logging_source}
  INCLUDE_DIRS
  ${fortran_dir}/dpmjetIII-19.3/include/phojet
  ${fortran_dir}/dpmjetIII-19.3/include/dpmjet
  ${fortran_dir}/dpmjetIII-19.3/include/pythia
  ${fortran_dir}/dpmjetIII-19.3/include/flinclude
  ${NUMPY_INCLUDE_DIRS}
  COMPILE_DEFS
  ${chromo_definitions}
)


##### extra optional models follow below ######

if(BUILD_sib23c00)
  f2py_add_module(_sib23c00
    FUNCTIONS
    ${SIBYLL_COMMON_FUNCTIONS} gasdev
    SOURCES
    ${fortran_dir}/sibyll/sibyll2.3c00.f
    ${fortran_dir}/sibyll/sibyll_init.fpp
    ${logging_source}
    ${rangen_source}
    INCLUDE_DIRS
    ${NUMPY_INCLUDE_DIRS}
    COMPILE_DEFS
    ${chromo_definitions}
  )
endif()

if(BUILD_sib23c02)
  f2py_add_module(_sib23c02
    FUNCTIONS
    ${SIBYLL_COMMON_FUNCTIONS} gasdev
    SOURCES
    ${fortran_dir}/sibyll/sibyll2.3c02.f
    ${fortran_dir}/sibyll/sibyll_init.fpp
    ${logging_source}
    ${rangen_source}
    INCLUDE_DIRS
    ${NUMPY_INCLUDE_DIRS}
    COMPILE_DEFS
    ${chromo_definitions}
  )
endif()

if(BUILD_sib23c03)
  f2py_add_module(_sib23c03
    FUNCTIONS
    ${SIBYLL_COMMON_FUNCTIONS} gasdev
    SOURCES
    ${fortran_dir}/sibyll/sibyll2.3c03.f
    ${fortran_dir}/sibyll/sibyll_init.fpp
    ${logging_source}
    ${rangen_source}
    INCLUDE_DIRS
    ${NUMPY_INCLUDE_DIRS}
    COMPILE_DEFS
    ${chromo_definitions}
  )
endif()

if(BUILD_dev_dpmjetIII193)

  set(dev_source ${BUILD_dev_dpmjetIII193})

  file(GLOB dev_dpmjetIII193_sources
    ${dev_source}/src/phojet/*.f
    ${dev_source}/src/pythia/*.f
    ${dev_source}/src/dpmjet/*.f
    ${dev_source}/common/*.f
  )

  list(FILTER dev_dpmjetIII193_sources EXCLUDE REGEX DT_RNDM\.f)
  list(FILTER dev_dpmjetIII193_sources EXCLUDE REGEX DT_RNDMST\.f)
  list(FILTER dev_dpmjetIII193_sources EXCLUDE REGEX DT_RNDMTE\.f)
  list(FILTER dev_dpmjetIII193_sources EXCLUDE REGEX PYR\.f)

  f2py_add_module(_dev_dpmjetIII193
    FUNCTIONS
    pho_event dt_init dt_kkinc idt_icihad dt_xsglau pycomp dt_initjs
    dt_inucas idt_ipdgha dt_evtout pho_init pho_setpar poevt1 poevt2
    pho_pname pho_pmass pho_setmdl pho_setpdf pycomp pho_xsect pho_borncs
    pho_harmci pho_fitout pho_mcini pho_ptcut pytune pho_rregpar pho_sregpar
    pho_prevnt ipho_pdg2id ipho_id2pdg pho_harint pho_harxto pho_harxpt
    pho_setpcomb dt_phoxs dt_xshn dt_flahad dt_title pho_ghhias
    ${chromo_functions}
    SOURCES
    ${dev_dpmjetIII193_sources}
    ${dev_source}/common/dummies.f
    ${rangen_source}
    ${logging_source}
    INCLUDE_DIRS
    ${dev_source}/include/phojet
    ${dev_source}/include/dpmjet
    ${dev_source}/include/pythia
    ${dev_source}/include/flinclude
    ${NUMPY_INCLUDE_DIRS}
    COMPILE_DEFS
    ${chromo_definitions}
  )

endif()