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
cmake_print_variables(PYTHON_LIBRARIES)
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
set(fortran_dir ${CMAKE_SOURCE_DIR}/src/fortran)
set(f2py_dir ${CMAKE_SOURCE_DIR}/src/f2py)

set(impy_logging impy_openlogfile impy_closelogfile)
set(logging_source ${fortran_dir}/logging.f)
set(rangen_source rangen.f)
set(rangen_fpp_source ${fortran_dir}/rangen.fpp)
set(rangen_sib21_source rangen_sib21.f)
set(sibyll_init_source sibyll_init.f)
set(sibyll_init_sib21_source sibyll_init_sib21.f)

#
set(sophia_jetset74dp_source jetset74dp.f)
set(sophia_eventgen_source eventgen.f)
set(pythia6_pythia6428_source pythia-6.4.28.f)
set(epos_random_source epos-random.f)

# WARNING this is not platform-independent and needs to be adjusted
# when new fortran compilers should be supported
add_custom_command(
  OUTPUT ${rangen_source}
  COMMAND ${CMAKE_Fortran_COMPILER}
    -E ${fortran_dir}/rangen.fpp
    -o ${rangen_source}
)
add_custom_command(
  OUTPUT ${rangen_sib21_source}
  COMMAND ${CMAKE_Fortran_COMPILER}
    -E ${fortran_dir}/rangen.fpp
    -o ${rangen_sib21_source}
    -DSIBYLL_SP
)
add_custom_command(
  OUTPUT ${sibyll_init_source}
  COMMAND ${CMAKE_Fortran_COMPILER}
    -E ${fortran_dir}/sibyll/sibyll_init.fpp
    -o ${sibyll_init_source}
    -DSIBYLL_TABLE_LENGTH=99
)
add_custom_command(
  OUTPUT ${sibyll_init_sib21_source}
  COMMAND ${CMAKE_Fortran_COMPILER}
    -E ${fortran_dir}/sibyll/sibyll_init.fpp
    -o ${sibyll_init_sib21_source}
    -DSIBYLL_TABLE_LENGTH=49
)


add_custom_command(
  OUTPUT ${pythia6_pythia6428_source}
  COMMAND ${CMAKE_Fortran_COMPILER}
    -E -cpp ${fortran_dir}/pythia6/pythia-6.4.28.f
    -o ${pythia6_pythia6428_source}
    -DIMPY
)

add_custom_command(
  OUTPUT ${sophia_jetset74dp_source}
  COMMAND ${CMAKE_Fortran_COMPILER}
    -E -cpp ${fortran_dir}/sophia/jetset74dp.f
    -o ${sophia_jetset74dp_source}
    -DIMPY
)

add_custom_command(
  OUTPUT ${sophia_eventgen_source}
  COMMAND ${CMAKE_Fortran_COMPILER}
    -E -cpp ${fortran_dir}/sophia/eventgen.f
    -o ${sophia_eventgen_source}
    -DIMPY
)


### eposlhc
file(GLOB eposlhc_sources ${fortran_dir}/EPOS-LHC/sources/*.f)
list(FILTER eposlhc_sources EXCLUDE REGEX epos_example\.f)
# list(FILTER eposlhc_sources EXCLUDE REGEX epos-random\.f)


add_custom_command(
  OUTPUT ${epos_random_source}
  COMMAND ${CMAKE_Fortran_COMPILER}
    -E -cpp ${fortran_dir}/EPOS-LHC/sources/epos-random.f
    -o ${epos_random_source}
    -DIMPY
  COMMAND cp  ${fortran_dir}/EPOS-LHC/sources/epos.inc epos.inc  
)

f2py_add_module(_eposlhc
  FUNCTIONS
  aaset ainit aepos afinal hepmcstore
  getcharge idtrafo initializeepos initeposevt
  setstable setunstable xsection init_rmmard ${impy_logging}
  SOURCES
  # ${f2py_dir}/_eposlhc.pyf
  # ${f2py_dir}/_eposlhcmodule.c
  # ${f2py_dir}/_eposlhc-f2pywrappers.f
  ${eposlhc_sources}
  ${logging_source}
  ${rangen_source}
  INTERFACE_SOURCES
  ${fortran_dir}/EPOS-LHC/sources/epos-bas-lhc.f 
  ${fortran_dir}/EPOS-LHC/sources/epos-ids-lhc.f 
  ${fortran_dir}/EPOS-LHC/sources/epos_interface.f
  ${fortran_dir}/EPOS-LHC/sources/epos-random.f
  ${rangen_source}
  ${logging_source}
  INCLUDE_DIRS
  ${fortran_dir}/EPOS-LHC/sources
  COMPILE_DEFS
  IMPY
)


### sib21
set(SIBYLL_COMMON_FUNCTIONS
  sibyll sibyll_ini sib_sigma_hp sib_sigma_hair sib_sigma_hnuc
  int_nuc decsib decpar sibini sib_list isib_pid2pdg
  isib_pdg2pid pdg_ini init_rmmard ${impy_logging})

f2py_add_module(_sib21
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep1
  SOURCES
  # ${f2py_dir}/_sib21.pyf
  # ${f2py_dir}/_sib21module.c
  # ${f2py_dir}/_sib21-f2pywrappers.f
  ${fortran_dir}/sibyll/sibyll_21.f
  ${fortran_dir}/sibyll/sib21aux.f
  ${logging_source}
  ${rangen_sib21_source}
  ${sibyll_init_sib21_source}
  COMPILE_DEFS
  IMPY
)

### sib23
f2py_add_module(_sib23
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  # ${f2py_dir}/_sib23.pyf
  # ${f2py_dir}/_sib23module.c
  # ${f2py_dir}/_sib23-f2pywrappers.f
  ${fortran_dir}/sibyll/sibyll2.3.f
  ${sibyll_init_source}
  ${logging_source}
  ${rangen_source}
  COMPILE_DEFS
  IMPY
)

### sib23c00
f2py_add_module(_sib23c00
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  # ${f2py_dir}/_sib23c00.pyf
  # ${f2py_dir}/_sib23c00module.c
  # ${f2py_dir}/_sib23c00-f2pywrappers.f
  ${fortran_dir}/sibyll/sibyll2.3c00.f
  ${sibyll_init_source}
  ${logging_source}
  ${rangen_source}
  COMPILE_DEFS
  IMPY
)

### sib23c01
f2py_add_module(_sib23c01
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  # ${f2py_dir}/_sib23c01.pyf
  # ${f2py_dir}/_sib23c01module.c
  # ${f2py_dir}/_sib23c01-f2pywrappers.f
  ${fortran_dir}/sibyll/sibyll2.3c01.f
  ${sibyll_init_source}
  ${logging_source}
  ${rangen_source}
  COMPILE_DEFS
  IMPY
)

### sib23c02
f2py_add_module(_sib23c02
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  # ${f2py_dir}/_sib23c02.pyf
  # ${f2py_dir}/_sib23c02module.c
  # ${f2py_dir}/_sib23c02-f2pywrappers.f
  ${fortran_dir}/sibyll/sibyll2.3c02.f
  ${sibyll_init_source}
  ${logging_source}
  ${rangen_source}
  COMPILE_DEFS
  IMPY
)

### sib23c03
f2py_add_module(_sib23c03
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  # ${f2py_dir}/_sib23c03.pyf
  # ${f2py_dir}/_sib23c03module.c
  # ${f2py_dir}/_sib23c03-f2pywrappers.f
  ${fortran_dir}/sibyll/sibyll2.3c03.f
  ${sibyll_init_source}
  ${logging_source}
  ${rangen_source}
  COMPILE_DEFS
  IMPY
)

### sib23d
f2py_add_module(_sib23d
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  # ${f2py_dir}/_sib23d.pyf
  # ${f2py_dir}/_sib23dmodule.c
  # ${f2py_dir}/_sib23d-f2pywrappers.f
  ${fortran_dir}/sibyll/sibyll2.3d.f
  ${sibyll_init_source}
  ${logging_source}
  ${rangen_source}
  COMPILE_DEFS
  IMPY
)

### qgs01
f2py_add_module(_qgs01
  FUNCTIONS
  cqgsini sectnu xxaini psconf xxreg psaini chepevt
  init_rmmard ${impy_logging}
  SOURCES
  # ${f2py_dir}/_qgs01.pyf
  # ${f2py_dir}/_qgs01module.c
  # ${f2py_dir}/_qgs01-f2pywrappers.f
  ${fortran_dir}/qgsjet/qgsjet01d.f
  ${fortran_dir}/qgsjet/impy_qgs1.f
  ${logging_source}
  ${rangen_source}
  COMPILE_DEFS
  IMPY
)

### qgsII03
set(QGSJETII_COMMON_FUNCTIONS
  cqgsini qgsect qgini qgconf qgreg chepevt init_rmmard ${impy_logging})

f2py_add_module(_qgsII03
  FUNCTIONS
  ${QGSJETII_COMMON_FUNCTIONS}
  SOURCES
  # ${f2py_dir}/_qgsII03.pyf
  # ${f2py_dir}/_qgsII03module.c
  # ${f2py_dir}/_qgsII03-f2pywrappers.f
  ${fortran_dir}/qgsjet/qgsjet-II-03.f
  ${fortran_dir}/qgsjet/impy_qgsII.f
  ${logging_source}
  ${rangen_source}
  COMPILE_DEFS
  IMPY
)

### qgsII04
f2py_add_module(_qgsII04
  FUNCTIONS
  ${QGSJETII_COMMON_FUNCTIONS}
  SOURCES
  # ${f2py_dir}/_qgsII04.pyf
  # ${f2py_dir}/_qgsII04module.c
  # ${f2py_dir}/_qgsII04-f2pywrappers.f
  ${fortran_dir}/qgsjet/qgsjet-II-04.f
  ${fortran_dir}/qgsjet/impy_qgsII.f
  ${logging_source}
  ${rangen_source}
  COMPILE_DEFS
  IMPY
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
  ${fortran_dir}/urqmd-3.4/impy_urqmd.f
  ${logging_source}
  ${rangen_source}
)
f2py_add_module(_urqmd34
  FUNCTIONS
  urqmd init uinit set0 params uounit strini loginit
  loadwtab norm_init output cascinit nucrad urqini partname
  chepevt ptsigtot init_rmmard ${impy_logging}
  SOURCES
  # ${f2py_dir}/_urqmd34.pyf
  # ${f2py_dir}/_urqmd34module.c
  # ${f2py_dir}/_urqmd34-f2pywrappers.f
  ${urqmd34_sources}
  COMPILE_DEFS
  IMPY
)
# to find terms.mod ?
target_include_directories(_urqmd34 PRIVATE ${fortran_dir}/urqmd-3.4/sources)

### pythia6
f2py_add_module(_pythia6
  FUNCTIONS
  pyinit pyexec pytune pylist pyevnt pyevnw pystat pyedit
  pyhepc pychge pycomp pyk init_rmmard ${impy_logging}
  SOURCES
  # ${f2py_dir}/_pythia6.pyf
  # ${f2py_dir}/_pythia6module.c
  # ${f2py_dir}/_pythia6-f2pywrappers.f
  # ${fortran_dir}/pythia6/pythia-6.4.28.f
  ${pythia6_pythia6428_source}
  ${logging_source}
  ${rangen_source}
  COMPILE_DEFS
  IMPY
)

### sophia
f2py_add_module(_sophia
  FUNCTIONS
  eventgen print_event crossection initial
  icon_pdg_sib init_rmmard toevt ${impy_logging}
  # INTERFACE_SOURCES
  # ${sophia_eventgen_source}
  # ${fortran_dir}/sophia/eventgen.f
  # ${rangen_source}
  # ${logging_source}
  SOURCES
  # ${f2py_dir}/_sophia.pyf
  # ${f2py_dir}/_sophiamodule.c
  # ${f2py_dir}/_sophia-f2pywrappers.f
  ${fortran_dir}/sophia/SOPHIA20.f
  ${fortran_dir}/sophia/eventgen.f
  ${fortran_dir}/sophia/sampling.f
  ${fortran_dir}/sophia/inpoutput.f
  # ${fortran_dir}/sophia/jetset74dp.f
  ${sophia_jetset74dp_source}
  ${fortran_dir}/sophia/impy_sophia.f
  ${logging_source}
  ${rangen_source}
  COMPILE_DEFS
  IMPY
)

### dpmjet306

set(dpmjet_dpmjet306f dpmjet3.0-6.f)

add_custom_command(
  OUTPUT ${dpmjet_dpmjet306f}
  COMMAND ${CMAKE_Fortran_COMPILER}
    -E -cpp ${fortran_dir}/dpmjet3.0-6/sources/dpmjet3.0-6.f
    -o ${dpmjet_dpmjet306f}
    -DIMPY
)


file(GLOB dpmjet306_sources ${fortran_dir}/dpmjet3.0-6/sources/*.f)
list(FILTER dpmjet306_sources EXCLUDE REGEX dpmjet3.0-6\.f)

f2py_add_module(_dpmjet306
  FUNCTIONS
  pho_event dt_init dt_kkinc idt_icihad dt_xsglau
  pycomp dt_initjs dt_rndmst dt_rndm dt_inucas idt_ipdgha
  pho_init pho_setpar pho_pname pho_pmass pho_setmdl
  pho_setpdf pycomp pho_xsect pho_borncs pho_harmci pho_fitout
  pho_mcini pho_ptcut pytune pho_rregpar pho_sregpar pho_prevnt
  ipho_pdg2id ipho_id2pdg pho_harint pho_harxto pho_harxpt
  pho_setpcomb dt_phoxs dt_xshn dt_flahad init_rmmard
  ${impy_logging}
  SOURCES
  # ${f2py_dir}/_dpmjet306.pyf
  # ${f2py_dir}/_dpmjet306module.c
  # ${f2py_dir}/_dpmjet306-f2pywrappers.f
  ${dpmjet306_sources}
  ${dpmjet_dpmjet306f}
  ${rangen_source}
  ${logging_source}
  COMPILE_DEFS
  IMPY
)



### dpmjetIII191
file(GLOB dpmjetIII191_sources
  ${fortran_dir}/dpmjetIII-19.1/src/phojet/*.f
  ${fortran_dir}/dpmjetIII-19.1/src/pythia/*.f
  ${fortran_dir}/dpmjetIII-19.1/src/dpmjet/*.f
  ${fortran_dir}/dpmjetIII-19.1/common/*.f
)

list(FILTER dpmjetIII191_sources EXCLUDE REGEX DT_RNDM\.f)
list(FILTER dpmjetIII191_sources EXCLUDE REGEX DT_RNDMST\.f)
list(FILTER dpmjetIII191_sources EXCLUDE REGEX DT_RNDMTE\.f)
list(FILTER dpmjetIII191_sources EXCLUDE REGEX PYR\.f)
# cmake does not like filenames with parantheses, but some Fortran files
# include those. As a workaround, we generate modified source files and
# renamed copies of the headers. Not elegant, but works. 🤷‍♂️
execute_process(
  COMMAND
  ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/dpmjetIII191_workaround.py
  ${CMAKE_CURRENT_BINARY_DIR} ${fortran_dir}/dpmjetIII-19.1/include
  ${dpmjetIII191_sources}
  OUTPUT_VARIABLE dpmjetIII191_modded_sources
)


f2py_add_module(_dpmjetIII191
  FUNCTIONS
  pho_event dt_init dt_kkinc
  idt_icihad dt_xsglau pycomp dt_initjs dt_rndmst dt_rndm dt_inucas idt_ipdgha dt_evtout
  pho_init pho_setpar poevt1 poevt2 pho_pname pho_pmass pho_setmdl
  pho_setpdf pycomp pho_xsect pho_borncs pho_harmci pho_fitout pho_mcini pho_ptcut
  pytune pho_rregpar pho_sregpar pho_prevnt ipho_pdg2id ipho_id2pdg pho_harint
  pho_harxto pho_harxpt pho_setpcomb
  dt_phoxs dt_xshn dt_flahad dt_title pho_ghhias init_rmmard
  ${impy_logging}
  SOURCES
  # ${f2py_dir}/_dpmjetIII191.pyf
  # ${f2py_dir}/_dpmjetIII191module.c
  # ${f2py_dir}/_dpmjetIII191-f2pywrappers.f
  # ${dpmjetIII191_modded_sources}
  ${dpmjetIII191_modded_sources}
  ${rangen_fpp_source}
  INCLUDE_DIRS
  ${fortran_dir}/dpmjetIII-19.1/include/phojet
  ${fortran_dir}/dpmjetIII-19.1/include/dpmjet
  ${fortran_dir}/dpmjetIII-19.1/include/pythia
  ${CMAKE_CURRENT_BINARY_DIR}/include
  COMPILE_DEFS
  IMPY
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
list(APPEND dpmjetIII193_sources ${rangen_fpp_source})


f2py_add_module(_dpmjetIII193
  FUNCTIONS
  pho_event dt_init dt_kkinc
  idt_icihad dt_xsglau pycomp dt_initjs dt_rndmst dt_rndm dt_inucas idt_ipdgha dt_evtout
  pho_init pho_setpar poevt1 poevt2 pho_pname pho_pmass pho_setmdl
  pho_setpdf pycomp pho_xsect pho_borncs pho_harmci pho_fitout pho_mcini pho_ptcut
  pytune pho_rregpar pho_sregpar pho_prevnt ipho_pdg2id ipho_id2pdg pho_harint
  pho_harxto pho_harxpt pho_setpcomb
  dt_phoxs dt_xshn dt_flahad dt_title pho_ghhias dt_getptn init_rmmard
  ${impy_logging}
  # INTERFACE_SOURCES
  # ${dpmjetIII193_interface_sources}
  SOURCES
  # ${f2py_dir}/_dpmjetIII193.pyf
  # ${f2py_dir}/_dpmjetIII193module.c
  # ${f2py_dir}/_dpmjetIII193-f2pywrappers.f
  ${dpmjetIII193_sources}
  INCLUDE_DIRS
  ${fortran_dir}/dpmjetIII-19.3/include/phojet
  ${fortran_dir}/dpmjetIII-19.3/include/dpmjet
  ${fortran_dir}/dpmjetIII-19.3/include/pythia
  ${fortran_dir}/dpmjetIII-19.3/include/flinclude
  COMPILE_DEFS
  IMPY
)