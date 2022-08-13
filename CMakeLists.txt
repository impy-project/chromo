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
set(fortran_src_dir ${CMAKE_SOURCE_DIR}/src)

set(impy_logging impy_openlogfile impy_closelogfile)
set(logging_source ${fortran_src_dir}/logging.f)
set(rangen_source rangen.f)
set(rangen_sib21_source rangen_sib21.f)

# WARNING this is not platform-independent and needs to be adjusted
# when new fortran compilers should be supported
add_custom_command(
  OUTPUT ${rangen_source}
  COMMAND ${CMAKE_Fortran_COMPILER}
    -E ${fortran_src_dir}/rangen.fpp
    -o ${rangen_source}
)
add_custom_command(
  OUTPUT ${rangen_sib21_source}
  COMMAND ${CMAKE_Fortran_COMPILER}
    -E ${fortran_src_dir}/rangen.fpp
    -o ${rangen_sib21_source}
    -DSIBYLL_SP
)

### eposlhc
file(GLOB eposlhc_sources
  ${fortran_src_dir}/EPOS-LHC/sources/*.f
  ${fortran_src_dir}/EPOS-LHC/eposlhc.pyf
)
list(FILTER eposlhc_sources EXCLUDE REGEX epos_example\.f)

f2py_add_module(_eposlhc
  FUNCTIONS
  aaset ainit aepos afinal hepmcstore
  getcharge idtrafo initializeepos initeposevt
  setstable setunstable xsection ${impy_logging}
  SOURCES
  ${fortran_src_dir}/EPOS-LHC/eposlhc.pyf
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
  ${fortran_src_dir}/sibyll/sib21.pyf
  ${fortran_src_dir}/sibyll/sibyll_21.f
  ${fortran_src_dir}/sibyll/sib21aux.f
  ${fortran_src_dir}/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_sib21_source}
)

### sib23
f2py_add_module(_sib23
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${fortran_src_dir}/sibyll/sib23.pyf
  ${fortran_src_dir}/sibyll/sibyll2.3.f
  ${fortran_src_dir}/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_source}
)

### sib23c00
f2py_add_module(_sib23c00
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${fortran_src_dir}/sibyll/sib23c00.pyf
  ${fortran_src_dir}/sibyll/sibyll2.3c00.f
  ${fortran_src_dir}/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_source}
)

### sib23c01
f2py_add_module(_sib23c01
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${fortran_src_dir}/sibyll/sib23c01.pyf
  ${fortran_src_dir}/sibyll/sibyll2.3c01.f
  ${fortran_src_dir}/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_source}
)

### sib23c02
f2py_add_module(_sib23c02
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${fortran_src_dir}/sibyll/sib23c02.pyf
  ${fortran_src_dir}/sibyll/sibyll2.3c02.f
  ${fortran_src_dir}/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_source}
)

### sib23c03
f2py_add_module(_sib23c03
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${fortran_src_dir}/sibyll/sib23c03.pyf
  ${fortran_src_dir}/sibyll/sibyll2.3c03.f
  ${fortran_src_dir}/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_source}
)

### sib23d
f2py_add_module(_sib23d
  FUNCTIONS
  ${SIBYLL_COMMON_FUNCTIONS} sibhep3
  SOURCES
  ${fortran_src_dir}/sibyll/sib23d.pyf
  ${fortran_src_dir}/sibyll/sibyll2.3d.f
  ${fortran_src_dir}/sibyll/sibyll_init.f
  ${logging_source}
  ${rangen_source}
)

### qgs01
f2py_add_module(_qgs01
  FUNCTIONS
  cqgsini sectnu xxaini psconf xxreg psaini chepevt
  ${impy_logging}
  SOURCES
  ${fortran_src_dir}/qgsjet/qgs01.pyf
  ${fortran_src_dir}/qgsjet/qgsjet01d.f
  ${fortran_src_dir}/qgsjet/impy_qgs1.f
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
  ${fortran_src_dir}/qgsjet/qgsII-03.pyf
  ${fortran_src_dir}/qgsjet/qgsjet-II-03.f
  ${fortran_src_dir}/qgsjet/impy_qgsII.f
  ${logging_source}
  ${rangen_source}
)

### qgsII04
f2py_add_module(_qgsII04
  FUNCTIONS
  ${QGSJETII_COMMON_FUNCTIONS}
  SOURCES
  ${fortran_src_dir}/qgsjet/qgsII-04.pyf
  ${fortran_src_dir}/qgsjet/qgsjet-II-04.f
  ${fortran_src_dir}/qgsjet/impy_qgsII.f
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
list(TRANSFORM urqmd34_sources PREPEND ${fortran_src_dir}/urqmd-3.4/sources/)
list(APPEND urqmd34_sources
  ${fortran_src_dir}/urqmd-3.4/impy_urqmd.f
  ${logging_source}
  ${rangen_source}
)
f2py_add_module(_urqmd34
  FUNCTIONS
  urqmd init uinit set0 params uounit strini loginit
  loadwtab norm_init output cascinit nucrad urqini partname
  chepevt ptsigtot init_rmmard ${impy_logging}
  SOURCES
  # Using pre-generated pyf file results in this error
  #  Undefined symbols for architecture x86_64:
  #  "_f2pyinitterms_", referenced from:
  #  _f2py_init_terms in _urqmd34module.c.o
  ${urqmd34_sources}
)

### pythia6
f2py_add_module(_pythia6
  FUNCTIONS
  pyinit pyexec pytune pylist pyevnt pyevnw pystat pyedit
  pyhepc pychge pycomp pyk ${impy_logging}
  SOURCES
  ${fortran_src_dir}/pythia6/pythia6.pyf
  ${fortran_src_dir}/pythia6/pythia-6.4.28.f
  ${logging_source}
)

### sophia
f2py_add_module(_sophia
  FUNCTIONS
  eventgen print_event crossection initial
  icon_pdg_sib init_rmmard toevt ${impy_logging}
  SOURCES
  ${fortran_src_dir}/sophia/sophia.pyf
  ${fortran_src_dir}/sophia/SOPHIA20.f
  ${fortran_src_dir}/sophia/eventgen.f
  ${fortran_src_dir}/sophia/sampling.f
  ${fortran_src_dir}/sophia/inpoutput.f
  ${fortran_src_dir}/sophia/jetset74dp.f
  ${fortran_src_dir}/sophia/impy_sophia.f
  ${logging_source}
  ${rangen_source}
)

### dpmjet306
file(GLOB dpmjet306_sources ${fortran_src_dir}/dpmjet3.0-6/sources/*.f)
f2py_add_module(_dpmjet306
  FUNCTIONS
  pho_event dt_init dt_kkinc idt_icihad dt_xsglau
  pycomp dt_initjs dt_rndmst dt_rndm dt_inucas idt_ipdgha
  pho_init pho_setpar pho_pname pho_pmass pho_setmdl
  pho_setpdf pycomp pho_xsect pho_borncs pho_harmci pho_fitout
  pho_mcini pho_ptcut pytune pho_rregpar pho_sregpar pho_prevnt
  ipho_pdg2id ipho_id2pdg pho_harint pho_harxto pho_harxpt
  pho_setpcomb dt_phoxs dt_xshn dt_flahad
  ${impy_logging}
  SOURCES
  ${fortran_src_dir}/dpmjet3.0-6/dpmjet306.pyf
  ${dpmjet306_sources}
  ${logging_source}
)

### dpmjetIII191 does not compile in CMake
if (FALSE)
  # make[2]: *** No rule to make target `../src/dpmjetIII-19.1/include/flinclude/(COMPUT)', needed by `CMakeFiles/libdpmjetIII191.dir/src/dpmjetIII-19.1/common/oauxfi.f.o'.  Stop.
  # This is silly, because ../src/dpmjetIII-19.1/include/flinclude/(COMPUT) exists
  # and does not need to be created
  file(GLOB dpmjetIII191_sources
    ${fortran_src_dir}/dpmjetIII-19.1/src/phojet/*.f
    ${fortran_src_dir}/dpmjetIII-19.1/src/pythia/*.f
    ${fortran_src_dir}/dpmjetIII-19.1/src/dpmjet/*.f
    ${fortran_src_dir}/dpmjetIII-19.1/common/*.f
  )
  # sources need preprocessing
  set_source_files_properties(
    ${dpmjetIII191_sources}
    PROPERTIES
    Fortran_PREPROCESS True
  )
  add_library(libdpmjetIII191 STATIC
    ${dpmjetIII191_sources}
  )
  target_include_directories(libdpmjetIII191
    PUBLIC
    ${fortran_src_dir}/dpmjetIII-19.1/include/phojet
    ${fortran_src_dir}/dpmjetIII-19.1/include/dpmjet
    ${fortran_src_dir}/dpmjetIII-19.1/include/pythia
    ${fortran_src_dir}/dpmjetIII-19.1/include/flinclude
  )
  f2py_add_module(_dpmjetIII191
    FUNCTIONS
    pho_event dt_init dt_kkinc
    idt_icihad dt_xsglau pycomp dt_initjs dt_rndmst dt_rndm dt_inucas idt_ipdgha dt_evtout
    pho_init pho_setpar poevt1 poevt2 pho_pname pho_pmass pho_setmdl
    pho_setpdf pycomp pho_xsect pho_borncs pho_harmci pho_fitout pho_mcini pho_ptcut
    pytune pho_rregpar pho_sregpar pho_prevnt ipho_pdg2id ipho_id2pdg pho_harint
    pho_harxto pho_harxpt pho_setpcomb
    dt_phoxs dt_xshn dt_flahad dt_title pho_ghhias
    ${impy_logging}
    SOURCES
    ${fortran_src_dir}/dpmjetIII-19.1/common/dpmjetIII191.pyf
  )
  target_link_libraries(_dpmjetIII191 PRIVATE libdpmjetIII191)
endif()