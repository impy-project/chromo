# Generate a python module from fortran code (optionally) with numpy.f2py.
#
# f2py is only used to generate source code. All compilation is done
# with cmake, to make sure that a consistent single pair of C and Fortran
# compilers is used. numpy.f2py can also compile the code directly, but it
# selects the wrong compilers on some machines.
#
# f2py generates an interface file (.pyf), which a user can edit. If you want
# to use a pre-generated .pyf file, pass it as part of the SOURCES.
#
# f2py further generates two files, <target>module.c and 
# <target>-f2pywrappers.f, which only change if the interface file changes.
# If you want to use pre-generated files, pass them as part of the SOURCES.
#
# This function expects the following variables to be set:
#
#   PYTHON_EXECUTABLE
#   PYTHON_LIBRARIES
#   F2PY_INCLUDE_DIR
#   PYTHON_MODULE_EXTENSION
#   PYTHON_MODULE_PREFIX
#
# The following arguments are accepted:
#
# f2py_add_module(target_name FUNCTIONS [...] INCLUDE_DIRS [...] 
#                 INTERFACE_SOURCES [...] SOURCES [...])
#
# Keyword           : Value
# ------------------------------------------------------------------------------------
# FUNCTIONS         : List function names that should be exposed.
# INCLUDE_DIRS      : Any include directories that the fortran code needs. We cannot
#                     simply use target_include_directories because f2py also needs
#                     these files, too, if interface files are to be generated.
# INTERFACE_SOURCES : List of sources that f2py needs to generate the .pyf file. This
#                     is a subset of SOURCES. This is optional. If it is not supplied,
#                     SOURCES is used.
# SOURCES           : List of sources that are needed to compile the module.
#                     These can be C files, Fortran files, object files. A single .pyf
#                     file is also accepted.
#
# The module generates a target with target_name. You can manipulate this target like
# any other target in cmake to change its properties, for example, to set special
# compiler flags.
#
# A log file is generated as a side effect with name <target>.log. The log is placed
# in the current build directory.

function (f2py_add_module target_name)

  cmake_parse_arguments(F2PY_ADD_MODULE
    ""
    ""
    "FUNCTIONS;INCLUDE_DIRS;INTERFACE_SOURCES;SOURCES;COMPILE_DEFS"
    ${ARGN})

  if (NOT F2PY_ADD_MODULE_INTERFACE_SOURCES)
    set(F2PY_ADD_MODULE_INTERFACE_SOURCES ${F2PY_ADD_MODULE_SOURCES})
  endif()

  # filter out pyf file if included in sources
  set(F2PY_ADD_MODULE_PYF_FILE ${F2PY_ADD_MODULE_INTERFACE_SOURCES})
  list(FILTER F2PY_ADD_MODULE_PYF_FILE INCLUDE REGEX "\\.pyf$")
  list(FILTER F2PY_ADD_MODULE_INTERFACE_SOURCES EXCLUDE REGEX "\\.pyf$")

  # filter out *module.c and -f2pywrappers.f files if included in sources
  set(F2PY_ADD_MODULE_GEN_1 ${F2PY_ADD_MODULE_INTERFACE_SOURCES})
  list(FILTER F2PY_ADD_MODULE_GEN_1 INCLUDE REGEX "${target_name}module\\.c$")
  set(F2PY_ADD_MODULE_GEN_2 ${F2PY_ADD_MODULE_INTERFACE_SOURCES})
  list(FILTER F2PY_ADD_MODULE_GEN_2 INCLUDE REGEX "${target_name}-f2pywrappers.*")

  set(F2PY_ADD_MODULE_LOG_FILE ${CMAKE_CURRENT_BINARY_DIR}/${target_name}.log)

  set(F2PY_ADD_MODULE_INC)
  if (F2PY_ADD_MODULE_INCLUDE_DIRS)
    STRING(JOIN ":" _joined_dirs ${F2PY_ADD_MODULE_INCLUDE_DIRS})
    set(F2PY_ADD_MODULE_INC --include-paths ${_joined_dirs})
  endif()
  
  if (NOT F2PY_ADD_MODULE_PYF_FILE)
    # Set directory name to output generated file *.pyf files
    # *module.c and *-f2pywrappers.f
    set(model_out ${CMAKE_CURRENT_BINARY_DIR}/${target_name}_out)
    # Set directory name to output processed source files
    # needed to *.pyf file generation
    set(pyf_sources "CMakeFiles/${target_name}.dir/pyf_sources")
    
    set(F2PY_ADD_MODULE_PYF_FILE ${model_out}/${target_name}.pyf)
    file(WRITE ${F2PY_ADD_MODULE_LOG_FILE} "f2py_add_module: Generating ${F2PY_ADD_MODULE_PYF_FILE}\n")

    # Definitions for source files processing
    set(fortran_defs)
    foreach(_def ${F2PY_ADD_MODULE_COMPILE_DEFS})
      STRING(APPEND fortran_defs "-D${_def} ")
    endforeach()

    # Target to make directories
    set(mkdir_pyf_sources "mkdir${target_name}_pyf")
    add_custom_target(${mkdir_pyf_sources}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${pyf_sources}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${model_out})
    
    # Source files processing for *.pyf
    set(processed_files)
    foreach(_file ${F2PY_ADD_MODULE_INTERFACE_SOURCES})
      string(REGEX REPLACE "\.f$" "\.f\.proc" proc_file ${_file})
      string(REGEX REPLACE "\.fpp$" "\.f\.proc"  proc_file ${proc_file})
      get_filename_component(barename ${proc_file} NAME)
      set(proc_file ${pyf_sources}/${barename})

      add_custom_command(
        OUTPUT ${proc_file}
        COMMAND ${CMAKE_Fortran_COMPILER}
        -E -cpp ${_file} ${fortran_defs} -o ${proc_file}
        DEPENDS ${mkdir_pyf_sources}
      )
      list(APPEND processed_files ${proc_file})
    endforeach()

    # *.pyf file generation
    add_custom_command(
      OUTPUT
      ${F2PY_ADD_MODULE_PYF_FILE}

      COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py
        -m ${target_name}
        -h ${F2PY_ADD_MODULE_PYF_FILE}
        --overwrite-signature only: ${F2PY_ADD_MODULE_FUNCTIONS} :
        ${F2PY_ADD_MODULE_INC}
        ${processed_files}
        >> ${F2PY_ADD_MODULE_LOG_FILE} 2>&1

      DEPENDS ${processed_files}
    )
  else()
    file(WRITE ${F2PY_ADD_MODULE_LOG_FILE} "f2py_add_module: Use existing ${F2PY_ADD_MODULE_PYF_FILE}\n")
    message(STATUS "f2py_add_module: Use existing ${F2PY_ADD_MODULE_PYF_FILE}")
  endif()


  if (F2PY_ADD_MODULE_GEN_1 AND F2PY_ADD_MODULE_GEN_2)
    file(APPEND ${F2PY_ADD_MODULE_LOG_FILE} "f2py_add_module: Use existing ${F2PY_ADD_MODULE_GEN_1} ${F2PY_ADD_MODULE_GEN_2}\n")
    message(STATUS "f2py_add_module: Use existing ${F2PY_ADD_MODULE_GEN_1}")
    message(STATUS "f2py_add_module: Use existing ${F2PY_ADD_MODULE_GEN_2}")
  else()
    set(F2PY_ADD_MODULE_GEN_1 ${model_out}/${target_name}module.c)
    set(F2PY_ADD_MODULE_GEN_2 ${model_out}/${target_name}-f2pywrappers.f)

    add_custom_command(
      OUTPUT ${F2PY_ADD_MODULE_GEN_1} ${F2PY_ADD_MODULE_GEN_2}

      COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py
        ${F2PY_ADD_MODULE_PYF_FILE}
        ${F2PY_ADD_MODULE_INC}
        >> ${F2PY_ADD_MODULE_LOG_FILE} 2>&1
      WORKING_DIRECTORY ${model_out}
      DEPENDS ${F2PY_ADD_MODULE_SOURCES} ${F2PY_ADD_MODULE_PYF_FILE}
    )
  endif()

  add_library(${target_name} MODULE
    ${f2py_source}
    ${F2PY_ADD_MODULE_GEN_1}
    ${F2PY_ADD_MODULE_GEN_2}
    ${F2PY_ADD_MODULE_SOURCES}
  )
  if (PYTHON_LIBRARIES) # may not be available (e.g. on manylinux)
    target_link_libraries(${target_name} PRIVATE ${PYTHON_LIBRARIES})
  endif()
  if (F2PY_ADD_MODULE_INCLUDE_DIRS)
    target_include_directories(${target_name} PRIVATE ${F2PY_ADD_MODULE_INCLUDE_DIRS})
  endif()

  if (F2PY_ADD_MODULE_COMPILE_DEFS)
    target_compile_definitions(${target_name} PRIVATE ${F2PY_ADD_MODULE_COMPILE_DEFS})
  endif()  
  target_compile_options(${target_name} PRIVATE -cpp)
  set_property(TARGET ${target_name} PROPERTY SUFFIX ${PYTHON_MODULE_EXTENSION})
  # must be a string, so that empty string works correcty
  set_property(TARGET ${target_name} PROPERTY PREFIX "${PYTHON_MODULE_PREFIX}")


endfunction()