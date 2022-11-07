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
#   f2py_source
#   f2py_dir
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
# COMPILE_DEFS      : List of definitions for compiler. If not empty the list is past
#                     to target_compile_definitions and to fortran_defs variable used
#                     for preprocessing. The definitions are intended for preprocessing
#                     of source files, i.e. it is implied that -D will prepended to each
#                     definition
#
# The module generates a target with target_name. You can manipulate this target like
# any other target in cmake to change its properties, for example, to set special
# compiler flags.
#
# A log file is generated as a side effect with name <target>.log. The log is placed
# in the current build directory.
#

function (f2py_add_module target_name)

  cmake_parse_arguments(F2PY_ADD_MODULE
    ""
    ""
    "FUNCTIONS;INCLUDE_DIRS;INTERFACE_SOURCES;SOURCES;COMPILE_DEFS"
    ${ARGN})


  # f2py files that exist in f2py source directory
  # In case of absence the files are regenerated
  set(pyf_file ${f2py_dir}/${target_name}.pyf)
  set(modulec_file ${f2py_dir}/${target_name}module.c)
  set(f2pywrap_file ${f2py_dir}/${target_name}-f2pywrappers.f)

  if (NOT F2PY_ADD_MODULE_INTERFACE_SOURCES)
    set(F2PY_ADD_MODULE_INTERFACE_SOURCES ${F2PY_ADD_MODULE_SOURCES})
  endif()

  set(log_file ${CMAKE_CURRENT_BINARY_DIR}/${target_name}.log)

  if (F2PY_ADD_MODULE_INCLUDE_DIRS)
    STRING(JOIN ":" _joined_dirs ${F2PY_ADD_MODULE_INCLUDE_DIRS})
    set(f2py_include_paths --include-paths ${_joined_dirs})
  endif()

  if (EXISTS ${pyf_file})
    file(WRITE ${log_file} "f2py_add_module: Use existing ${pyf_file}\n")
    message(STATUS "f2py_add_module: Use existing ${pyf_file}")
  else()
    file(WRITE ${log_file} "f2py_add_module: Generating ${pyf_file}\n")
    message(STATUS "f2py_add_module: Generating ${pyf_file}")
    # Definitions for source files processing
    set(fortran_defs)
    foreach(_def ${F2PY_ADD_MODULE_COMPILE_DEFS})
      STRING(APPEND fortran_defs "-D${_def} ")
    endforeach()
    
    # Source files processing for *.pyf
    set(processed_files)
    foreach(src_file ${F2PY_ADD_MODULE_INTERFACE_SOURCES})
      get_filename_component(src_filename ${src_file} NAME)
      set(proc_file CMakeFiles/${target_name}.dir/${src_filename})

      add_custom_command(
        OUTPUT ${proc_file}
        COMMAND ${CMAKE_Fortran_COMPILER}
        -E -cpp ${src_file} ${fortran_defs} -o ${proc_file}
        DEPENDS ${src_file}
      )
      list(APPEND processed_files ${proc_file})
    endforeach()

    get_filename_component(pyf_file_new ${pyf_file} NAME)
    
    # Generate in binary directory and copy to source directory.
    add_custom_command(
      OUTPUT ${pyf_file}
      COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py
        -m ${target_name}
        -h ${pyf_file_new}
        --overwrite-signature only: ${F2PY_ADD_MODULE_FUNCTIONS} :
        ${f2py_include_paths}
        ${processed_files}
        >> ${log_file} 2>&1
      COMMAND ${CMAKE_COMMAND} -E copy
        ${pyf_file_new} ${f2py_dir}
      DEPENDS ${processed_files}
    )

  endif()
  
  if ((EXISTS ${modulec_file}) AND (EXISTS ${f2pywrap_file}))
    file(APPEND ${log_file} "f2py_add_module: Use existing ${modulec_file} ${f2pywrap_file}\n")
    message(STATUS "f2py_add_module: Use existing ${modulec_file}")
    message(STATUS "f2py_add_module: Use existing ${f2pywrap_file}")
  else()
    message(STATUS "f2py_add_module: Generating ${modulec_file}")
    message(STATUS "f2py_add_module: Generating ${f2pywrap_file}")

    get_filename_component(modulec_file_new ${modulec_file} NAME)
    get_filename_component(f2pywrap_file_new ${f2pywrap_file} NAME)

    # Generate in binary directory and copy to source directory.
    add_custom_command(
      OUTPUT ${modulec_file} ${f2pywrap_file}
      COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py
        ${pyf_file}
        ${f2py_include_paths}
        >> ${log_file} 2>&1 
      COMMAND ${CMAKE_COMMAND} -E copy
        ${modulec_file_new} ${f2pywrap_file_new} ${f2py_dir}
      DEPENDS ${F2PY_ADD_MODULE_SOURCES} ${pyf_file}
    )

  endif()

  add_library(${target_name} MODULE
    ${f2py_source}
    ${modulec_file}
    ${f2pywrap_file}
    ${F2PY_ADD_MODULE_SOURCES}
  )

  if (PYTHON_LIBRARIES) # may not be available (e.g. on manylinux)
    target_link_libraries(${target_name} PRIVATE ${PYTHON_LIBRARIES})
  endif()
  if (F2PY_ADD_MODULE_INCLUDE_DIRS)
    target_include_directories(${target_name}
    PRIVATE ${F2PY_ADD_MODULE_INCLUDE_DIRS})
  endif()

  if (F2PY_ADD_MODULE_COMPILE_DEFS)
    target_compile_definitions(${target_name}
    PRIVATE ${F2PY_ADD_MODULE_COMPILE_DEFS})
  endif()  
  target_compile_options(${target_name} PRIVATE -cpp)
  # Link dll statically in Windows
  # It can be potentially a PROBLEM!!!
  # But this is the only way found to build a working library on Windows.
  # However, on Linux and MacOS, it throws an error at the linking stage.
  if (WIN32)
    target_link_libraries(${target_name} PUBLIC "-static")
  endif()
  set_property(TARGET ${target_name} PROPERTY SUFFIX ${PYTHON_MODULE_EXTENSION})
  # must be a string, so that empty string works correcty
  set_property(TARGET ${target_name} PROPERTY PREFIX "${PYTHON_MODULE_PREFIX}")

endfunction()