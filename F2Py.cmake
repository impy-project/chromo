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
    "FUNCTIONS;INCLUDE_DIRS;INTERFACE_SOURCES;SOURCES"
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
    set(F2PY_ADD_MODULE_PYF_FILE ${CMAKE_CURRENT_BINARY_DIR}/${target_name}.pyf)
    file(WRITE ${F2PY_ADD_MODULE_LOG_FILE} "f2py_add_module: Generating ${F2PY_ADD_MODULE_PYF_FILE}\n")

    add_custom_command(
      OUTPUT
      ${F2PY_ADD_MODULE_PYF_FILE}

      COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py
        -m ${target_name}
        -h ${F2PY_ADD_MODULE_PYF_FILE}
        --overwrite-signature only: ${F2PY_ADD_MODULE_FUNCTIONS} :
        ${F2PY_ADD_MODULE_INC}
        ${F2PY_ADD_MODULE_INTERFACE_SOURCES}
        >> ${F2PY_ADD_MODULE_LOG_FILE} 2>&1

      DEPENDS ${F2PY_ADD_MODULE_INTERFACE_SOURCES}
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
    set(F2PY_ADD_MODULE_GEN_1 ${CMAKE_CURRENT_BINARY_DIR}/${target_name}module.c)
    set(F2PY_ADD_MODULE_GEN_2 ${CMAKE_CURRENT_BINARY_DIR}/${target_name}-f2pywrappers.f)

    add_custom_command(
      OUTPUT ${F2PY_ADD_MODULE_GEN_1} ${F2PY_ADD_MODULE_GEN_2}

      COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py
        ${F2PY_ADD_MODULE_PYF_FILE}
        ${F2PY_ADD_MODULE_INC}
        >> ${F2PY_ADD_MODULE_LOG_FILE} 2>&1

      DEPENDS ${F2PY_ADD_MODULE_SOURCES} ${F2PY_ADD_MODULE_PYF_FILE}
    )
  endif()

  add_library(${target_name} MODULE
    ${F2PY_INCLUDE_DIR}/fortranobject.c
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
  set_property(TARGET ${target_name} PROPERTY SUFFIX ${PYTHON_MODULE_EXTENSION})
  # must be a string, so that empty string works correcty
  set_property(TARGET ${target_name} PROPERTY PREFIX "${PYTHON_MODULE_PREFIX}")

endfunction()