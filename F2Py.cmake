# Generate a python module from fortran code with numpy.f2py.
#
# f2py is only used to generate source code, all compilation is done
# with cmake, to make sure that a consistent single pair of C and Fortran
# compilers is used. numpy.f2py can also compile the code directly, but it
# selects the wrong compilers on some machines.
#
# f2py generates signature files (.pyf), which a user can edit. If you want
# to use a pre-generated .pyf file, pass it as part of the INTERFACE_SOURCES.
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
# f2py_add_module(target_name FUNCTIONS [...] INTERFACE_SOURCES [...] SOURCES [...])
#
# Keyword           : Value
# ------------------------------------------------------------------------------------
# FUNCTIONS         : List function names that should be exposed.
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
# A log file is generated as a side effect with name ${target_name}.log. The log is
# placed in the current build directory.

function (f2py_add_module target_name)

  cmake_parse_arguments(F2PY_ADD_MODULE
    ""
    ""
    "FUNCTIONS;INTERFACE_SOURCES;SOURCES" 
    ${ARGN})

  if (NOT F2PY_ADD_MODULE_INTERFACE_SOURCES)
    set(F2PY_ADD_MODULE_INTERFACE_SOURCES ${F2PY_ADD_MODULE_SOURCES})
  endif()

  set(F2PY_ADD_MODULE_PYF_FILE ${F2PY_ADD_MODULE_INTERFACE_SOURCES})
  list(FILTER F2PY_ADD_MODULE_PYF_FILE INCLUDE REGEX "\.pyf")
  list(FILTER F2PY_ADD_MODULE_INTERFACE_SOURCES EXCLUDE REGEX "\.pyf")
  
  # clear log file if it already exists
  file(WRITE ${target_name}.log "\n")

  if (NOT F2PY_ADD_MODULE_PYF_FILE)
    # generate new .pyf file
    set(F2PY_ADD_MODULE_PYF_FILE ${target_name}.pyf)

    add_custom_command(
      OUTPUT
      ${F2PY_ADD_MODULE_PYF_FILE}

      COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py
        -m ${target_name}
        -h ${F2PY_ADD_MODULE_PYF_FILE}
        --overwrite-signature only: ${F2PY_ADD_MODULE_FUNCTIONS} :
        ${F2PY_ADD_MODULE_INTERFACE_SOURCES}
        >> ${target_name}.log 2>&1

      DEPENDS ${F2PY_ADD_MODULE_INTERFACE_SOURCES}
    )
  endif()

  cmake_print_variables(F2PY_ADD_MODULE_PYF_FILE)
  cmake_print_variables(F2PY_ADD_MODULE_INTERFACE_SOURCES)

  add_custom_command(
    OUTPUT
    ${target_name}module.c
    ${target_name}-f2pywrappers.f

    COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py
      ${F2PY_ADD_MODULE_PYF_FILE}
      >> ${target_name}.log 2>&1

    DEPENDS ${F2PY_ADD_MODULE_SOURCES} ${F2PY_ADD_MODULE_PYF_FILE}
  )

  add_library(${target_name} MODULE
    ${F2PY_INCLUDE_DIR}/fortranobject.c
    ${target_name}module.c
    ${target_name}-f2pywrappers.f
    ${F2PY_ADD_MODULE_SOURCES}
  )
  target_link_libraries(${target_name} PRIVATE ${PYTHON_LIBRARIES})
  set_property(TARGET ${target_name} PROPERTY SUFFIX ${PYTHON_MODULE_EXTENSION})
  # must be a string, so that empty string works correcty
  set_property(TARGET ${target_name} PROPERTY PREFIX "${PYTHON_MODULE_PREFIX}")

endfunction()