function (f2py_add_module target_name)

  cmake_parse_arguments(F2PY_ADD_MODULE
    ""
    ""
    "SIGNATURE;INTERFACE_SOURCES;SOURCES" 
    ${ARGN})

  if (NOT F2PY_ADD_MODULE_INTERFACE_SOURCES)
    set(F2PY_ADD_MODULE_INTERFACE_SOURCES ${F2PY_ADD_MODULE_SOURCES})
  endif()

  add_custom_command(
    OUTPUT
    ${target_name}module.c
    ${target_name}-f2pywrappers.f

    COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py
      -m ${target_name}
      -h ${target_name}.pyf
      --overwrite-signature only: ${F2PY_ADD_MODULE_SIGNATURE} :
      ${F2PY_ADD_MODULE_INTERFACE_SOURCES}
      > ${target_name}.log 2>&1

    COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py
      ${target_name}.pyf
      >> ${target_name}.log 2>&1

    DEPENDS ${F2PY_ADD_MODULE_SOURCES} ${F2PY_ADD_MODULE_INTERFACE_SOURCES}
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