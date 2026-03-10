INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}")
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_BINARY_DIR}")
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/urqmd34")
FILE(GLOB epos_files
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/addpart.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/coload.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/detbal.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/getmass.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/hepnam.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/iso.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/output.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/scatter.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/uaaa.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/angdis.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/boxprg.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/dwidth.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/getspin.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/init.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/ityp2pdg.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/numrec.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/paulibl.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/siglookup.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/upmerge.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/anndec.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/cascinit.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/dectim.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/error.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/hepchg.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/input.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/jdecay2.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/proppot.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/string.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/upyth.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/blockres.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/delpart.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/hepcmp.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/make22.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/saveinfo.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/tabinit.f
  ${CMAKE_CURRENT_SOURCE_DIR}/urqmd34/whichres.f
  ${CMAKE_CURRENT_SOURCE_DIR}/al.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/*.f
  ${CMAKE_CURRENT_BINARY_DIR}/*.f)



IF (CRMC_STATIC)
  add_library(CrmcBasic STATIC ${epos_files}) #always needed
ELSE (CRMC_STATIC)
  add_library(CrmcBasic SHARED ${epos_files}) #always needed
  add_library(Epos      SHARED ${base_files})
  INSTALL (TARGETS Epos
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib/static
    )
  target_link_libraries(Epos CrmcBasic)
ENDIF (CRMC_STATIC)
INSTALL (TARGETS CrmcBasic
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static)

