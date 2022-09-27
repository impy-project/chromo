if(DEFINED ENV{IMPY_DEV_EXTRA_MODELS})
  ### sib23c00
  f2py_add_module(_sib23c00
    FUNCTIONS
    ${SIBYLL_COMMON_FUNCTIONS} sibhep3
    SOURCES
    ${fortran_dir}/sibyll/sibyll2.3c00.f
    ${fortran_dir}/sibyll/sibyll_init.fpp
    ${logging_source}
    ${rangen_source}
    COMPILE_DEFS
    ${impy_definitions}
    SIBYLL_TABLE_LENGTH=99
  )

  ### sib23c02
  f2py_add_module(_sib23c02
    FUNCTIONS
    ${SIBYLL_COMMON_FUNCTIONS} sibhep3
    SOURCES
    ${fortran_dir}/sibyll/sibyll2.3c02.f
    ${fortran_dir}/sibyll/sibyll_init.fpp
    ${logging_source}
    ${rangen_source}
    COMPILE_DEFS
    ${impy_definitions}
    SIBYLL_TABLE_LENGTH=99
  )

  ### sib23c03
  f2py_add_module(_sib23c03
    FUNCTIONS
    ${SIBYLL_COMMON_FUNCTIONS} sibhep3
    SOURCES
    ${fortran_dir}/sibyll/sibyll2.3c03.f
    ${fortran_dir}/sibyll/sibyll_init.fpp
    ${logging_source}
    ${rangen_source}
    COMPILE_DEFS
    ${impy_definitions}
    SIBYLL_TABLE_LENGTH=99
  )
  
endif()


if(DEFINED ENV{IMPY_DEV_DPMJETIII193})

  set (dev_source $ENV{IMPY_DEV_DPMJETIII193})

  message("IMPY_DEV_DPMJETIII193 = $ENV{IMPY_DEV_DPMJETIII193}")

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
    init_rmmard ${impy_logging}
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
    COMPILE_DEFS
    ${impy_definitions}
  )
  
endif()