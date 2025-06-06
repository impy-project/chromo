
project('chromo', 'c', 'cpp', 'fortran',
  license: 'BSD-3-Clause',
  version: run_command([
    'python3', meson.project_source_root() + '/scripts/get_version.py'
    ], check: true).stdout().strip(),
  meson_version: '>=1.4.0',
  default_options: [
    'buildtype=release',
    'cpp_std=c++14',
    'fortran_std=legacy',
    'warning_level=1',
    'optimization=3'
  ]
)
add_project_arguments('-Wno-uninitialized', '-fPIC',
                      '-fno-second-underscore', '-Wno-unused', '-Wno-conversion',
                      '-Wno-unused-dummy-argument', '-Wno-character-truncation', 
                      '-Wno-integer-division','-Wno-surprising', '-Wno-tabs',
                      '-cpp', '-UNDEBUG',
                      language: 'fortran')
add_project_arguments('-fPIC',
                      '-Wno-unused-const-variable',
                      language: 'C')

# Directory for utilities
scripts_dir = meson.project_source_root() + '/scripts/'

# Update submodules if needed
run_command('git', 'submodule', 'update', '--init', '--recursive',
  check: true, capture: true)

# Find Python and its development headers
py = import('python').find_installation()
py_dep = py.dependency()

# Find numpy and its headers
numpy_dep = dependency('numpy')
numpy_inc = run_command(py, '-c', 'import numpy; print(numpy.get_include())', check: true).stdout().strip()

# Find f2py include directory
f2py_inc = run_command(py, '-c', 
  'import numpy.f2py; from pathlib import Path; print(Path(numpy.f2py.__file__).parent / "src")',
  check: true).stdout().strip()

# F2py source files
f2py_fortranobject_c = f2py_inc / 'fortranobject.c'

# Add include directories
chromo_inc_dirs = [
  numpy_inc,
  f2py_inc,
  'src/fortran',
]

# Add compile definitions
add_project_arguments('-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', language: ['c', 'cpp'])

# Common source files
fortran_dir = 'src/fortran'
logging_source = ['src/fortran/logging.f']
rangen_source = ['src/fortran/rangen.c', 'src/fortran/rangen.fpp']
normal_source = ['src/fortran/normal.c']
chromo_common_functions = 'chromo_openlogfile,chromo_closelogfile,npyrng'
chromo_common_flags = ['-DCHROMO']

# SIBYLL-specific definitions
sibyll_sources_base = [
  'src/fortran/sibyll/sibyll_init.fpp',
  logging_source,
  rangen_source,
]
sibyll_common_functions = [
    'sibyll',
    'sibyll_ini',
    'sib_sigma_hp',
    'sib_sigma_hair',
    'sib_sigma_hnuc',
    'int_nuc',
    'decsib',
    'decpar',
    'sibini',
    'sibhep',
    'sib_list',
    'isib_pid2pdg',
    'isib_pdg2pid',
    'pdg_ini'
] + chromo_common_functions

# QGSJET-specific definitions

qgsjetII_common_functions = [
    'cqgsini',
    'qgsect',
    'qgini',
    'qgconf',
    'qgreg',
    'chepevt',
    'qgcrossc',
    'cqgshh_ha_cs'
] + chromo_common_functions

# DPMJET and PHOJET-specific definitions


# Paramterize all targets. The models must have an _ underline prefix
# to avoid conflicts with the Python module names.
models = {
  '_sib21': {
    'sources': sibyll_sources_base + ['src/fortran/sibyll/sibyll_21.f', 'src/fortran/sibyll/sib21aux.f'],
    'functions': ','.join(sibyll_common_functions + ['spgasdev', 'glauber']),
    'flags' : ['-DSIBYLL_21'],
  },
  '_sib23': {
    'sources': sibyll_sources_base + ['src/fortran/sibyll/sibyll2.3.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_sib23c00': {
    'sources': sibyll_sources_base + ['src/fortran/sibyll/sibyll2.3c00.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_sib23c01': {
    'sources': sibyll_sources_base + ['src/fortran/sibyll/sibyll2.3c01.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_sib23c02': {
    'sources': sibyll_sources_base + ['src/fortran/sibyll/sibyll2.3c02.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_sib23c03': {
    'sources': sibyll_sources_base + ['src/fortran/sibyll/sibyll2.3c03.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_sib23d': {
    'sources': sibyll_sources_base + ['src/fortran/sibyll/sibyll2.3d.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_sib23d_star': {
    'sources': sibyll_sources_base + ['src/fortran/sibyll/sibyll2.3d-star-p02.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_qgs01': {
    'sources': [
      'src/fortran/qgsjet/qgsjet01d.f',
      'src/fortran/qgsjet/chromo_qgs1.f',
      logging_source,
      rangen_source
    ],
    'functions': ','.join([
      'cqgsini','sectnu','xxaini','psconf','xxreg','psaini','chepevt','xxfz',
      'cqgshh_ha_cs','crossc'] + chromo_common_functions),
  },
  '_qgsII03': {
    'sources': [
      'src/fortran/qgsjet/qgsjet-II-03.f',
      'src/fortran/qgsjet/chromo_qgsII.f',
      logging_source,
      rangen_source
    ],
    'functions': ','.join(qgsjetII_common_functions)
  },
  '_qgsII04': {
    'sources': [
      'src/fortran/qgsjet/qgsjet-II-04.f',
      'src/fortran/qgsjet/chromo_qgsII.f',
      logging_source,
      rangen_source
    ],
    'functions': ','.join(qgsjetII_common_functions)
  },
  '_qgsIII': {
    'sources': [
      'src/fortran/qgsjet/qgsjet-III.f',
      'src/fortran/qgsjet/chromo_qgsIII.f',
      logging_source,
      rangen_source
    ],
    'functions': ','.join(qgsjetII_common_functions + ['qgran'])
  },
}

# Read from pyproject.toml which models to build
enabled_models = run_command([
  py, scripts_dir + 'get_models.py', '--enabled'], 
  check: true).stdout().strip().split('\n')
disabled_models = run_command([
  py, scripts_dir + 'get_models.py', '--disabled'], 
  check: true).stdout().strip().split('\n')

foreach f2py_name, model_info : models
  if f2py_name.substring(1) not in enabled_models
    message('Skipping model: ' + f2py_name + ' (disabled)')
    continue
  endif

  message('Building model: ' + f2py_name)
  output = [f2py_name + 'module.c', f2py_name + '-f2pywrappers.f', f2py_name + '.pyf']
  includes = chromo_inc_dirs + model_info.get('include_directories', [])
  include_args = []
  foreach p : includes
    include_args += '-I' + (meson.project_source_root() / p)
  endforeach

  foreach f : files(model_info.get('interface_sources', model_info['sources']))
    message(f.full_path())
  endforeach
  flags = chromo_common_flags + model_info.get('flags', [])
  message('Flags: ' + ' '.join(flags))

  interface_sources = files(model_info.get('interface_sources', model_info['sources']))
  sources = files(model_info['sources'])

  # Custom target to generate f2py wrapper
  f2py_target = custom_target(f2py_name,
    output: output,
    input: interface_sources,
    command: [
      py, scripts_dir + 'generate_f2py.py',
      f2py_name,
      model_info['functions'],
      meson.global_build_root() / 'meson-logs',
      ' '.join(include_args),
      ' '.join(flags),
      '@OUTDIR@',
      '@INPUT@'
    ],
    build_by_default: true
  )

  foreach f : output
    message('F2PY output: ' + f)
  endforeach
  # Create the sib21 extension module
  py.extension_module(f2py_name,
    sources: sources + f2py_target + [f2py_fortranobject_c] + [normal_source],
    include_directories: include_directories(includes),
    dependencies: [py_dep, numpy_dep],
    c_args: ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION'] + flags,
    fortran_args: flags,
    install: true,
    subdir: 'chromo/models',
  )
endforeach




# EPOS-LHC sources
# epos_sources = [
#   files('src/fortran/epos/sources/epos-bas-lhc.f'),
#   files('src/fortran/epos/sources/epos-hnb-lhc.f'),
#   files('src/fortran/epos/sources/epos-dro-lhc.f'),
#   files('src/fortran/epos/sources/epos-con-lhc.f'),
#   files('src/fortran/epos/sources/epos-dky-lhc.f'),
#   files('src/fortran/epos/sources/epos-int-lhc.f'),
#   files('src/fortran/epos/sources/epos-jps-lhc.f'),
#   files('src/fortran/epos/sources/epos-omg-lhc.f'),
#   files('src/fortran/epos/sources/epos-qsh-lhc.f'),
#   files('src/fortran/epos/sources/epos-rsh-lhc.f'),
#   files('src/fortran/epos/sources/epos-sem-lhc.f'),
#   files('src/fortran/epos/sources/epos-uti-lhc.f'),
#   files('src/fortran/epos/sources/epos-xan-lhc.f'),
#   files('src/fortran/epos/sources/epos-ids-lhc.f'),
#   files('src/fortran/epos/sources/epos_interface.f'),
#   files('src/fortran/epos/epos-random-dummy.f'),
#   logging_source,
#   rangen_source,
# ]

# Custom target to generate f2py wrapper for EPOS
# epos_f2py_target = custom_target('epos_f2py',
#   output: ['_eposlhcmodule.c', '_eposlhc-f2pywrappers.f', '_eposlhc.pyf'],
#   input: epos_sources[0],  # Use first source as input for dependency tracking
#   command: [
#     py, files('generate_f2py.py'),
#     '_eposlhc',
#     'aaset,ainit,aepos,alist,afinal,hepmcstore,getcharge,idtrafo,initepos,initeposevt,xsection,rmmard,ranfgt,ranfst',
#     '@INPUT@',
#     '',
#     'src/fortran',
#     '@OUTDIR@'
#   ],
#   build_by_default: true
# )

# # Create the EPOS-LHC extension module
# py.extension_module('_eposlhc',
#   sources: epos_sources + epos_f2py_target + [f2py_fortranobject_c] + [normal_source],
#   include_directories: [inc_dirs],
#   dependencies: [py_dep, numpy_dep],
#   c_args: ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', ],
#   fortran_args: ['-DCHROMO'],
#   install: true,
#   subdir: 'chromo/models',
# )

# DPMJET III-19.1 sources - DPMJET-specific include directories
# dpmjet_inc_dirs = [
#   numpy_inc,
#   f2py_inc,
#   'src/fortran',
#   'src/fortran/dpmjetIII-19.1/include/dpmjet',
#   'src/fortran/dpmjetIII-19.1/include/phojet',
#   'src/fortran/dpmjetIII-19.1/include/flinclude',
#   'src/fortran/dpmjetIII-19.1/include/pythia',
# ]

# # Get all DPMJET source files
# dpmjet191_sources = [
#   files('src/fortran/dpmjetIII-19.1/common/dummies.f'),
#   logging_source, 
#   rangen_source
# ]
# # Files to exclude from the build
# dpmjet_exclude_files = [
#   'src/fortran/dpmjetIII-19.1/src/dpmjet/DT_RNDM.f',
#   'src/fortran/dpmjetIII-19.1/src/dpmjet/DT_RNDMST.f',
#   'src/fortran/dpmjetIII-19.1/src/dpmjet/DT_RNDMTE.f',
#   'src/fortran/dpmjetIII-19.1/src/dpmjet/PYR.f',
#   'src/fortran/dpmjetIII-19.1/src/exe/DPMJET.f',
#   'src/fortran/dpmjetIII-19.1/src/exe/PHOJET.f',
#   'src/fortran/dpmjetIII-19.1/src/exe/pho_aux.f'
# ]

# # Add all .f files from the dpmjetIII-19.1/src subdirectories, excluding specific files
# foreach f : run_command('sh', '-c', 'ls src/fortran/dpmjetIII-19.1/src/**/*.f', check:true, capture:true).stdout().split()
#   if f != '' and f not in dpmjet_exclude_files
#     dpmjet191_sources += files(f)
#   endif
# endforeach

# dpmjetIII191_srcs_arg = ','.join(dpmjet191_source_paths)
# dpmjetIII191_include_args = []
# foreach p : dpmjet_inc_dirs
#   dpmjetIII191_include_args += '-I' + (meson.project_source_root() / p)
# endforeach

# # Custom target to generate f2py wrapper for DPMJET
# dpmjet191_f2py_target = custom_target('dpmjet191_f2py',
#   output: ['_dpmjetIII191module.c', '_dpmjetIII191-f2pywrappers.f', '_dpmjetIII191.pyf'],
#   input: dpmjet191_sources,
#   command: [
#     py, files('generate_f2py.py'),
#     '_dpmjetIII191',
#     'pho_event,dt_init,dt_kkinc,dt_xsglau,pycomp,dt_initjs,dt_inucas,dt_evtout,pho_init,pho_setpar,pho_pname,pho_pmass,pho_setmdl,pho_setpdf,pho_xsect',
#     meson.global_build_root() / 'meson-logs',
#     ' '.join(dpmjetIII191_include_args),
#     '@OUTDIR@',
#     '@INPUT@',
#   ],
#   build_by_default: true,
#   capture: false
# )

# Create the DPMJET III-19.1 extension module
# py.extension_module('_dpmjetIII191',
#   sources: dpmjet191_sources + dpmjet191_f2py_target + [f2py_fortranobject_c] + [normal_source],
#   include_directories: include_directories(dpmjet_inc_dirs),
#   dependencies: [py_dep, numpy_dep],
#   c_args: ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', ],
#   fortran_args: ['-DCHROMO'],
#   install: true,
#   subdir: 'chromo/models',
# )

message('Multi-model meson build configuration complete with SIBYLL-21, EPOS-LHC, and DPMJET III-19.1')
